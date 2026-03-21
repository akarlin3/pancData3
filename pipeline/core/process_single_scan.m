function [result, b0_ref_out, gtvp_ref_out, gtvn_ref_out] = process_single_scan(ctx)
    % PROCESS_SINGLE_SCAN Process one fraction x repeat scan for a patient
    %   Handles DICOM conversion, mask saving, dose resampling, volume loading,
    %   model fitting, DIR registration, DnCNN/IVIMnet loading, and biomarker
    %   extraction. Returns a result struct with all outputs.
    %
    %   ctx — struct with all scan context (see caller for fields)
    %   b0_ref_out / gtvp_ref_out / gtvn_ref_out — updated Fx1 references
    %     (non-empty only when fi==1)
    %
    % ANALYTICAL RATIONALE — PER-SCAN PROCESSING PIPELINE
    %   This function encapsulates the complete processing chain for a single
    %   DWI acquisition (one patient, one fraction, one repeat). The pipeline
    %   proceeds through five stages:
    %
    %   Stage 1: DICOM Conversion — Raw scanner output (DICOM) is converted
    %     to NIfTI for standardized spatial metadata and multi-b-value volume
    %     organization. GTV masks (physician-drawn tumor contours) and RT dose
    %     maps are also converted/resampled to the DWI coordinate system.
    %
    %   Stage 2: Volume Loading — NIfTI DWI volumes are loaded along with
    %     b-value metadata. DnCNN-denoised and IVIMnet results (pre-computed
    %     externally) are loaded when available. B-values are sorted ascending
    %     to ensure consistent array indexing across all patients.
    %
    %   Stage 3: Model Fitting — IVIM biexponential and ADC monoexponential
    %     models are fit to the signal decay curves within the combined GTV mask.
    %     Fitting is restricted to mask voxels for computational efficiency.
    %
    %   Stage 4: Deformable Image Registration (DIR) — For fractions 2-5 and
    %     post-RT, the parameter maps and dose maps are warped to the Fx1
    %     (baseline) anatomy using displacement fields from imregdemons.
    %     This ensures voxel-level longitudinal comparisons are spatially valid.
    %
    %   Stage 5: Biomarker Extraction — Voxel-level parameter values within
    %     the GTV mask are extracted as 1D vectors for downstream statistical
    %     analysis (summary metrics, distributions, survival modeling).

    fi = ctx.fi;    % fraction index (1-5 = on-treatment, 6 = post-RT)
    rpi = ctx.rpi;  % repeat scan index (1-6, for Fx1 repeatability acquisitions)
    % Initialize reference outputs as empty; only populated when fi==1
    b0_ref_out = [];
    gtvp_ref_out = [];
    gtvn_ref_out = [];

    % Initialize result with NaN defaults so that missing/failed scans
    % propagate as NaN through downstream nanmean/nanstd computations
    % rather than causing indexing errors or silent zeros.
    result = struct();
    result.bad_dwi_list = {};       % cell array of DICOM paths that failed conversion
    result.adc_mean = nan;          % mean ADC within GTVp (mm^2/s)
    result.adc_kurtosis = nan;      % histogram kurtosis of ADC distribution
    result.d_mean = nan;            % mean IVIM D within GTVp (mm^2/s)
    result.d_kurtosis = nan;        % histogram kurtosis of D distribution
    result.d_mean_dncnn = nan;      % mean D from DnCNN-denoised fitting
    result.d_mean_ivimnet = nan;    % mean D from IVIMnet neural network fitting
    result.dmean_gtvp = nan;        % mean RT dose inside GTVp (Gy)
    result.dmean_gtvn = nan;        % mean RT dose inside GTVn (Gy)
    result.d95_gtvp = nan;          % dose to 95% of GTVp volume (Gy)
    result.d95_gtvn = nan;          % dose to 95% of GTVn volume (Gy)
    result.v50gy_gtvp = nan;        % fraction of GTVp receiving >= 50 Gy (%)
    result.v50gy_gtvn = nan;        % fraction of GTVn receiving >= 50 Gy (%)

    % Build standardised naming IDs for this scan.
    % Naming convention: {fraction}_{modality}{repeat_index}
    % e.g., "fx1_dwi1", "fx3_gtv2", "post_dose_on_dwi1"
    if fi <= ctx.n_rtdose_cols
        fx_id = ['fx' int2str(fi)];   % on-treatment fractions: "fx1"..."fx5"
    else
        fx_id = 'post';               % post-RT follow-up scan
    end
    scanID    = [fx_id '_dwi' int2str(rpi)];         % DWI NIfTI filename stem
    gtvname   = [fx_id '_gtv' int2str(rpi)];         % primary GTV mask filename stem
    gtvn_name = [fx_id '_gtvn' int2str(rpi)];        % nodal GTV mask filename stem
    dosename  = [fx_id '_dose_on_dwi' int2str(rpi)]; % dose map resampled to DWI grid

    % NIfTI output directory: all converted volumes for this patient are stored
    % in a flat 'nii' subdirectory within the patient folder.
    outloc = fullfile(ctx.basefolder, 'nii');
    if ~isfolder(outloc)
        [mk_status, mk_msg] = mkdir(outloc);
        if mk_status ~= 1
            error('process_single_scan:mkdirFailed', ...
                'Cannot create output folder %s: %s', outloc, mk_msg);
        end
    end

    bad_dwi_found = 0;  % flag: set to 1 if any conversion/loading error occurs
    bad_list = {};       % accumulates DICOM paths that failed

    % --- Convert DWI DICOMs to NIfTI using dcm2niix ---
    if ~isempty(ctx.dicomloc)
        bad_dwi_found_flag = convert_dicom(ctx.dicomloc, outloc, scanID, ctx.dcm2nii_call, fx_id);
        if bad_dwi_found_flag
            bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
            bad_dwi_found = 1;
        end
    end

    % --- Save GTVp mask as NIfTI for consistency with DWI volumes ---
    % GTV masks are stored as .mat files (from contour export) and must be
    % converted to NIfTI to match the DWI coordinate system. The rot90
    % inverse rotation (-1) converts from MATLAB's display convention back
    % to the NIfTI storage convention, ensuring spatial alignment when the
    % mask is later loaded alongside the DWI volume (which undergoes rot90
    % forward rotation at load time). safe_load_mask prevents arbitrary
    % code execution from untrusted .mat files by validating variable types.
    if ~isempty(ctx.struct_file)
        if ~exist(fullfile(outloc, [gtvname '.nii.gz']),'file')
            gtv_mask_raw = safe_load_mask(ctx.struct_file, 'Stvol3d');
            if ~isempty(gtv_mask_raw)
                niftiwrite(rot90(double(gtv_mask_raw),-1),fullfile(outloc, gtvname),'Compressed',true);
            else
                fprintf('Warning: Failed to load GTVp mask from %s\n', ctx.struct_file);
            end
        end
    end

    % --- Save GTVn (nodal) mask as NIfTI, if present ---
    % Nodal GTV masks (GTVn) are handled identically to primary (GTVp).
    % Not all patients have nodal disease involvement, so this is conditional.
    % When present, GTVn allows separate analysis of nodal vs primary tumor
    % diffusion response to RT, which may differ due to different tissue
    % composition and vascularity of lymph nodes vs pancreatic parenchyma.
    if ~isempty(ctx.struct_file_gtvn)
        if ~exist(fullfile(outloc, [gtvn_name '.nii.gz']),'file')
            gtvn_mask_raw = safe_load_mask(ctx.struct_file_gtvn, 'Stvol3d');
            if ~isempty(gtvn_mask_raw)
                niftiwrite(rot90(double(gtvn_mask_raw),-1),fullfile(outloc, gtvn_name),'Compressed',true);
            else
                fprintf('Warning: Failed to load GTVn mask from %s\n', ctx.struct_file_gtvn);
            end
        end
    end

    % --- Resample RT dose onto DWI geometry and save as NIfTI ---
    % The RT dose map is computed on the treatment planning CT grid, which
    % has different resolution, orientation, and field-of-view than the DWI
    % MRI. To enable voxel-level dose-diffusion correlation (e.g., "does
    % higher local dose at a voxel correlate with greater ADC increase?"),
    % the dose must be resampled onto the DWI geometry. This is done using
    % the b=0 DWI DICOM headers, which contain the spatial metadata (image
    % position, orientation, pixel spacing) needed to define the DWI grid.
    % b=0 images are used because they have the highest SNR and most
    % faithful spatial encoding (no diffusion-related geometric distortion).
    if ~isempty(ctx.dicomdoseloc) && ~isempty(ctx.dicomloc)
        if ~exist(fullfile(outloc, [dosename '.nii.gz']),'file')
            % Identify b=0 DICOM slices by reading the DiffusionBValue tag.
            % These slices define the spatial coordinate system for dose
            % resampling via sample_rtdose_on_image().
            dicom_files = dir(fullfile(ctx.dicomloc, '*.dcm'));
            b0list = cell(0, 1);
            b0count = 0;
            n_dicom_files = length(dicom_files);
            for bi = 1:n_dicom_files
                text_progress_bar(bi, n_dicom_files, 'Scanning DICOMs for b=0');
                data_tmp = dicominfo(fullfile(dicom_files(bi).folder, dicom_files(bi).name), 'UseDictionaryVR', true);
                if ~isfield(data_tmp, 'DiffusionBValue'), continue; end
                if data_tmp.DiffusionBValue == 0
                    b0count = b0count+1;
                    b0list{b0count,1} = fullfile(dicom_files(bi).folder, dicom_files(bi).name);
                end
            end
            rtdose_dicom = dir(fullfile(ctx.dicomdoseloc, '*.dcm'));
            if isempty(rtdose_dicom)
                warning('process_single_scan:noDose', 'No DICOM dose file found in %s', ctx.dicomdoseloc);
                rtdosefile = '';
            else
                if numel(rtdose_dicom) > 1
                    warning('process_single_scan:multipleDose', ...
                        'Multiple DICOM dose files found in %s; using first.', ctx.dicomdoseloc);
                end
                rtdosefile = fullfile(rtdose_dicom(1).folder, rtdose_dicom(1).name);
            end
            if b0count > 0 && ~isempty(rtdosefile)
                dose_sampled = sample_rtdose_on_image(b0list,rtdosefile);
                niftiwrite(rot90(dose_sampled,-1),fullfile(outloc, dosename),'Compressed',true);
            else
                if b0count == 0
                    warning('process_single_scan:noB0', 'No b=0 images found in %s. Skipping dose resampling.', ctx.dicomloc);
                end
                if isempty(rtdosefile)
                    warning('process_single_scan:noDoseFile', 'No dose DICOM found in %s. Skipping dose resampling.', ctx.dicomdoseloc);
                end
            end
        end
    end

    % --- Load NIfTI DWI volume and extract b-values ---
    % The 4D DWI volume has dimensions (x, y, z, b-value). rot90 converts
    % from NIfTI storage convention (radiological) to MATLAB display
    % convention. Voxel volume is computed in cm^3 (dimensions in mm / 10)
    % for consistency with clinical reporting of GTV volumes in cc.
    havedwi = 0;        % flag: DWI volume successfully loaded and validated
    dwi = [];           % 4D DWI volume (x, y, z, b-value)
    bvalues = [];       % b-value vector sorted ascending
    i_sort = [];        % sort indices mapping original b-value order to ascending
    dwi_vox_vol = nan;  % single voxel volume in cm^3
    dwi_dims = [];      % voxel dimensions in mm [dx, dy, dz]
    if exist(fullfile(outloc, [scanID '.nii.gz']),'file')
        dwi_info = niftiinfo(fullfile(outloc, [scanID '.nii.gz']));
        dwi = rot90(niftiread(dwi_info));
        dwi_dims = dwi_info.PixelDimensions(1:3);  % voxel dimensions in mm
        dwi_vox_vol = prod(dwi_dims*0.1);           % convert mm to cm, then volume in cm^3
        fprintf('...Loaded %s. ',fullfile(outloc, [scanID '.nii.gz']));
        havedwi = 1;

        % Load b-values from the sidecar file generated by dcm2niix.
        % B-values are sorted ascending because the IVIM segmented fit
        % assumes ordered b-values to separate low-b (perfusion-sensitive)
        % from high-b (diffusion-dominated) data points. The same sort
        % order is applied to the 4th dimension of the DWI volume so that
        % dwi(:,:,:,1) always corresponds to b=0 (the S0 reference for
        % ADC fitting) and dwi(:,:,:,end) corresponds to the highest
        % b-value (most diffusion-weighted, most sensitive to cellularity).
        bval_file = fullfile(outloc, [scanID '.bval']);
        if exist(bval_file,'file')
            fid = fopen(bval_file);
            if fid < 0
                warning('process_single_scan:bvalOpenFailed', ...
                    'Cannot open bval file %s for %s. Skipping.', bval_file, scanID);
                havedwi = 0;
                bvalues = [];
            end
            if havedwi
                tline = fgetl(fid);
                fclose(fid);
                if ~ischar(tline)
                    warning('process_single_scan:bvalEmpty', ...
                        'Empty bval file %s for %s. Skipping.', bval_file, scanID);
                    havedwi = 0;
                    bvalues = [];
                else
                    bvalues = sscanf(tline, '%f');
                end
            end
            % Validate volume count BEFORE sorting to prevent indexing
            % errors when the bval file has more entries than the NIfTI
            % has volumes (e.g., truncated acquisition, protocol mismatch).
            if size(dwi, 4) ~= numel(bvalues)
                warning('process_single_scan:bvalVolumeMismatch', ...
                    'DWI has %d volumes but bval file has %d entries for %s. Skipping.', ...
                    size(dwi, 4), numel(bvalues), scanID);
                havedwi = 0;
                bvalues = [];
                if bad_dwi_found==0
                    bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
                    bad_dwi_found = 1;
                end
            else
                [~,i_sort] = sort(bvalues,'ascend');
                bvalues = bvalues(i_sort);
                dwi = double(dwi(:,:,:,i_sort));
                fprintf('loaded bvalues\n');
            end
        else
            fprintf('bvalue file not found!\n');
            havedwi = 0;
            if bad_dwi_found==0
                bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
                bad_dwi_found = 1;
            end
        end

        % Validate minimum volume count for DWI processing.
        if isempty(bvalues) && size(dwi, 4) < 2
            fprintf('DWI does not have enough volumes: %s — skipping\n', mat2str(size(dwi)));
            havedwi = 0;
            if bad_dwi_found==0
                bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
            end
        end
    end

    % --- Load DnCNN-denoised DWI (deep learning denoising) ---
    % DnCNN (Denoising Convolutional Neural Network) reduces Rician noise
    % in DWI images before conventional IVIM model fitting. Pancreatic DWI
    % has inherently low SNR due to: (1) small organ size, (2) respiratory
    % and peristaltic motion, (3) signal attenuation at high b-values.
    % Denoising improves parameter estimation stability, especially for
    % the highly variable D* parameter. The denoised volumes are either
    % pre-computed (cached as NIfTI) or computed on-the-fly as a fallback.
    havedenoised = 0;
    dwi_dncnn = [];
    if havedwi==1
        dncnnid = [scanID '_dncnn.nii.gz'];
        dncnn_file = fullfile(ctx.basefolder, 'dncnn', dncnnid);
        if exist(dncnn_file,'file')
            dncnn_info = niftiinfo(dncnn_file);
            dwi_dncnn = rot90(niftiread(dncnn_info));
            % Validate that the DnCNN cache has the same number of b-value
            % volumes as the raw DWI.  A mismatch (e.g., cache generated
            % from a different protocol or truncated file) would cause
            % i_sort indexing to exceed the 4th dimension.
            if size(dwi_dncnn, 4) ~= numel(i_sort)
                warning('process_single_scan:dncnnSizeMismatch', ...
                    'DnCNN cache has %d volumes but DWI has %d b-values. Skipping cached DnCNN for %s.', ...
                    size(dwi_dncnn, 4), numel(i_sort), scanID);
                dwi_dncnn = [];
            else
                % Normalise denoised signal to [0,1] for IVIM fitting.
                % mat2gray maps the volume's [min,max] to [0,1], which is
                % required because the cached DnCNN output may have different
                % intensity scaling than the raw DWI.
                dwi_dncnn = double(mat2gray(dwi_dncnn(:,:,:,i_sort)));
                havedenoised=1;
            end
        else
            % Load GTVp/GTVn masks for the fallback (may already exist on disk)
            gtv_mask_for_dncnn = [];
            gtvn_mask_for_dncnn = [];
            gtvp_filepath = fullfile(outloc, [gtvname '.nii.gz']);
            if exist(gtvp_filepath, 'file')
                gtv_mask_for_dncnn = rot90(niftiread(niftiinfo(gtvp_filepath)));
            end
            gtvn_filepath = fullfile(outloc, [gtvn_name '.nii.gz']);
            if exist(gtvn_filepath, 'file')
                gtvn_mask_for_dncnn = rot90(niftiread(niftiinfo(gtvn_filepath)));
            end
            use_gpu_dncnn = isfield(ctx, 'use_gpu') && ctx.use_gpu;
            [dwi_dncnn, havedenoised] = compute_dncnn_fallback(dwi, i_sort, gtv_mask_for_dncnn, gtvn_mask_for_dncnn, use_gpu_dncnn);
        end
    end

    % --- Load IVIMnet deep-learning fit results (pre-computed) ---
    % IVIMnet is a neural network trained to estimate IVIM parameters (D, f, D*)
    % directly from multi-b-value signal curves, bypassing the ill-conditioned
    % nonlinear fitting that makes conventional D* estimation unreliable.
    % IVIMnet results are pre-computed externally (Python/TensorFlow) and
    % stored as .mat files containing 3D parameter maps. Unlike DnCNN (which
    % denoises the raw DWI before conventional fitting), IVIMnet replaces
    % the fitting step entirely, providing an alternative estimation approach
    % for methodological comparison.
    haveivimnet = 0;
    D_ivimnet = []; f_ivimnet = []; Dstar_ivimnet = [];
    if havedwi==1
        ivimid = [scanID '_ivimnet.mat'];
        ivimnet_file = fullfile(ctx.basefolder, 'ivimnet', ivimid);
        if exist(ivimnet_file,'file')
            tmp = load(ivimnet_file);
            % Apply rot90 at load time, consistent with how DWI and DnCNN
            % volumes are rotated from NIfTI to MATLAB orientation (line 112/154).
            % IVIMnet .mat files store maps in NIfTI convention.
            D_ivimnet = rot90(tmp.D_ivimnet);
            f_ivimnet = rot90(tmp.f_ivimnet);
            Dstar_ivimnet = rot90(tmp.Dstar_ivimnet);
            haveivimnet=1;
        end
    end

    if havedwi
        dwi_size = size(dwi);
    else
        dwi_size = [0 0 0 0];
    end

    % --- Load GTVp mask and validate spatial dimensions ---
    havegtvp = 0; gtv_mask = [];   % primary tumor (GTVp) mask
    havegtvn = 0; gtvn_mask = [];  % nodal disease (GTVn) mask
    if havedwi
        gtvp_filepath = fullfile(outloc, [gtvname '.nii.gz']);
        [havegtvp, gtv_mask] = load_nifti_mask(gtvp_filepath, dwi_size, '', 'gtvp');

        % --- Load GTVn (nodal) mask and validate dimensions ---
        gtvn_filepath = fullfile(outloc, [gtvn_name '.nii.gz']);
        [havegtvn, gtvn_mask] = load_nifti_mask(gtvn_filepath, dwi_size, '*NODAL* ', 'gtvn');
    end

    % --- Load resampled RT dose map ---
    havedose = 0;
    dose_map = [];
    if exist(fullfile(outloc, [dosename '.nii.gz']),'file')
        dose_info = niftiinfo(fullfile(outloc, [dosename '.nii.gz']));
        dose_map = rot90(niftiread(dose_info));
        fprintf('...Loaded %s\n',fullfile(outloc, [dosename '.nii.gz']));
        havedose = 1;
    end

    % --- Fit ADC and IVIM models ---
    % Model fitting is performed on the UNION of GTVp and GTVn masks.
    % Using the union ensures that both primary tumor and nodal disease
    % voxels receive fitted parameters in a single pass, avoiding redundant
    % computation. The combined mask is passed to fit_models() which
    % restricts fitting to only masked voxels for efficiency (~100-2000
    % tumor voxels out of ~1.3M total volume voxels).
    %
    % Models are fit independently on:
    %   1. Standard (raw) DWI — baseline comparison
    %   2. DnCNN-denoised DWI — noise-reduced estimation
    % IVIMnet parameters are pre-computed externally and loaded above.
    % Pre-initialize parameter maps to empty: populated by fit_models() below
    d_map = []; f_map = []; dstar_map = []; adc_map = [];               % Standard DWI
    d_map_dncnn = []; f_map_dncnn = []; dstar_map_dncnn = []; adc_map_dncnn = [];  % DnCNN-denoised
    if havedwi && (havegtvn || havegtvp)
        % Build the union mask of GTVp and GTVn for model fitting.
        % logical() + logical() = OR operation, ensuring all tumor voxels
        % (primary and nodal) are included in a single fitting pass.
        mask_ivim = false(size(dwi,1),size(dwi,2),size(dwi,3));
        if havegtvp, mask_ivim = logical(mask_ivim + logical(gtv_mask)); end
        if havegtvn, mask_ivim = logical(mask_ivim + logical(gtvn_mask)); end

        % --- Motion artifact detection (optional) ---
        % When enabled, flag DWI volumes with excessive motion corruption
        % before model fitting. Flagged volumes can be excluded to improve
        % parameter estimation accuracy.
        if isfield(ctx, 'exclude_motion_volumes') && ctx.exclude_motion_volumes
            try
                motion = detect_motion_artifacts(dwi, bvalues, mask_ivim);
                if any(motion.flagged)
                    n_flagged = sum(motion.flagged);
                    fprintf('    ⚠️  Motion: %d/%d volumes flagged, excluding from fit\n', ...
                        n_flagged, numel(motion.flagged));
                    % Remove flagged volumes from DWI and b-values
                    keep_vols = ~motion.flagged;
                    dwi = dwi(:,:,:,keep_vols);
                    bvalues = bvalues(keep_vols);
                    if havedenoised == 1
                        dwi_dncnn = dwi_dncnn(:,:,:,keep_vols);
                    end
                end
            catch me_motion
                fprintf('    💡 Motion detection skipped: %s\n', me_motion.message);
            end
        end

        % opts.bthr: b-value threshold (typically 100 s/mm^2) for IVIM
        % segmented fitting. b < bthr captures perfusion, b >= bthr captures diffusion.
        opts = [];
        opts.bthr = ctx.ivim_bthr;
        % Pass GPU flag to fit_models for GPU-accelerated ADC WLS fitting.
        if isfield(ctx, 'use_gpu')
            opts.use_gpu = ctx.use_gpu;
        end

        [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts);
        if havedenoised==1
            % Fit the same models on denoised data. This produces a matched
            % set of parameter maps that can be directly compared to standard
            % maps to quantify the effect of denoising on parameter estimates.
            [d_map_dncnn, f_map_dncnn, dstar_map_dncnn, adc_map_dncnn] = fit_models(dwi_dncnn, bvalues, mask_ivim, opts);
        end
    end

    % --- Clear DWI volume variables to reduce memory usage ---
    % Model fitting is complete; clear large DWI volumes immediately to
    % reduce peak memory usage before visualization/metric computation
    clear dwi dwi_dncnn;

    % --- Deformable Image Registration (DIR) ---
    % DIR warps later-fraction (Fx2-Fx5, post) images to the Fx1 baseline
    % coordinate system using the b=0 images as anatomical references.
    % This is critical because:
    %   1. Patient setup varies between fractions (different table position,
    %      organ filling, respiratory phase)
    %   2. The tumor may deform, shrink, or shift during the treatment course
    %   3. Voxel-level longitudinal comparisons (e.g., "did ADC at voxel X
    %      increase between Fx1 and Fx3?") require spatial correspondence
    %
    % The displacement field (D_forward) maps each voxel in Fx1 space to its
    % corresponding location in the current fraction's space. Warping
    % parameter maps by -D_forward brings them INTO Fx1 space. The GTV mask
    % is also warped to define the tumor boundary in later fractions, which
    % may differ from the original Fx1 contour due to tumor deformation.
    gtv_mask_for_dvh = gtv_mask;
    dose_map_dvh     = dose_map;
    D_forward_cur    = [];

    if havedwi && havegtvp
        % Use dwi_info to recreate b0_current from disk to avoid keeping dwi in memory
        b0_current = rot90(niftiread(dwi_info));
        b0_current = double(b0_current(:,:,:,i_sort));
        b0_current = b0_current(:,:,:,1);  % b=0 volume: highest-SNR anatomical reference
        
        if fi == 1
            % At Fx1 (baseline), save the b=0 volume and GTV masks as
            % reference images for subsequent fraction registrations.
            % No warping is needed at baseline — it IS the reference frame.
            b0_ref_out = b0_current;
            gtvp_ref_out = gtv_mask;
            if havegtvn
                gtvn_ref_out = gtvn_mask;
            end
        elseif fi > 1 && ~isempty(ctx.b0_fx1_ref) && ~isempty(ctx.gtv_mask_fx1_ref)
            dir_cache_file = fullfile(ctx.dataloc, ctx.id_j, 'nii', ...
                sprintf('dir_field_rpi%d_fx%d.mat', rpi, fi));
            if exist(dir_cache_file, 'file')
                tmp_dir = load(dir_cache_file, 'gtv_mask_warped', 'D_forward', 'ref3d');
                gtv_mask_for_dvh = tmp_dir.gtv_mask_warped;
                if isfield(tmp_dir, 'D_forward'),  D_forward_cur = tmp_dir.D_forward;  end
                fprintf('  [DIR] Loaded cached warped mask + D_forward for Fx%d rpi%d\n', fi, rpi);
            else
                fprintf('  [DIR] Running imregdemons for Fx%d rpi%d...\n', fi, rpi);
                [gtv_mask