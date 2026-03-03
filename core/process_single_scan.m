function [result, b0_ref_out, gtvp_ref_out, gtvn_ref_out] = process_single_scan(ctx)
    % PROCESS_SINGLE_SCAN Process one fraction x repeat scan for a patient
    %   Handles DICOM conversion, mask saving, dose resampling, volume loading,
    %   model fitting, DIR registration, DnCNN/IVIMnet loading, and biomarker
    %   extraction. Returns a result struct with all outputs.
    %
    %   ctx — struct with all scan context (see caller for fields)
    %   b0_ref_out / gtvp_ref_out / gtvn_ref_out — updated Fx1 references
    %     (non-empty only when fi==1)

    fi = ctx.fi;
    rpi = ctx.rpi;
    b0_ref_out = [];
    gtvp_ref_out = [];
    gtvn_ref_out = [];

    % Initialize result with NaN defaults
    result = struct();
    result.bad_dwi_list = {};
    result.adc_mean = nan;
    result.adc_kurtosis = nan;
    result.d_mean = nan;
    result.d_kurtosis = nan;
    result.d_mean_dncnn = nan;
    result.d_mean_ivimnet = nan;
    result.dmean_gtvp = nan;
    result.dmean_gtvn = nan;
    result.d95_gtvp = nan;
    result.d95_gtvn = nan;
    result.v50gy_gtvp = nan;
    result.v50gy_gtvn = nan;

    % Build standardised naming IDs for this scan
    if fi <= ctx.n_rtdose_cols
        fx_id = ['fx' int2str(fi)];
    else
        fx_id = 'post';
    end
    scanID    = [fx_id '_dwi' int2str(rpi)];
    gtvname   = [fx_id '_gtv' int2str(rpi)];
    gtvn_name = [fx_id '_gtvn' int2str(rpi)];
    dosename  = [fx_id '_dose_on_dwi' int2str(rpi)];

    outloc = fullfile(ctx.basefolder, 'nii');
    if ~isfolder(outloc), mkdir(outloc); end

    bad_dwi_found = 0;
    bad_list = {};

    % --- Convert DWI DICOMs to NIfTI using dcm2niix ---
    if ~isempty(ctx.dicomloc)
        bad_dwi_found_flag = convert_dicom(ctx.dicomloc, outloc, scanID, ctx.dcm2nii_call, fx_id);
        if bad_dwi_found_flag
            bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
            bad_dwi_found = 1;
        end
    end

    % --- Save GTVp mask as NIfTI for consistency with DWI volumes ---
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
    if ~isempty(ctx.dicomdoseloc) && ~isempty(ctx.dicomloc)
        if ~exist(fullfile(outloc, [dosename '.nii.gz']),'file')
            dicom_files = dir(fullfile(ctx.dicomloc, '*.dcm'));
            b0list = cell(1);
            b0count = 0;
            for bi = 1:length(dicom_files)
                data_tmp = dicominfo(fullfile(dicom_files(bi).folder, dicom_files(bi).name), 'UseDictionaryVR', true);
                if data_tmp.DiffusionBValue == 0
                    b0count = b0count+1;
                    b0list{b0count,1} = fullfile(dicom_files(bi).folder, dicom_files(bi).name);
                end
            end
            rtdose_dicom = dir(fullfile(ctx.dicomdoseloc, '*.dcm'));
            rtdosefile = fullfile(rtdose_dicom.folder, rtdose_dicom.name);
            dose_sampled = sample_rtdose_on_image(b0list,rtdosefile);
            niftiwrite(rot90(dose_sampled,-1),fullfile(outloc, dosename),'Compressed',true);
        end
    end

    % --- Load NIfTI DWI volume and extract b-values ---
    havedwi = 0;
    dwi = [];
    bvalues = [];
    i_sort = [];
    dwi_vox_vol = nan;
    dwi_dims = [];
    if exist(fullfile(outloc, [scanID '.nii.gz']),'file')
        dwi_info = niftiinfo(fullfile(outloc, [scanID '.nii.gz']));
        dwi = rot90(niftiread(dwi_info));
        dwi_dims = dwi_info.PixelDimensions(1:3);
        dwi_vox_vol = prod(dwi_dims*0.1);
        fprintf('...Loaded %s. ',fullfile(outloc, [scanID '.nii.gz']));
        havedwi = 1;

        bval_file = fullfile(outloc, [scanID '.bval']);
        if exist(bval_file,'file')
            fid = fopen(bval_file);
            tline = fgetl(fid);
            fclose(fid);
            bvalues = sscanf(tline, '%f');
            [~,i_sort] = sort(bvalues,'ascend');
            bvalues = bvalues(i_sort);
            dwi = double(dwi(:,:,:,i_sort));
            fprintf('loaded bvalues\n');
        else
            fprintf('bvalue file not found!\n');
            havedwi = 0;
            if bad_dwi_found==0
                bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
                bad_dwi_found = 1;
            end
        end

        if size(dwi,4)~=4
            fprintf('DWI does not have expected dimensions found: %s skipping\n',mat2str(size(dwi)))
            havedwi = 0;
            if bad_dwi_found==0
                bad_list{end+1} = ctx.dicomloc; %#ok<AGROW>
            end
        end
    end

    % --- Load DnCNN-denoised DWI (deep learning denoising) ---
    havedenoised = 0;
    dwi_dncnn = [];
    if havedwi==1
        dncnnid = [scanID '_dncnn.nii.gz'];
        dncnn_file = fullfile(ctx.basefolder, 'dncnn', dncnnid);
        if exist(dncnn_file,'file')
            dncnn_info = niftiinfo(dncnn_file);
            dwi_dncnn = rot90(niftiread(dncnn_info));
            dwi_dncnn = double(mat2gray(dwi_dncnn(:,:,:,i_sort)));
            havedenoised=1;
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
            [dwi_dncnn, havedenoised] = compute_dncnn_fallback(dwi, i_sort, gtv_mask_for_dncnn, gtvn_mask_for_dncnn);
        end
    end

    % --- Load IVIMnet deep-learning fit results (pre-computed) ---
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
    havegtvp = 0; gtv_mask = [];
    havegtvn = 0; gtvn_mask = [];
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
    d_map = []; f_map = []; dstar_map = []; adc_map = [];
    d_map_dncnn = []; f_map_dncnn = []; dstar_map_dncnn = []; adc_map_dncnn = [];
    if havedwi && (havegtvn || havegtvp)
        mask_ivim = false(size(dwi,1),size(dwi,2),size(dwi,3));
        if havegtvp, mask_ivim = logical(mask_ivim + logical(gtv_mask)); end
        if havegtvn, mask_ivim = logical(mask_ivim + logical(gtvn_mask)); end

        opts = [];
        opts.bthr = ctx.ivim_bthr;

        [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts);
        if havedenoised==1
            [d_map_dncnn, f_map_dncnn, dstar_map_dncnn, adc_map_dncnn] = fit_models(dwi_dncnn, bvalues, mask_ivim, opts);
        end
    end

    % --- Deformable Image Registration (DIR) ---
    gtv_mask_for_dvh = gtv_mask;
    dose_map_dvh     = dose_map;
    D_forward_cur    = [];

    if havedwi && havegtvp
        b0_current = dwi(:,:,:,1);
        if fi == 1
            % Return Fx1 references to caller
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
                [gtv_mask_warped, D_forward_cur, ref3d_cur] = ...
                    apply_dir_mask_propagation(ctx.b0_fx1_ref, b0_current, ctx.gtv_mask_fx1_ref);
                if ~isempty(gtv_mask_warped)
                    gtv_mask_for_dvh = gtv_mask_warped;
                    parsave_dir_cache(dir_cache_file, gtv_mask_warped, D_forward_cur, ref3d_cur);
                    fprintf('  [DIR] Done. Warped mask + D_forward saved.\n');
                else
                    fprintf('  [DIR] Registration failed for Fx%d rpi%d. Falling back to rigid dose/mask.\n', fi, rpi);
                end
            end
        end
    end

    % --- Warp native-space parameter maps to baseline geometry -----------
    % All DWI types are warped consistently so longitudinal comparisons
    % use the same spatial reference frame.
    can_warp = (fi > 1) && ~isempty(D_forward_cur) && ~isempty(ctx.b0_fx1_ref);
    if can_warp
        % Standard maps
        adc_map = imwarp(adc_map, -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
        d_map   = imwarp(d_map,   -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
        f_map   = imwarp(f_map,   -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
        dstar_map = imwarp(dstar_map, -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);

        % DnCNN maps
        if havedenoised
            d_map_dncnn     = imwarp(d_map_dncnn,     -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
            f_map_dncnn     = imwarp(f_map_dncnn,     -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
            dstar_map_dncnn = imwarp(dstar_map_dncnn, -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
            adc_map_dncnn   = imwarp(adc_map_dncnn,   -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
        end

        % IVIMnet maps
        if haveivimnet
            D_ivimnet     = imwarp(D_ivimnet,     -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
            f_ivimnet     = imwarp(f_ivimnet,     -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
            Dstar_ivimnet = imwarp(Dstar_ivimnet, -D_forward_cur, 'Interp', 'linear', 'FillValues', nan);
        end

        % Warp dose map to baseline geometry so DVH metrics are in the same
        % reference frame as warped parameter maps.  Without this, dose stays
        % in native Fx2+ space while parameters are in Fx1 space, creating a
        % spatial mismatch in dose-parameter correlations.
        if havedose
            dose_map = imwarp(dose_map, -D_forward_cur, 'Interp', 'linear', 'FillValues', 0);
        end

        fprintf('  [DIR] Warped all parameter maps and dose to baseline geometry.\n');
    end

    % Build maps struct for extract_biomarkers
    maps = struct('adc_map', adc_map, 'd_map', d_map, 'f_map', f_map, 'dstar_map', dstar_map);
    meta = struct('id', ctx.id_j, 'mrn', ctx.mrn_j, 'lf', ctx.pat_lf, ...
        'immuno', ctx.pat_immuno, 'fi', fi, 'rpi', rpi, 'vox_vol', dwi_vox_vol);

    % DnCNN maps struct
    dncnn_maps = struct();
    if havedenoised
        dncnn_maps.d_map_dncnn     = d_map_dncnn;
        dncnn_maps.f_map_dncnn     = f_map_dncnn;
        dncnn_maps.dstar_map_dncnn = dstar_map_dncnn;
        dncnn_maps.adc_map_dncnn   = adc_map_dncnn;
    end

    % IVIMnet maps struct (already rotated at load time, lines 182-184)
    ivimnet_maps = struct();
    if haveivimnet
        ivimnet_maps.D_ivimnet     = D_ivimnet;
        ivimnet_maps.f_ivimnet     = f_ivimnet;
        ivimnet_maps.Dstar_ivimnet = Dstar_ivimnet;
    end

    % Determine the appropriate mask for biomarker extraction.
    % When parameter maps are warped to Fx1 (baseline) geometry (can_warp),
    % ALL DWI types must use the Fx1 reference mask for consistency.
    % Previously DnCNN used the Fx1 mask while standard maps used the
    % native Fx2+ mask, causing a spatial mismatch that invalidated
    % cross-DWI-type longitudinal comparisons.
    biomarker_mask_p = gtv_mask;
    biomarker_mask_n = gtvn_mask;
    if can_warp && ~isempty(ctx.gtv_mask_fx1_ref)
        biomarker_mask_p = ctx.gtv_mask_fx1_ref;
        if ~isempty(ctx.gtvn_mask_fx1_ref)
            biomarker_mask_n = ctx.gtvn_mask_fx1_ref;
        end
    end

    % DnCNN masks now use the same mask as all other DWI types
    dncnn_mask_p = biomarker_mask_p;
    dncnn_mask_n = biomarker_mask_n;

    % --- Extract biomarkers for GTVp ---
    % Initialize with empty struct matching init_scan_structs fields
    [empty_p, empty_n] = init_scan_structs(1, 1);
    result.data_gtvp = empty_p;
    result.data_gtvn = empty_n;

    if havegtvp
        result.data_gtvp = extract_biomarkers(biomarker_mask_p, maps, meta, dncnn_maps, dncnn_mask_p, ivimnet_maps);

        result.adc_mean = nanmean(adc_map(biomarker_mask_p==1));
        % NOTE: Histogram kurtosis of trace-average ADC — NOT valid DKI.
        result.adc_kurtosis = kurtosis(adc_map(biomarker_mask_p==1));
        result.d_mean = nanmean(d_map(biomarker_mask_p==1));
        % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
        result.d_kurtosis = kurtosis(d_map(biomarker_mask_p==1));

        if havedenoised
            result.d_mean_dncnn = nanmean(d_map_dncnn(dncnn_mask_p==1));
        end
        if haveivimnet
            result.d_mean_ivimnet = nanmean(D_ivimnet(biomarker_mask_p==1));
        end
    end

    % --- DVH for GTVp ---
    % Always use the same mask for dose as for biomarkers so that
    % dose_vector and adc_vector have identical lengths.  When maps are
    % warped to Fx1 the dose map is already in that frame; otherwise the
    % native dose map is resampled onto the native biomarker mask.
    dvh_mask_p = biomarker_mask_p;
    dvh_dose   = dose_map_dvh;
    if can_warp && havedose
        dvh_dose   = dose_map;          % dose already warped to Fx1 above
    end
    if havedose && havegtvp
        result.dmean_gtvp = nanmean(dvh_dose(dvh_mask_p==1));
        [dvhparams, dvh_values] = dvh(dvh_dose, dvh_mask_p, dwi_dims, 2000, 'Dperc',95,'Vperc',50,'Normalize',true);
        result.d95_gtvp = dvhparams.("D95% (Gy)");
        result.v50gy_gtvp = dvhparams.("V50Gy (%)");

        result.data_gtvp.dose_vector = dvh_dose(dvh_mask_p==1);
        result.data_gtvp.dvh = dvh_values;
        result.data_gtvp.d95 = dvhparams.("D95% (Gy)");
        result.data_gtvp.v50gy = dvhparams.("V50Gy (%)");
    end

    % --- Extract biomarkers for GTVn ---
    if havegtvn
        % GTVn uses the same IVIMnet maps as GTVp (already rotated once at
        % lines 296-298).  The legacy code incorrectly applied rot90 a
        % second time, causing a 180-deg orientation mismatch.
        ivimnet_maps_n = ivimnet_maps;
        result.data_gtvn = extract_biomarkers(biomarker_mask_n, maps, meta, dncnn_maps, dncnn_mask_n, ivimnet_maps_n);
    end

    % --- DVH for GTVn ---
    % Always use the same mask for dose as for biomarkers (mirrors GTVp).
    dvh_mask_n = biomarker_mask_n;
    dvh_dose_n = dose_map;
    if havedose && havegtvn
        result.dmean_gtvn = nanmean(dvh_dose_n(dvh_mask_n==1));
        [dvhparams, dvh_values] = dvh(dvh_dose_n, dvh_mask_n, dwi_dims, 2000, 'Dperc',95,'Vperc',50,'Normalize',true);
        result.d95_gtvn = dvhparams.("D95% (Gy)");
        result.v50gy_gtvn = dvhparams.("V50Gy (%)");

        result.data_gtvn.dose_vector = dvh_dose_n(dvh_mask_n==1);
        result.data_gtvn.dvh = dvh_values;
        result.data_gtvn.d95 = dvhparams.("D95% (Gy)");
        result.data_gtvn.v50gy = dvhparams.("V50Gy (%)");
    end

    result.bad_dwi_list = bad_list;
end

% =========================================================================
%  Local helper functions
% =========================================================================

function [have_mask, mask_data] = load_nifti_mask(filepath, dwi_size, message_prefix, mask_name)
    % LOAD_NIFTI_MASK Helper function to load a NIfTI mask and validate dimensions
    have_mask = 0;
    mask_data = [];

    if exist(filepath, 'file')
        info = niftiinfo(filepath);
        mask_data = rot90(niftiread(info));
        fprintf('...%sLoaded %s\n', message_prefix, filepath);
        have_mask = 1;

        mask_size = size(mask_data);
        if sum(mask_size ~= dwi_size(1:3)) > 0
            have_mask = 0;
            mask_data = [];
            fprintf('size mismatch. excluding %s\n', mask_name);
        end
    end
end

function bio = extract_biomarkers(mask, maps, meta, dncnn_maps, dncnn_mask, ivimnet_maps)
    % EXTRACT_BIOMARKERS Extract voxel-level biomarkers within a GTV mask
    %   Returns a struct with all fields matching init_scan_structs layout.

    % Start from template to ensure all fields exist
    [tmp, ~] = init_scan_structs(1, 1);
    bio = tmp;

    mask_idx = (mask == 1);

    bio.adc_vector = maps.adc_map(mask_idx);
    bio.d_vector   = maps.d_map(mask_idx);
    bio.f_vector   = maps.f_map(mask_idx);
    bio.dstar_vector = maps.dstar_map(mask_idx);
    bio.ID = meta.id;
    bio.MRN = meta.mrn;
    bio.LF = meta.lf;
    bio.Immuno = meta.immuno;
    bio.Fraction = meta.fi;
    bio.Repeatability_index = meta.rpi;
    bio.vox_vol = meta.vox_vol;

    % DnCNN-denoised vectors
    if isfield(dncnn_maps, 'd_map_dncnn') && ~isempty(dncnn_maps.d_map_dncnn)
        dm = (dncnn_mask == 1);
        bio.d_vector_dncnn     = dncnn_maps.d_map_dncnn(dm);
        bio.f_vector_dncnn     = dncnn_maps.f_map_dncnn(dm);
        bio.dstar_vector_dncnn = dncnn_maps.dstar_map_dncnn(dm);
        bio.adc_vector_dncnn   = dncnn_maps.adc_map_dncnn(dm);
    end

    % IVIMnet vectors
    if isfield(ivimnet_maps, 'D_ivimnet') && ~isempty(ivimnet_maps.D_ivimnet)
        bio.d_vector_ivimnet     = ivimnet_maps.D_ivimnet(mask_idx);
        bio.f_vector_ivimnet     = ivimnet_maps.f_ivimnet(mask_idx);
        bio.dstar_vector_ivimnet = ivimnet_maps.Dstar_ivimnet(mask_idx);
    end
end

function [dwi_dncnn, havedenoised] = compute_dncnn_fallback(dwi, i_sort, gtv_mask, gtvn_mask)
    % COMPUTE_DNCNN_FALLBACK On-the-fly DnCNN denoising when cache is missing
    havedenoised = 0;
    dwi_dncnn = [];
    fprintf('  [DnCNN] Cache missing. Executing deep learning denoising on CPU...\n');
    try
        loaded_model = load(fullfile(fileparts(mfilename('fullpath')), '..', 'dependencies', 'dncnn_model.mat'), 'net');
        dncnn_net = loaded_model.net;

        dwi_cpu = single(dwi);

        if ~isempty(gtv_mask)
            mask_cpu = single(gtv_mask);
        elseif ~isempty(gtvn_mask)
            mask_cpu = single(gtvn_mask);
        else
            mask_cpu = ones(size(dwi, 1), size(dwi, 2), size(dwi, 3), 'single');
        end

        dwi_dncnn_cpu = zeros(size(dwi_cpu), 'single');
        for b_idx = 1:size(dwi_cpu, 4)
            dwi_dncnn_cpu(:,:,:,b_idx) = apply_dncnn_symmetric(dwi_cpu(:,:,:,b_idx), mask_cpu, dncnn_net, 15);
        end

        dwi_dncnn = double(mat2gray(double(dwi_dncnn_cpu(:,:,:,i_sort))));
        havedenoised = 1;
        fprintf('  [DnCNN] Deep learning denoising completed.\n');
    catch CPU_ME
        fprintf('  [DnCNN] CPU Computation failed: %s\n', CPU_ME.message);
        fprintf('  [DnCNN] Proceeding without denoised data.\n');
    end
end
