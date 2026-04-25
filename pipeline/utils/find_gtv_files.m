function [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, pat_name, opts)
% FIND_GTV_FILES Locate GTVp and (optionally) GTVn mask files for a scan.
%
%   [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, pat_name)
%   [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, pat_name, opts)
%
%   Inputs:
%       fxfolder - Absolute path to the fraction folder containing mask files
%       dwii     - Repeatability index (1-based scan repeat within the fraction)
%       pat_name - Patient identifier string; the token 'two' in the name
%                  signals that this patient has both a primary pancreatic
%                  tumor (GTVp) and lymph node metastases (GTVn)
%       opts     - (Optional) struct with discovery flags:
%                  .process_gtvn (logical, default true) — when false, the
%                     dual-GTV branch is bypassed even for "twoGTVs"
%                     patients: the function takes the single-GTV path,
%                     gtvn_path is always '', and the broad '*GTV*'
%                     fallback is allowed for primary discovery.
%
%   Outputs:
%       gtvp_path - Full path to the primary GTV mask file ('' if not found)
%       gtvn_path - Full path to the nodal GTV mask file ('' if not found;
%                   always '' for single-GTV patients)
%
% --- Analytical Rationale ---
% Pancreatic cancer patients undergoing adaptive radiotherapy may have:
%   (a) A single primary tumor (GTVp only) -- the majority of cases
%   (b) A primary tumor plus involved lymph nodes (GTVp + GTVn)
%
% The distinction matters because diffusion parameters (ADC, D, f, D*) can
% differ substantially between the primary pancreatic mass and nodal disease
% due to differences in tissue cellularity, perfusion, and fibrosis. Mixing
% voxels from both volumes into a single analysis would introduce biological
% heterogeneity that obscures treatment response signals. Therefore, GTVp and
% GTVn are tracked separately when both are present.
%
% The 'two' token in the patient name is a clinical data convention at MSK
% that flags patients with dual target volumes. When only GTVp is present,
% a broader wildcard '*GTV*' is used as a final fallback; this fallback is
% deliberately excluded from the dual-GTV branch to prevent the broad
% pattern from accidentally matching a nodal mask as the primary tumor.
%
% --- Candidate-Folder Search ---
% Some sites place GTV .mat files directly in the Fx folder, while others
% nest them in a dedicated subfolder such as GTV/, GTVtimepoint<N>/,
% GTV_tp<N>/, GTVtp<N>/, GTVfx<N>/, or GTVfraction<N>/. To handle both,
% this function builds a list of candidate folders to search and delegates
% pattern matching to discover_gtv_file for each. The first non-empty
% result wins, searching the Fx folder first (highest priority).

    gtvp_path = '';
    gtvn_path = '';

    if nargin < 4 || isempty(opts) || ~isstruct(opts)
        opts = struct();
    end
    if ~isfield(opts, 'process_gtvn')
        opts.process_gtvn = true;
    end
    process_gtvn = logical(opts.process_gtvn);

    candidate_folders = local_candidate_gtv_folders(fxfolder, dwii);

    if ~contains(pat_name, 'two') || ~process_gtvn
        % --- Single-GTV Patient (Primary Tumor Only) ---
        % Search patterns are ordered from most specific to least specific:
        %   1. *GTV_MR   - Standard MSK MR-contoured GTV naming
        %   2. *GTVp     - Explicit primary GTV label
        %   3. *GTV_panc* - Pancreas-specific GTV label (some sites)
        %   4. *GTV*     - Broad fallback; safe here because there is no
        %                  nodal GTV to confuse it with
        % The date portion of filenames is deliberately avoided in patterns
        % because DICOM StudyDate strings embedded in filenames would create
        % fragile, patient-specific search patterns.
        gtvp_patterns = {'*GTV_MR', '*GTVp', '*GTV_panc*', '*GTV*'};
        for k = 1:numel(candidate_folders)
            gtvp_path = discover_gtv_file(candidate_folders{k}, gtvp_patterns, dwii);
            if ~isempty(gtvp_path); break; end
        end
    else
        % --- Dual-GTV Patient (Primary Tumor + Nodal Metastasis) ---
        % The broad '*GTV*' fallback is intentionally omitted from GTVp
        % patterns here. If it were included, it could match a GTVn file
        % (e.g., 'GTVn1.mat'), incorrectly labeling nodal tissue as the
        % primary tumor and corrupting downstream dose-response analysis.
        gtvp_patterns = {'*GTV_MR', '*GTVp', '*GTV_panc*'};
        for k = 1:numel(candidate_folders)
            gtvp_path = discover_gtv_file(candidate_folders{k}, gtvp_patterns, dwii);
            if ~isempty(gtvp_path); break; end
        end

        % Nodal GTV naming conventions vary by contouring physician:
        %   - GTV_LN  = "Gross Tumor Volume - Lymph Node"
        %   - GTVn    = Standard ICRU-style nodal GTV label
        %   - GTV_node = Descriptive alternative used at some sites
        % These patterns are distinct enough from GTVp patterns to avoid
        % cross-contamination.
        gtvn_patterns = {'*GTV*LN', '*GTVn', '*GTV_node*'};
        for k = 1:numel(candidate_folders)
            gtvn_path = discover_gtv_file(candidate_folders{k}, gtvn_patterns, dwii);
            if ~isempty(gtvn_path); break; end
        end
    end
end

function folders = local_candidate_gtv_folders(fxfolder, dwii)
% Return an ordered cell array of folders in which to search for GTV masks.
% The Fx folder itself is always first (highest priority). Any GTV* subfolder
% whose trailing numeric index matches dwii is added, as are non-indexed
% GTV subfolders (e.g., "GTV/", "GTVmasks/") which act as shared containers.
    folders = {fxfolder};
    if exist(fxfolder, 'dir') ~= 7
        return;
    end

    dd = dir(fxfolder);
    if isempty(dd); return; end
    isdir_mask = [dd.isdir];
    names = {dd.name};
    keep = isdir_mask & ~strcmp(names, '.') & ~strcmp(names, '..');
    dd = dd(keep);

    for k = 1:numel(dd)
        nm = dd(k).name;
        nm_lower = lower(nm);
        if ~startsWith(nm_lower, 'gtv')
            continue;
        end
        % Extract trailing integer index if present (e.g., GTVtimepoint1 -> 1)
        tok = regexp(nm_lower, '(\d+)$', 'tokens', 'once');
        if isempty(tok)
            % Non-indexed GTV container: include unconditionally
            folders{end+1} = fullfile(fxfolder, nm); %#ok<AGROW>
        else
            % Indexed container: include only if the embedded index matches dwii
            if str2double(tok{1}) == dwii
                folders{end+1} = fullfile(fxfolder, nm); %#ok<AGROW>
            end
        end
    end
end
