function [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, pat_name)
% FIND_GTV_FILES Locate GTVp and (optionally) GTVn mask files for a scan.
%
%   [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, pat_name)
%
%   Inputs:
%       fxfolder - Absolute path to the fraction folder containing mask files
%       dwii     - Repeatability index (1-based scan repeat within the fraction)
%       pat_name - Patient identifier string; the token 'two' in the name
%                  signals that this patient has both a primary pancreatic
%                  tumor (GTVp) and lymph node metastases (GTVn)
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

    gtvp_path = '';
    gtvn_path = '';

    if ~contains(pat_name, 'two')
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
        gtvp_path = discover_gtv_file(fxfolder, gtvp_patterns, dwii);
    else
        % --- Dual-GTV Patient (Primary Tumor + Nodal Metastasis) ---
        % The broad '*GTV*' fallback is intentionally omitted from GTVp
        % patterns here. If it were included, it could match a GTVn file
        % (e.g., 'GTVn1.mat'), incorrectly labeling nodal tissue as the
        % primary tumor and corrupting downstream dose-response analysis.
        gtvp_patterns = {'*GTV_MR', '*GTVp', '*GTV_panc*'};
        gtvp_path = discover_gtv_file(fxfolder, gtvp_patterns, dwii);

        % Nodal GTV naming conventions vary by contouring physician:
        %   - GTV_LN  = "Gross Tumor Volume - Lymph Node"
        %   - GTVn    = Standard ICRU-style nodal GTV label
        %   - GTV_node = Descriptive alternative used at some sites
        % These patterns are distinct enough from GTVp patterns to avoid
        % cross-contamination.
        gtvn_patterns = {'*GTV*LN', '*GTVn', '*GTV_node*'};
        gtvn_path = discover_gtv_file(fxfolder, gtvn_patterns, dwii);
    end
end
