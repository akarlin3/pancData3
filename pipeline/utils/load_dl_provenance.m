function dl_provenance = load_dl_provenance(manifest_file, dl_provenance_workspace)
% LOAD_DL_PROVENANCE Loads the DL provenance manifest for rigor audit
%   Encapsulates the logic to load deep learning training ids to ensure
%   no data leakage between training and testing sets.
%
% Syntax:
%   dl_provenance = load_dl_provenance(manifest_file)
%   dl_provenance = load_dl_provenance(manifest_file, dl_provenance_workspace)
%
% Inputs:
%   manifest_file             - (string/char) Full path to the dl_validation_manifest.mat
%   dl_provenance_workspace   - (Optional struct) Pre-loaded provenance struct.
%                                Replaces the old evalin('caller'/'base') pattern,
%                                which violated the project's explicit parameter-passing convention.
%
% Outputs:
%   dl_provenance - (struct) Structure containing dncnn_train_ids and ivimnet_train_ids
%
%   Analytical Rationale — Why DL Provenance Tracking is Critical:
%   ---------------------------------------------------------------
%   The dnCNN and IVIMnet deep learning models used for DWI denoising and
%   parameter estimation were trained on a subset of patients from the same
%   institution.  If any of those training patients appear in the
%   downstream statistical analysis (survival modeling, predictive
%   performance evaluation), the results are contaminated by circular
%   reasoning: the DL model has already "seen" those patients' data during
%   its training, so it will produce artificially clean/accurate parameter
%   estimates for them, biasing the analysis toward showing that DL
%   denoising improves outcomes.
%
%   This function loads the manifest of patient IDs used to train each DL
%   model, enabling the downstream pipeline to either (a) exclude those
%   patients from analysis, or (b) flag the analysis as potentially
%   contaminated.  The manifest_loaded flag distinguishes "no manifest
%   file found" (dangerous — leakage guard is inactive) from "manifest
%   found but lists are empty" (benign — no patients to exclude, e.g.,
%   when using a pre-trained model from a different institution).
%

    dl_provenance = struct();

    % Prefer explicit parameter passing over workspace introspection.
    % This follows the project's strict no-global-state convention:
    % all data flows through function arguments, never through evalin()
    % or assignin() which create hidden dependencies and make the code
    % difficult to test and reason about.
    if nargin >= 2 && ~isempty(dl_provenance_workspace) && isstruct(dl_provenance_workspace)
        dl_provenance = dl_provenance_workspace;
        fprintf('  Using DL provenance manifest from explicit argument.\n');
    elseif nargin > 0 && ~isempty(manifest_file) && exist(manifest_file, 'file')
        % Load only the 'dl_provenance' variable by name to avoid
        % deserializing other potentially unsafe variables in the manifest
        % file.  This is a defense-in-depth measure complementing
        % safe_load_mask for non-mask .mat files.
        tmp = load(manifest_file, 'dl_provenance');
        if isfield(tmp, 'dl_provenance')
            dl_provenance = tmp.dl_provenance;
        end
    end

    % Provide defaults if not loaded.  Empty cell arrays mean "no patients
    % to exclude" — this is the safe default because it does not suppress
    % any patients from the analysis.  The separate manifest_loaded flag
    % (below) alerts callers when the manifest was not found, so they can
    % decide whether to proceed with unguarded analysis or halt.
    %
    % Two separate lists are maintained because the DnCNN denoiser and
    % IVIMnet parameter estimator may have been trained on different
    % patient subsets (e.g., DnCNN trained on all patients, IVIMnet
    % trained on a held-out subset).
    if ~isfield(dl_provenance, 'dncnn_train_ids')
        dl_provenance.dncnn_train_ids = {};
    end
    if ~isfield(dl_provenance, 'ivimnet_train_ids')
        dl_provenance.ivimnet_train_ids = {};
    end

    % Flag whether the manifest was actually loaded so callers can detect
    % when the leakage guard is inactive (empty IDs + no manifest = no protection).
    % Distinguish "manifest file not found" (dangerous — no leakage protection)
    % from "manifest found but empty" (benign — no DL patients to guard against).
    if ~isfield(dl_provenance, 'manifest_loaded')
        manifest_file_existed = (nargin >= 2 && ~isempty(dl_provenance_workspace) && isstruct(dl_provenance_workspace)) || ...
            (nargin > 0 && ~isempty(manifest_file) && exist(manifest_file, 'file'));
        dl_provenance.manifest_loaded = manifest_file_existed;
    end
end
