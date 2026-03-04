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

    dl_provenance = struct();

    % Prefer explicit parameter passing over workspace introspection
    if nargin >= 2 && ~isempty(dl_provenance_workspace) && isstruct(dl_provenance_workspace)
        dl_provenance = dl_provenance_workspace;
        fprintf('  Using DL provenance manifest from explicit argument.\n');
    elseif nargin > 0 && ~isempty(manifest_file) && exist(manifest_file, 'file')
        % Load structurally to prevent unsafe deserialization of untrusted .mat files
        tmp = load(manifest_file, 'dl_provenance');
        if isfield(tmp, 'dl_provenance')
            dl_provenance = tmp.dl_provenance;
        end
    end

    % Provide defaults if not loaded
    if ~isfield(dl_provenance, 'dncnn_train_ids')
        dl_provenance.dncnn_train_ids = {};
    end
    if ~isfield(dl_provenance, 'ivimnet_train_ids')
        dl_provenance.ivimnet_train_ids = {};
    end

    % Flag whether the manifest was actually loaded so callers can detect
    % when the leakage guard is inactive (empty IDs + no manifest = no protection).
    if ~isfield(dl_provenance, 'manifest_loaded')
        dl_provenance.manifest_loaded = ~isempty(fieldnames(dl_provenance)) && ...
            (~isempty(dl_provenance.dncnn_train_ids) || ~isempty(dl_provenance.ivimnet_train_ids));
    end
end
