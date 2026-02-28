function dl_provenance = load_dl_provenance(manifest_file)
% LOAD_DL_PROVENANCE Loads the DL provenance manifest for rigor audit
%   Encapsulates the logic to load deep learning training ids to ensure
%   no data leakage between training and testing sets.
%
% Syntax:
%   dl_provenance = load_dl_provenance(manifest_file)
%
% Inputs:
%   manifest_file - (string/char) Full path to the dl_validation_manifest.mat
%
% Outputs:
%   dl_provenance - (struct) Structure containing dncnn_train_ids and ivimnet_train_ids
%

    dl_provenance = struct();

    % Check if passed via caller workspace
    if evalin('caller', 'exist(''dl_provenance_workspace'', ''var'')')
        dl_provenance = evalin('caller', 'dl_provenance_workspace');
        fprintf('  Using DL provenance manifest from workspace.\n');
    % Fallback to base workspace (often used in interactive sessions or tests)
    elseif evalin('base', 'exist(''dl_provenance_workspace'', ''var'')')
        dl_provenance = evalin('base', 'dl_provenance_workspace');
        fprintf('  Using DL provenance manifest from workspace.\n');
    elseif nargin > 0 && exist(manifest_file, 'file')
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
end
