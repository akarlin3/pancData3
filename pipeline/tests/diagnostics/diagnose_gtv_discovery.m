function diagnose_gtv_discovery()
% DIAGNOSE_GTV_DISCOVERY  Explain why gtv_locations is mostly empty.
%
% Takes each patient folder from config_struct.dataloc and asks
% find_gtv_files() what it would return for every Fx1 repeat (dwii=1..6),
% then lists the actual GTV*.mat files sitting in the Fx1 folder so the
% user can see the naming mismatch.
%
% This runs in seconds — no DICOM conversion, no model fitting.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'utils'));

    cfg = parse_config('config.json');
    dataloc = cfg.dataloc;
    fprintf('\n=== GTV Discovery Diagnostic ===\n');
    fprintf('dataloc: %s\n\n', dataloc);

    patlist = dir(dataloc);
    patlist = patlist([patlist.isdir] & ~startsWith({patlist.name}, '.'));

    % Look specifically at Fx1 folders.
    fx_search = {'Fx1', 'fx1', 'FX1'};

    fprintf('%-22s %-6s  %s\n', 'Patient', 'dwii', 'find_gtv_files result / actual files in Fx1');
    fprintf('%s\n', repmat('-', 1, 100));

    n_pat_ok = 0;
    n_pat_any = 0;
    for j = 1:numel(patlist)
        pat = patlist(j).name;

        % Find the Fx1 subfolder.
        fxfolder = '';
        for fs = 1:numel(fx_search)
            cand = fullfile(dataloc, pat, fx_search{fs});
            if exist(cand, 'dir'); fxfolder = cand; break; end
        end
        if isempty(fxfolder)
            continue;  % not an imaging patient folder
        end
        n_pat_any = n_pat_any + 1;

        % Ask find_gtv_files about each repeat index dwii=1..6.
        per_dwii = cell(6, 1);
        any_hit = false;
        for dwii = 1:6
            try
                [gp, ~] = find_gtv_files(fxfolder, dwii, pat);
            catch
                gp = '';
            end
            per_dwii{dwii} = gp;
            if ~isempty(gp); any_hit = true; end
        end
        if any_hit; n_pat_ok = n_pat_ok + 1; end

        % List the actual GTV*.mat files sitting in the Fx1 folder.
        actual = dir(fullfile(fxfolder, '*GTV*.mat'));
        actual_names = {actual.name};

        % Also list GTV files one level deeper (in DWI1/, DWI2/, etc.).
        deeper = dir(fullfile(fxfolder, '*', '*GTV*.mat'));
        if ~isempty(deeper)
            deeper_names = cell(numel(deeper), 1);
            for dd = 1:numel(deeper)
                rel = strrep(deeper(dd).folder, fxfolder, '');
                if startsWith(rel, filesep); rel = rel(2:end); end
                deeper_names{dd} = fullfile(rel, deeper(dd).name);
            end
            actual_names = [actual_names, deeper_names']; %#ok<AGROW>
        end

        for dwii = 1:6
            fprintf('%-22s %-6d  %s\n', pat, dwii, iif(isempty(per_dwii{dwii}), '(none)', per_dwii{dwii}));
        end
        if isempty(actual_names)
            fprintf('%-22s %-6s  %s\n', '', 'files', '(no *GTV*.mat in Fx1 folder)');
        else
            fprintf('%-22s %-6s  %s\n', '', 'files', strjoin(actual_names, ', '));
        end
        fprintf('%s\n', repmat('-', 1, 100));
    end

    fprintf('\n=== Summary ===\n');
    fprintf('Patients with Fx1 folder:                    %d\n', n_pat_any);
    fprintf('Patients where find_gtv_files found >=1 GTV: %d / %d\n', n_pat_ok, n_pat_any);
    fprintf(['\nIf the "actual files in Fx1 folder" column shows GTV*.mat files ' ...
             'but find_gtv_files returns (none), the naming convention does not\n' ...
             'match the patterns in find_gtv_files.m:\n' ...
             '  single-GTV:  *GTV_MR, *GTVp, *GTV_panc*, *GTV*\n' ...
             '  dual-GTV:    *GTV_MR, *GTVp, *GTV_panc*  (broad fallback removed)\n' ...
             'The indexed pattern expects e.g. "GTV1.mat", "GTV2.mat"; the\n' ...
             'unindexed fallback expects exactly "GTV.mat".\n']);
end

function s = iif(c, a, b)
    if c; s = a; else; s = b; end
end
