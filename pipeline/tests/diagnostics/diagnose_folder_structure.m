function diagnose_folder_structure()
% DIAGNOSE_FOLDER_STRUCTURE  Show the actual directory layout under dataloc.
%
% Prints the subdirectory tree (up to 3 levels deep) for every patient
% folder, so the mismatch between the expected ./<P##>/Fx1/ layout and
% the real layout becomes obvious.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'utils'));
    cfg = parse_config('config.json');
    dataloc = cfg.dataloc;
    fprintf('\n=== Folder Structure Diagnostic ===\n');
    fprintf('dataloc: %s\n\n', dataloc);

    top = dir(dataloc);
    top = top([top.isdir] & ~startsWith({top.name}, '.'));

    fprintf('Top-level entries under dataloc (%d):\n', numel(top));
    for i = 1:numel(top); fprintf('  %s\n', top(i).name); end
    fprintf('\n---\n');

    % For each patient folder, list up to 2 levels of subfolders.
    for i = 1:numel(top)
        pat = top(i).name;
        patpath = fullfile(dataloc, pat);
        fprintf('\n[%s]\n', pat);
        lvl1 = dir(patpath);
        lvl1 = lvl1([lvl1.isdir] & ~startsWith({lvl1.name}, '.'));
        if isempty(lvl1)
            fprintf('  (no subdirectories)\n');
            continue;
        end
        for a = 1:numel(lvl1)
            fprintf('  %s/\n', lvl1(a).name);
            lvl2_path = fullfile(patpath, lvl1(a).name);
            lvl2 = dir(lvl2_path);
            lvl2 = lvl2([lvl2.isdir] & ~startsWith({lvl2.name}, '.'));
            for b = 1:min(numel(lvl2), 10)
                fprintf('    %s/\n', lvl2(b).name);
            end
            if numel(lvl2) > 10
                fprintf('    ... (%d more)\n', numel(lvl2) - 10);
            end
        end

        % Also list any *GTV*.mat files anywhere under this patient.
        gtvs = dir(fullfile(patpath, '**', '*GTV*.mat'));
        if ~isempty(gtvs)
            fprintf('  GTV .mat files found (any depth):\n');
            for g = 1:min(numel(gtvs), 15)
                rel = strrep(gtvs(g).folder, patpath, '');
                if startsWith(rel, filesep); rel = rel(2:end); end
                fprintf('    %s%s%s\n', rel, filesep, gtvs(g).name);
            end
            if numel(gtvs) > 15
                fprintf('    ... (%d more)\n', numel(gtvs) - 15);
            end
        else
            fprintf('  (no *GTV*.mat files anywhere)\n');
        end
    end
end
