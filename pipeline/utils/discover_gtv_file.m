function filepath = discover_gtv_file(folder, patterns, index)
% DISCOVER_GTV_FILE Locate a GTV mask file using flexible naming patterns.
%
% --- Analytical Rationale ---
% In pancreatic DWI studies, the Gross Tumor Volume (GTV) is manually contoured
% by radiation oncologists on each MRI scan. These contours are exported as 3D
% binary masks stored in .mat files. However, different clinical sites, physicians,
% and contouring software produce files with inconsistent naming conventions:
%   - Some use indexed names: GTV1.mat, GTV2.mat (one per scan repeat)
%   - Some append registration suffixes: GTV1_reg.mat, GTV1_deformed.mat
%   - Some use a single unindexed file: GTV.mat (shared across all repeats of
%     the same fraction, when the patient did not move between acquisitions)
%
% The 'index' parameter corresponds to the repeatability index (scan repeat
% within the same treatment fraction). Repeat scans are acquired back-to-back
% to assess measurement reproducibility of DWI parameters, which is essential
% for determining whether observed longitudinal changes exceed measurement noise.
%
% This function implements a priority-based search: indexed matches are preferred
% over unindexed ones, because indexed masks correspond to specific scan
% acquisitions and may reflect slight patient repositioning between repeats.

    % --- Normalize Input ---
    % Accept both single pattern strings and cell arrays for flexibility.
    % Callers like find_gtv_files pass cell arrays of naming conventions to
    % try in priority order.
    if ischar(patterns) || isstring(patterns)
        patterns = {char(patterns)};
    end
    filepath = '';

    % --- Initialize Empty Dir Struct Arrays ---
    % dir() on a nonexistent file returns an empty struct with the correct
    % fields (name, folder, date, bytes, isdir, datenum). This provides a
    % properly-typed empty array for concatenation below, avoiding errors
    % when appending dir() results with different struct field orders.
    gtv_search = dir(fullfile(folder, '__nonexistent__'));  % empty struct array with dir fields
    single_gtv_search = gtv_search;

    for p = 1:length(patterns)
        pat = patterns{p};

        % --- Indexed Search (Scan-Repeat-Specific Masks) ---
        % Search for exact index match first (e.g., GTV1.mat), then
        % index followed by a non-digit separator (e.g., GTV1_reg.mat).
        % Avoid [pat index '*.mat'] which for index=1 also matches
        % index 10, 11, 12, etc.
        % This precision matters because incorrectly pairing GTV masks with
        % DWI data would corrupt the voxel-level diffusion parameter extraction,
        % since the mask defines which voxels belong to tumor vs. normal tissue.
        tmp = dir(fullfile(folder, [pat int2str(index) '.mat']));
        if isempty(tmp)
            % Try the variant with a separator after the index number.
            % Registration or post-processing tools often append suffixes like
            % '_reg', '_warped', '_resampled' to indicate the mask has been
            % aligned to the DWI coordinate space.
            tmp = dir(fullfile(folder, [pat int2str(index) '_*.mat']));
        end
        if ~isempty(tmp); gtv_search = [gtv_search; tmp]; end

        % --- Unindexed Search (Single Mask for All Repeats) ---
        % Some sites provide only one GTV mask per fraction when the tumor
        % did not move between repeat acquisitions. This is the fallback
        % when no indexed mask is found.
        tmp = dir(fullfile(folder, [pat '.mat']));
        if ~isempty(tmp); single_gtv_search = [single_gtv_search; tmp]; end
    end

    % --- Priority Resolution ---
    % Prefer indexed (scan-specific) masks over unindexed (shared) masks.
    % When multiple patterns match, take the first result -- the patterns
    % are ordered from most-specific to least-specific by the caller, so
    % the first match is the best match.
    if isempty(gtv_search)
        % fallback: sometimes only 1 mask exists for all repeats
        if ~isempty(single_gtv_search)
            filepath = fullfile(single_gtv_search(1).folder, single_gtv_search(1).name);
        end
    else
        filepath = fullfile(gtv_search(1).folder, gtv_search(1).name);
    end
end