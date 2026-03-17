function aux_data = load_auxiliary_biomarkers(csv_path, id_list)
% LOAD_AUXILIARY_BIOMARKERS  Load non-DWI biomarker data from CSV.
%
%   Parses a CSV file with columns: patient_id, biomarker_name, timepoint, value.
%   Pivots to a patient × timepoint matrix per biomarker. Returns a struct
%   with one field per unique biomarker name.
%
% Inputs:
%   csv_path  - Path to the auxiliary biomarker CSV file
%   id_list   - Cell array of pipeline patient IDs for validation
%
% Outputs:
%   aux_data  - Struct with one field per biomarker, each a [n_patients x max_tp] matrix

    aux_data = struct();

    if ~exist(csv_path, 'file')
        error('load_auxiliary_biomarkers:fileNotFound', ...
            'Auxiliary biomarker CSV not found: %s', csv_path);
    end

    % Read CSV
    fid = fopen(csv_path, 'r');
    if fid < 0
        error('load_auxiliary_biomarkers:readError', ...
            'Cannot open CSV file: %s', csv_path);
    end

    % Read header
    header_line = fgetl(fid);
    if ~ischar(header_line)
        fclose(fid);
        return;
    end

    % Count data lines for pre-allocation (avoids repeated array growing)
    file_pos = ftell(fid);
    n_lines_est = 0;
    while ~feof(fid)
        fgetl(fid);
        n_lines_est = n_lines_est + 1;
    end
    fseek(fid, file_pos, 'bof');

    % Pre-allocate arrays
    raw_ids = cell(1, n_lines_est);
    raw_biomarkers = cell(1, n_lines_est);
    raw_timepoints = zeros(1, n_lines_est);
    raw_values = zeros(1, n_lines_est);
    line_num = 1;
    row_count = 0;

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(strtrim(line))
            continue;
        end
        line_num = line_num + 1;

        parts = strsplit(line, ',');
        if numel(parts) < 4
            warning('load_auxiliary_biomarkers:malformedLine', ...
                'Skipping malformed line %d', line_num);
            continue;
        end

        row_count = row_count + 1;
        raw_ids{row_count} = strtrim(parts{1});
        raw_biomarkers{row_count} = strtrim(parts{2});
        raw_timepoints(row_count) = str2double(strtrim(parts{3}));
        raw_values(row_count) = str2double(strtrim(parts{4}));
    end
    fclose(fid);

    % Trim to actual size
    raw_ids = raw_ids(1:row_count);
    raw_biomarkers = raw_biomarkers(1:row_count);
    raw_timepoints = raw_timepoints(1:row_count);
    raw_values = raw_values(1:row_count);

    if isempty(raw_ids)
        return;
    end

    % Get unique biomarkers
    unique_biomarkers = unique(raw_biomarkers, 'stable');
    max_tp = max(raw_timepoints);
    n_patients = numel(id_list);

    % Check for patient ID mismatches
    unique_csv_ids = unique(raw_ids, 'stable');
    unmatched = setdiff(unique_csv_ids, id_list);
    if ~isempty(unmatched)
        warning('load_auxiliary_biomarkers:unmatchedPatients', ...
            '%d patient IDs in CSV do not match pipeline cohort: %s', ...
            numel(unmatched), strjoin(unmatched(1:min(5,end)), ', '));
    end

    % Build patient index map
    pat_map = containers.Map();
    for p = 1:n_patients
        pat_map(id_list{p}) = p;
    end

    % Pivot per biomarker
    for b = 1:numel(unique_biomarkers)
        bname = unique_biomarkers{b};
        % Sanitize field name
        field_name = matlab.lang.makeValidName(bname);
        mat = nan(n_patients, max_tp);

        bm_idx = strcmp(raw_biomarkers, bname);
        bm_ids = raw_ids(bm_idx);
        bm_tps = raw_timepoints(bm_idx);
        bm_vals = raw_values(bm_idx);

        for i = 1:numel(bm_ids)
            pid = bm_ids{i};
            if pat_map.isKey(pid)
                row = pat_map(pid);
                tp = bm_tps(i);
                if tp >= 1 && tp <= max_tp
                    if ~isnan(mat(row, tp))
                        warning('load_auxiliary_biomarkers:duplicateEntry', ...
                            'Duplicate entry for %s/%s/tp%d — taking last value.', ...
                            pid, bname, tp);
                    end
                    mat(row, tp) = bm_vals(i);
                end
            end
        end

        aux_data.(field_name) = mat;
    end

    fprintf('  📁 Loaded %d auxiliary biomarkers from %s\n', ...
        numel(unique_biomarkers), csv_path);
    for b = 1:numel(unique_biomarkers)
        field_name = matlab.lang.makeValidName(unique_biomarkers{b});
        mat = aux_data.(field_name);
        n_valid = sum(~isnan(mat(:)));
        fprintf('    %s: %d values across %d patients\n', ...
            unique_biomarkers{b}, n_valid, sum(any(~isnan(mat), 2)));
    end
end
