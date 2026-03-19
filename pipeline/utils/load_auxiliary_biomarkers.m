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

    % Sanitize and validate header columns
    header_parts = strsplit(header_line, ',');
    sanitized_header = cell(size(header_parts));
    for i = 1:numel(header_parts)
        col_name = strtrim(header_parts{i});
        % Sanitize column name: alphanumeric + underscore only
        col_name = regexprep(col_name, '[^a-zA-Z0-9_]', '');
        if isempty(col_name)
            fclose(fid);
            error('load_auxiliary_biomarkers:invalidColumnName', ...
                'Column %d has invalid name after sanitization', i);
        end
        sanitized_header{i} = col_name;
    end
    
    % Validate expected column structure
    if numel(sanitized_header) < 4
        fclose(fid);
        error('load_auxiliary_biomarkers:invalidHeader', ...
            'CSV must have at least 4 columns: patient_id, biomarker_name, timepoint, value');
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

        % Sanitize string data
        patient_id = strtrim(parts{1});
        biomarker_name = strtrim(parts{2});
        
        % Sanitize patient ID and biomarker name
        patient_id = regexprep(patient_id, '[^a-zA-Z0-9_]', '');
        biomarker_name = regexprep(biomarker_name, '[^a-zA-Z0-9_]', '');
        
        if isempty(patient_id) || isempty(biomarker_name)
            warning('load_auxiliary_biomarkers:invalidData', ...
                'Skipping line %d: invalid patient_id or biomarker_name after sanitization', line_num);
            continue;
        end

        % Strict data type validation for numeric columns
        timepoint_str = strtrim(parts{3});
        value_str = strtrim(parts{4});
        
        % Validate timepoint is numeric and positive integer
        if ~isempty(regexp(timepoint_str, '^[0-9]+$', 'once'))
            timepoint_val = str2double(timepoint_str);
            if timepoint_val <= 0 || timepoint_val ~= floor(timepoint_val)
                warning('load_auxiliary_biomarkers:invalidTimepoint', ...
                    'Skipping line %d: timepoint must be positive integer', line_num);
                continue;
            end
        else
            warning('load_auxiliary_biomarkers:invalidTimepoint', ...
                'Skipping line %d: timepoint contains non-numeric characters', line_num);
            continue;
        end
        
        % Validate value is numeric (allow decimal, negative, scientific notation)
        if ~isempty(regexp(value_str, '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$', 'once'))
            value_val = str2double(value_str);
            if isnan(value_val) || isinf(value_val)
                warning('load_auxiliary_biomarkers:invalidValue', ...
                    'Skipping line %d: value is NaN or Inf', line_num);
                continue;
            end
        else
            warning('load_auxiliary_biomarkers:invalidValue', ...
                'Skipping line %d: value is not a valid number', line_num);
            continue;
        end

        row_count = row_count + 1;
        raw_ids{row_count} = patient_id;
        raw_biomarkers{row_count} = biomarker_name;
        raw_timepoints(row_count) = timepoint_val;
        raw_values(row_count) = value_val;
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