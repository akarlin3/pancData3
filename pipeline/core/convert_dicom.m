function bad_dwi_found = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id)
% CONVERT_DICOM Converts DWI DICOM files to NIfTI format using dcm2niix
%
%   bad_dwi_found = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id)
%   executes the dcm2niix system call and verifies that exactly 3 files
%   (.nii.gz, .bval, .bvec) were created.
%
% ANALYTICAL RATIONALE — DICOM-TO-NIFTI CONVERSION
%   DWI DICOM files store each slice and b-value as individual files with
%   vendor-specific metadata encoding. dcm2niix consolidates these into:
%     - .nii.gz: 4D volume (x, y, z, b-value) with standardized geometry
%     - .bval:   text file listing b-values in acquisition order
%     - .bvec:   text file listing diffusion gradient directions
%   NIfTI format provides consistent spatial metadata (affine transforms,
%   voxel dimensions) that is essential for subsequent operations like
%   mask overlay, dose resampling, and deformable registration.
%
%   The -z y flag requests gzip compression to reduce storage (~5-10x),
%   which is important when processing cohorts of 30+ patients with
%   multiple fractions and repeat scans.

    bad_dwi_found = 0;
    MAX_RETRIES = 2;
    
    % Save and suppress backtrace in warnings to keep console output clean
    % during batch processing. Restored at the end of the function.
    bt_state = warning('query', 'backtrace');
    warning('off', 'backtrace');

    % Lock file prevents parallel parfor workers from simultaneously
    % converting the same DICOM folder. Without this guard, two workers
    % processing the same patient could both invoke dcm2niix, leading to
    % file corruption or I/O errors on the shared output directory.
    lock_file = fullfile(outloc, [scanID '.lock']);
    % Only convert if the output NIfTI does not already exist (idempotent)
    % and no other worker holds the lock (parallel safety).
    if ~exist(fullfile(outloc, [scanID '.nii.gz']), 'file') && ~exist(lock_file, 'file')
        % Create lock file to prevent parallel workers from duplicating work.
        % Verify the lock was actually created — if it fails (permissions,
        % disk full), skip conversion to avoid a race condition.
        lock_fid = fopen(lock_file, 'w');
        if lock_fid < 0
            warning('convert_dicom:lockFailed', ...
                '⚠️ Cannot create lock file %s — skipping conversion for %s to avoid race condition.', ...
                lock_file, fx_id);
            bad_dwi_found = 1;
            return;
        end
        fclose(lock_fid);

        % Retry loop for conversion with validation
        conversion_success = false;
        for attempt = 1:MAX_RETRIES
            if attempt > 1
                warning('convert_dicom:retryAttempt', ...
                    '🔄 Retry attempt %d/%d for %s', attempt, MAX_RETRIES, fx_id);
                % Clean up any partial files from previous attempt
                cleanup_partial_files(outloc, scanID);
                pause(1); % Brief delay before retry
            end
            
            % Construct the dcm2niix command with shell-escaped arguments to
            % prevent injection attacks from paths containing spaces, quotes,
            % or special characters (common in clinical network share paths).
            nii_cmd = sprintf('%s -z y -f %s -o %s %s', ...
                escape_shell_arg(dcm2nii_call), ...
                escape_shell_arg(scanID), ...
                escape_shell_arg(outloc), ...
                escape_shell_arg(dicomloc));
            [exit_status, cmd_output] = system(nii_cmd);
            
            % Parse dcm2niix output for warnings and errors
            [has_warnings, has_errors, warning_msgs, error_msgs] = parse_dcm2niix_output(cmd_output);
            
            if exit_status ~= 0
                warning('convert_dicom:nonZeroExit', ...
                    '❌ dcm2niix returned exit code %d for %s (attempt %d):\n%s', ...
                    exit_status, fx_id, attempt, cmd_output);
                continue; % Try next attempt
            end
            
            % Report warnings but don't fail conversion unless they indicate corruption
            if has_warnings
                warning('convert_dicom:conversionWarnings', ...
                    '⚠️ dcm2niix warnings for %s:\n%s', fx_id, warning_msgs);
                % Check for corruption indicators in warnings
                if contains(warning_msgs, {'corrupt', 'invalid', 'truncated', 'bad header'}, 'IgnoreCase', true)
                    warning('convert_dicom:corruptionWarning', ...
                        '🔴 Possible DICOM corruption detected in warnings for %s', fx_id);
                    continue; % Try next attempt
                end
            end
            
            if has_errors
                warning('convert_dicom:conversionErrors', ...
                    '❌ dcm2niix errors for %s:\n%s', fx_id, error_msgs);
                continue; % Try next attempt
            end

            % Verify that all three expected output files were created and validate them
            [files_valid, validation_msg] = validate_output_files(outloc, scanID);
            
            if files_valid
                conversion_success = true;
                break; % Success - exit retry loop
            else
                warning('convert_dicom:validationFailed', ...
                    '❌ Output validation failed for %s (attempt %d): %s', ...
                    fx_id, attempt, validation_msg);
            end
        end
        
        if ~conversion_success
            warning('convert_dicom:allRetriesFailed', ...
                '❌ All conversion attempts failed for %s after %d retries', fx_id, MAX_RETRIES);
            bad_dwi_found = 1;
            % Clean up any partial files
            cleanup_partial_files(outloc, scanID);
        end
        
        % Clean up lock file after conversion completes (success or failure).
        % Lock removal happens unconditionally to avoid stale locks from
        % crashed workers that would permanently block reconversion attempts.
        if exist(lock_file, 'file'), delete(lock_file); end
    end
    % Restore the original backtrace warning state
    warning(bt_state.state, 'backtrace');
end

function [has_warnings, has_errors, warning_msgs, error_msgs] = parse_dcm2niix_output(cmd_output)
    % Parse dcm2niix command output for warnings and errors
    has_warnings = false;
    has_errors = false;
    warning_msgs = '';
    error_msgs = '';
    
    if isempty(cmd_output)
        return;
    end
    
    lines = splitlines(cmd_output);
    warning_lines = {};
    error_lines = {};
    
    for i = 1:length(lines)
        line = lines{i};
        if isempty(line), continue; end
        
        % Check for warning indicators
        if contains(line, {'warning', 'warn'}, 'IgnoreCase', true) || ...
           contains(line, {'Note:', 'WARNING:'}) || ...
           contains(line, {'corrupt', 'invalid', 'truncated'}, 'IgnoreCase', true)
            warning_lines{end+1} = line;
            has_warnings = true;
        end
        
        % Check for error indicators
        if contains(line, {'error', 'failed', 'unable', 'cannot'}, 'IgnoreCase', true) || ...
           contains(line, {'ERROR:', 'FATAL:'}) || ...
           startsWith(line, 'dcm2niix: error')
            error_lines{end+1} = line;
            has_errors = true;
        end
    end
    
    if has_warnings
        warning_msgs = strjoin(warning_lines, '\n');
    end
    
    if has_errors
        error_msgs = strjoin(error_lines, '\n');
    end
end

function [files_valid, validation_msg] = validate_output_files(outloc, scanID)
    % Validate that all expected DWI files exist and have reasonable content
    files_valid = false;
    validation_msg = '';
    
    nii_file = fullfile(outloc, [scanID '.nii.gz']);
    bval_file = fullfile(outloc, [scanID '.bval']);
    bvec_file = fullfile(outloc, [scanID '.bvec']);
    
    % Check file existence
    nii_exists = exist(nii_file, 'file');
    bval_exists = exist(bval_file, 'file');
    bvec_exists = exist(bvec_file, 'file');
    
    if ~(nii_exists && bval_exists && bvec_exists)
        missing_files = {};
        if ~nii_exists, missing_files{end+1} = '.nii.gz'; end
        if ~bval_exists, missing_files{end+1} = '.bval'; end
        if ~bvec_exists, missing_files{end+1} = '.bvec'; end
        validation_msg = sprintf('Missing files: %s', strjoin(missing_files, ', '));
        return;
    end
    
    % Check file sizes (should not be empty or suspiciously small)
    nii_info = dir(nii_file);
    bval_info = dir(bval_file);
    bvec_info = dir(bvec_file);
    
    if nii_info.bytes < 1000 % NIfTI should be at least 1KB
        validation_msg = sprintf('NIfTI file too small (%d bytes)', nii_info.bytes);
        return;
    end
    
    if bval_info.bytes < 5 % bval should have at least a few characters
        validation_msg = sprintf('bval file too small (%d bytes)', bval_info.bytes);
        return;
    end
    
    if bvec_info.bytes < 5 % bvec should have at least a few characters
        validation_msg = sprintf('bvec file too small (%d bytes)', bvec_info.bytes);
        return;
    end
    
    % Validate bval/bvec file content
    try
        % Check bval file contains numeric values
        bval_content = fileread(bval_file);
        bval_nums = str2num(bval_content); %#ok<ST2NM>
        if isempty(bval_nums) || any(bval_nums < 0) || any(bval_nums > 10000)
            validation_msg = 'bval file contains invalid values';
            return;
        end
        
        % Check bvec file has 3 rows of numeric values
        bvec_content = fileread(bvec_file);
        bvec_lines = splitlines(strtrim(bvec_content));
        bvec_lines = bvec_lines(~cellfun(@isempty, bvec_lines)); % Remove empty lines
        
        if length(bvec_lines) ~= 3
            validation_msg = sprintf('bvec file should have 3 rows, found %d', length(bvec_lines));
            return;
        end
        
        % Verify each bvec line contains numbers
        for i = 1:3
            bvec_row = str2num(bvec_lines{i}); %#ok<ST2NM>
            if isempty(bvec_row)
                validation_msg = sprintf('bvec file row %d contains no valid numbers', i);
                return;
            end
        end
        
        % Check that bval and bvec have consistent number of volumes
        num_bvals = length(bval_nums);
        num_bvecs = length(str2num(bvec_lines{1})); %#ok<ST2NM>
        if num_bvals ~= num_bvecs
            validation_msg = sprintf('bval/bvec mismatch: %d bvals vs %d bvecs', num_bvals, num_bvecs);
            return;
        end
        
    catch ME
        validation_msg = sprintf('Error validating bval/bvec files: %s', ME.message);
        return;
    end
    
    % Basic NIfTI header validation using MATLAB's niftiinfo if available
    try
        if exist('niftiinfo', 'file')
            nii_info_struct = niftiinfo(nii_file);
            if isempty(nii_info_struct.ImageSize) || any(nii_info_struct.ImageSize <= 0)
                validation_msg = 'NIfTI file has invalid image dimensions';
                return;
            end
        end
    catch ME
        validation_msg = sprintf('Error reading NIfTI header: %s', ME.message);
        return;
    end
    
    files_valid = true;
    validation_msg = 'All files validated successfully';
end

function cleanup_partial_files(outloc, scanID)
    % Clean up any partial or corrupted files from failed conversion
    files_to_clean = {
        fullfile(outloc, [scanID '.nii.gz']);
        fullfile(outloc, [scanID '.bval']);
        fullfile(outloc, [scanID '.bvec']);
        fullfile(outloc, [scanID '.json']); % dcm2niix sometimes creates this
    };
    
    for i = 1:length(files_to_clean)
        if exist(files_to_clean{i}, 'file')
            try
                delete(files_to_clean{i});
            catch
                % Ignore deletion errors - file might be locked or permissions issue
            end
        end
    end
end