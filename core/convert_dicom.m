function bad_dwi_found = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id)
% CONVERT_DICOM Converts DWI DICOM files to NIfTI format using dcm2niix
%
%   bad_dwi_found = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id)
%   executes the dcm2niix system call and verifies that exactly 3 files
%   (.nii.gz, .bval, .bvec) were created.

    bad_dwi_found = 0;
    bt_state = warning('query', 'backtrace');
    warning('off', 'backtrace');
    lock_file = fullfile(outloc, [scanID '.lock']);
    if ~exist(fullfile(outloc, [scanID '.nii.gz']), 'file') && ~exist(lock_file, 'file')
        % Create lock file to prevent parallel workers from duplicating work
        try fclose(fopen(lock_file, 'w')); catch; end
        nii_cmd = sprintf('%s -z y -f %s -o %s %s', ...
            escape_shell_arg(dcm2nii_call), ...
            escape_shell_arg(scanID), ...
            escape_shell_arg(outloc), ...
            escape_shell_arg(dicomloc));
        [exit_status, cmd_output] = system(nii_cmd);
        if exit_status ~= 0
            warning('convert_dicom:nonZeroExit', ...
                '❌ dcm2niix returned exit code %d for %s:\n%s', exit_status, fx_id, cmd_output);
            bad_dwi_found = 1;
        end

        % Verify the expected files were successfully created
        nii_exists = exist(fullfile(outloc, [scanID '.nii.gz']), 'file');
        bval_exists = exist(fullfile(outloc, [scanID '.bval']), 'file');
        bvec_exists = exist(fullfile(outloc, [scanID '.bvec']), 'file');

        if ~(nii_exists && bval_exists && bvec_exists)
            warning('convert_dicom:missingFiles', ...
                '❌ Expected DWI files not found for %s (need .nii.gz, .bval, .bvec).', fx_id);
            bad_dwi_found = 1;
        end
        % Clean up lock file
        if exist(lock_file, 'file'), delete(lock_file); end
    end
    warning(bt_state.state, 'backtrace');
end
