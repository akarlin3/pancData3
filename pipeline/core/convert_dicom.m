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

        % Construct the dcm2niix command with shell-escaped arguments to
        % prevent injection attacks from paths containing spaces, quotes,
        % or special characters (common in clinical network share paths).
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

        % Verify that all three expected output files were created.
        % DWI analysis requires all three: the volume (.nii.gz), b-values
        % (.bval) for model fitting, and gradient directions (.bvec) for
        % directional diffusion analysis. A missing .bval file makes IVIM
        % and ADC fitting impossible; a missing .bvec prevents directional
        % analysis (though this pipeline uses trace/isotropic DWI).
        nii_exists = exist(fullfile(outloc, [scanID '.nii.gz']), 'file');
        bval_exists = exist(fullfile(outloc, [scanID '.bval']), 'file');
        bvec_exists = exist(fullfile(outloc, [scanID '.bvec']), 'file');

        if ~(nii_exists && bval_exists && bvec_exists)
            warning('convert_dicom:missingFiles', ...
                '❌ Expected DWI files not found for %s (need .nii.gz, .bval, .bvec).', fx_id);
            bad_dwi_found = 1;
        end
        % Clean up lock file after conversion completes (success or failure).
        % Lock removal happens unconditionally to avoid stale locks from
        % crashed workers that would permanently block reconversion attempts.
        if exist(lock_file, 'file'), delete(lock_file); end
    end
    % Restore the original backtrace warning state
    warning(bt_state.state, 'backtrace');
end
