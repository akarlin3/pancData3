function bad_dwi_found = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id)
% CONVERT_DICOM Converts DWI DICOM files to NIfTI format using dcm2niix
%
%   bad_dwi_found = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id)
%   executes the dcm2niix system call and verifies that exactly 3 files
%   (.nii.gz, .bval, .bvec) were created.

    bad_dwi_found = 0;
    if ~exist(fullfile(outloc, [scanID '.nii.gz']), 'file')
        nii_cmd = sprintf('%s -z y -f %s -o %s %s', ...
            escape_shell_arg(dcm2nii_call), ...
            escape_shell_arg(scanID), ...
            escape_shell_arg(outloc), ...
            escape_shell_arg(dicomloc));
        system(nii_cmd);

        % Verify the expected files were successfully created
        nii_exists = exist(fullfile(outloc, [scanID '.nii.gz']), 'file');
        bval_exists = exist(fullfile(outloc, [scanID '.bval']), 'file');
        bvec_exists = exist(fullfile(outloc, [scanID '.bvec']), 'file');

        if ~(nii_exists && bval_exists && bvec_exists)
            fprintf('!!!! incorrect number of DWI files generated for %s... need to fix!\n', fx_id);
            bad_dwi_found = 1;
        end
    end
end
