function [have_mask, mask_data] = load_nifti_mask(filepath, dwi_size, message_prefix, mask_name)
    % LOAD_NIFTI_MASK Load a NIfTI mask and validate its spatial dimensions
    %   Extracted local helper from process_single_scan.m. Loads the mask at
    %   FILEPATH (rotated to MATLAB orientation) and verifies its first three
    %   dimensions match DWI_SIZE; on mismatch the mask is rejected.
    have_mask = 0;
    mask_data = [];
    if exist(filepath, 'file')
        info = niftiinfo(filepath);
        mask_data = rot90(niftiread(info));
        fprintf('...%sLoaded %s\n', message_prefix, filepath);
        have_mask = 1;
        mask_size = size(mask_data);
        if sum(mask_size ~= dwi_size(1:3)) > 0
            have_mask = 0;
            mask_data = [];
            fprintf('size mismatch. excluding %s\n', mask_name);
        end
    end
end
