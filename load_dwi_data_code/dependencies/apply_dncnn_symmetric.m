function img_isolated = apply_dncnn_symmetric(img, gtv_mask, net, expand_voxels)
% apply_dncnn_symmetric - Refactored to eliminate synthetic padding.
% 
% Extracts a bounding box expanded by exactly 'expand_voxels' in every 
% direction beyond the extremes of the GTV mask. This allows the spatial 
% kernels to process real surrounding healthy tissue gradients rather than 
% synthetic edge artifacts.
%
% Inputs:
%   img           - The original 2D or 3D raw, full-field MRI array
%   gtv_mask      - The Gross Tumor Volume mask (same dimensions as img)
%   net           - The trained DnCNN dlnetwork or SeriesNetwork
%   expand_voxels - (Optional) Number of absolute voxels to expand the 
%                   bounding box beyond the GTV extremes. Default is 15.
%
% Outputs:
%   img_isolated  - The denoised image, strictly masked to the GTV for 
%                   radiomic extraction, with surrounding values zeroed.

    if nargin < 4
        % Default expansion explicitly set to 15 voxels as required
        expand_voxels = 15; 
    end

    % 1. Get original dimensions and validate match
    orig_size = size(img);
    assert(isequal(orig_size, size(gtv_mask)), 'Image and GTV mask must have identical dimensions.');

    % 2. Extract bounding box from the raw, full-field MRI
    if ndims(img) == 2
        [row, col] = find(gtv_mask > 0);
        
        if isempty(row)
            img_isolated = img .* 0; % No GTV found
            return;
        end
        
        % Calculate expanded bounding box with rigid boundary clamping
        r_min = max(1, min(row) - expand_voxels);
        r_max = min(orig_size(1), max(row) + expand_voxels);
        c_min = max(1, min(col) - expand_voxels);
        c_max = min(orig_size(2), max(col) + expand_voxels);
        
        % Extract true anatomical region incorporating surrounding healthy tissue
        img_cropped = img(r_min:r_max, c_min:c_max);
        
    elseif ndims(img) == 3
        ind = find(gtv_mask > 0);
        if isempty(ind)
            img_isolated = img .* 0; 
            return;
        end
        [row, col, z] = ind2sub(orig_size, ind);
        
        % Calculate expanded bounding box with rigid boundary clamping
        r_min = max(1, min(row) - expand_voxels);
        r_max = min(orig_size(1), max(row) + expand_voxels);
        c_min = max(1, min(col) - expand_voxels);
        c_max = min(orig_size(2), max(col) + expand_voxels);
        z_min = max(1, min(z) - expand_voxels);
        z_max = min(orig_size(3), max(z) + expand_voxels);
        
        % Extract true anatomical region incorporating surrounding healthy tissue
        img_cropped = img(r_min:r_max, c_min:c_max, z_min:z_max);
    else
        error('Input image must be a 2D or 3D coordinate array.');
    end

    % 3. Extract the same cropped region from the GTV mask for mask-constrained normalization
    if ndims(img) == 2
        mask_cropped = gtv_mask(r_min:r_max, c_min:c_max);
    else
        mask_cropped = gtv_mask(r_min:r_max, c_min:c_max, z_min:z_max);
    end

    % 4. Calculate mu and sigma strictly usingvoxels within the GTV mask
    target_voxels = single(img_cropped(mask_cropped > 0));
    
    if ~isempty(target_voxels)
        mu_target = mean(target_voxels, 'all');
        sigma_target = std(target_voxels, 0, 'all');
        
        % Safety check: avoid divide-by-zero for uniform or very small regions
        if sigma_target < 1e-8
            sigma_target = 1.0;
        end
        
        % Apply target-derived normalization parameters to the entire bounding box
        img_cropped = (single(img_cropped) - mu_target) / sigma_target;
    end

    % 5. Pass this true, expanded anatomical region through the dnCNN
    if isa(net, 'dlnetwork') || isa(net, 'SeriesNetwork') || isa(net, 'DAGNetwork')
        img_denoised_cropped = predict(net, single(img_cropped));
        if isa(img_denoised_cropped, 'dlarray')
            img_denoised_cropped = extractdata(img_denoised_cropped);
        end
    else
        % Fallback for built-in or custom wrappers
        img_denoised_cropped = denoiseImage(img_cropped, net);
    end
    
    % 6. Matrix insertion back to full-field dimensions
    img_denoised_full = zeros(orig_size, class(img_denoised_cropped));
    
    if ndims(img) == 2
        img_denoised_full(r_min:r_max, c_min:c_max) = img_denoised_cropped;
    elseif ndims(img) == 3
        img_denoised_full(r_min:r_max, c_min:c_max, z_min:z_max) = img_denoised_cropped;
    end
    
    % 5. Apply the exact GTV mask after inference to strictly isolate tumor values
    img_isolated = img_denoised_full .* double(gtv_mask > 0);
end
