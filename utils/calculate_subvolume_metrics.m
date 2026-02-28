function [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, has_3d, gtv_mask_3d)
% CALCULATE_SUBVOLUME_METRICS — Computes dose metrics for a sub-volume.
%
%   Syntax:
%       [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, has_3d, gtv_mask_3d)
%
%   Inputs:
%       vector       - 1D array of diffusion parameters (e.g., ADC, D, f, D*).
%       threshold    - Threshold value to define the sub-volume (vector < threshold).
%       dose_vec     - 1D array of corresponding dose values.
%       has_3d       - Logical indicating if a 3D GTV mask is available.
%       gtv_mask_3d  - 3D logical mask of the GTV.
%
%   Outputs:
%       d95          - The dose covering 95% of the sub-volume.
%       v50          - The percentage of the sub-volume receiving at least 50 Gy.

    % Default return values
    d95 = NaN;
    v50 = NaN;

    se = strel('sphere', 1);
    min_cc_voxels = 10;
    min_subvol_voxels = 100;

    mask_1d = vector < threshold;

    if has_3d
        vol_3d = false(size(gtv_mask_3d));
        vol_3d(gtv_mask_3d == 1) = mask_1d;
        vol_3d = imclose(imopen(vol_3d, se), se);
        vol_3d = bwareaopen(vol_3d, min_cc_voxels);
        mask_1d = vol_3d(gtv_mask_3d == 1);
    end

    dose_sub = dose_vec(mask_1d);

    if ~isempty(dose_sub) && sum(mask_1d) >= min_subvol_voxels
        d95 = prctile(dose_sub, 5);
        v50 = sum(dose_sub >= 50) / length(dose_sub) * 100;
    end
end
