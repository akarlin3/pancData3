function [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, has_3d, gtv_mask_3d, direction)
% CALCULATE_SUBVOLUME_METRICS — Computes dose metrics for a sub-volume.
%
%   Syntax:
%       [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, has_3d, gtv_mask_3d)
%       [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, has_3d, gtv_mask_3d, direction)
%
%   Inputs:
%       vector       - 1D array of diffusion parameters (e.g., ADC, D, f, D*).
%       threshold    - Threshold value to define the sub-volume.
%       dose_vec     - 1D array of corresponding dose values.
%       has_3d       - Logical indicating if a 3D GTV mask is available.
%       gtv_mask_3d  - 3D logical mask of the GTV.
%       direction    - (Optional) 'below' (default) selects vector < threshold;
%                      'above' selects vector > threshold.
%
%   Outputs:
%       d95          - The dose covering 95% of the sub-volume.
%       v50          - The percentage of the sub-volume receiving at least 50 Gy.

    if nargin < 6 || isempty(direction), direction = 'below'; end

    % Default return values
    d95 = NaN;
    v50 = NaN;

    if exist('OCTAVE_VERSION', 'builtin')
        se = strel('arbitrary', ones(3, 3, 3));
    else
        se = strel('sphere', 1);
    end
    min_cc_voxels = 10;
    min_subvol_voxels = 100;

    if strcmpi(direction, 'above')
        mask_1d = vector > threshold;
    else
        mask_1d = vector < threshold;
    end

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
