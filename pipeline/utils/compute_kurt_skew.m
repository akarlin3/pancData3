function [kurt_val, skew_val] = compute_kurt_skew(v, min_vox_hist)
% COMPUTE_KURT_SKEW — Compute kurtosis and skewness with minimum sample guard.
% Returns NaN if the number of finite voxels is below min_vox_hist, because
% higher-order moments are unreliable with too few data points (kurtosis
% formally requires n >= 4, but practical stability needs more).
%
% Inputs:
%   v            - numeric vector of voxel values
%   min_vox_hist - minimum number of finite values required
%
% Outputs:
%   kurt_val - kurtosis (NaN if insufficient data)
%   skew_val - skewness (NaN if insufficient data)
kurt_val = NaN;
skew_val = NaN;
if numel(v) >= min_vox_hist
    v_finite = v(~isnan(v));
    if numel(v_finite) >= min_vox_hist
        kurt_val = kurtosis(v_finite);
        skew_val = skewness(v_finite);
    end
end
end
