function [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, has_3d, gtv_mask_3d, direction)
% CALCULATE_SUBVOLUME_METRICS — Computes dose metrics for a sub-volume
%   defined by a diffusion parameter threshold within the GTV.
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
%
%   Analytical Rationale — DWI-Guided Dose Assessment:
%   ---------------------------------------------------
%   This function addresses a central question in adaptive RT for
%   pancreatic cancer: "Is the most biologically aggressive part of the
%   tumor receiving adequate radiation dose?"
%
%   The sub-volume is defined by a diffusion parameter threshold:
%   - ADC < threshold identifies highly cellular (restricted diffusion)
%     regions — likely viable, treatment-resistant tumor tissue.
%   - f < threshold identifies hypovascular regions — poorly perfused
%     stroma that may be hypoxic and radioresistant.
%   - ADC > threshold (direction='above') identifies regions showing
%     treatment response (cell death increasing free water diffusion).
%
%   D95 (dose to 95% of the sub-volume) is the standard metric for
%   target coverage in radiation oncology — it answers "what is the
%   minimum dose that 95% of this sub-volume receives?"  A low D95
%   indicates cold spots within the biologically aggressive sub-region.
%
%   V50 (volume receiving >= 50 Gy) evaluates whether the sub-volume
%   receives the typical curative prescription dose for pancreatic SBRT.
%   50 Gy in 5 fractions (BED10 = 100 Gy) is a standard prescription.

    if nargin < 6 || isempty(direction), direction = 'below'; end

    % Default return values
    d95 = NaN;
    v50 = NaN;

    % Morphological structuring element for noise cleanup.  A sphere of
    % radius 1 (6-connected neighborhood in 3D) is used for opening/closing
    % to remove isolated voxels that pass the threshold due to noise rather
    % than genuine biological signal.  This is important because IVIM
    % fitting noise can produce spurious low-ADC voxels in healthy tissue,
    % inflating the apparent size of the aggressive sub-volume.
    if exist('OCTAVE_VERSION', 'builtin')
        % strel('sphere',1) produces a 7-voxel diamond (6-connected).
        % ones(3,3,3) is a 27-voxel cube (26-connected) — not equivalent.
        % Build the 6-connected diamond kernel manually for Octave compat.
        sphere_kernel = zeros(3, 3, 3);
        sphere_kernel(2,2,:) = 1; sphere_kernel(2,:,2) = 1; sphere_kernel(:,2,2) = 1;
        se = strel('arbitrary', sphere_kernel);
    else
        se = strel('sphere', 1);
    end
    % Minimum connected component size: clusters smaller than 10 voxels
    % are likely noise artifacts rather than coherent tumor sub-regions.
    % At typical DWI resolution (2-3 mm isotropic), 10 voxels ≈ 0.08-0.27
    % cm^3 — below the resolution limit for meaningful dose assessment.
    min_cc_voxels = 10;
    % Minimum sub-volume size for reliable dose statistics.  With fewer
    % than 100 voxels, D95 and V50 become highly sensitive to individual
    % voxel values and partial-volume effects at the GTV boundary.
    min_subvol_voxels = 100;

    if strcmpi(direction, 'above')
        mask_1d = vector > threshold;
    else
        mask_1d = vector < threshold;
    end

    if has_3d
        % When a 3D GTV mask is available, we can perform spatially-aware
        % sub-volume refinement.  The 1D threshold mask is embedded back
        % into the 3D GTV geometry to leverage spatial context:
        %
        % 1. imopen (erosion then dilation): removes isolated threshold-
        %    passing voxels that are not part of a coherent spatial cluster.
        %    These typically arise from IVIM fitting noise at voxels near
        %    the GTV boundary where partial-volume effects are worst.
        %
        % 2. imclose (dilation then erosion): fills small holes within
        %    contiguous sub-regions, producing smoother boundaries that
        %    better represent the underlying biological sub-volume.
        %
        % 3. bwareaopen: removes connected components smaller than
        %    min_cc_voxels, keeping only spatially coherent sub-regions
        %    that are large enough for meaningful dose assessment.
        %
        % After 3D cleanup, the refined mask is projected back to 1D for
        % dose extraction.
        n_gtv_voxels = sum(gtv_mask_3d(:) == 1);
        assert(n_gtv_voxels == numel(mask_1d), ...
            'calculate_subvolume_metrics:dimMismatch', ...
            'mask_1d length (%d) does not match GTV voxel count (%d).', ...
            numel(mask_1d), n_gtv_voxels);
        vol_3d = false(size(gtv_mask_3d));
        vol_3d(gtv_mask_3d == 1) = mask_1d;
        vol_3d = imclose(imopen(vol_3d, se), se);
        vol_3d = bwareaopen(vol_3d, min_cc_voxels);
        mask_1d = vol_3d(gtv_mask_3d == 1);
    end

    % Extract dose values for the sub-volume voxels
    dose_sub_raw = dose_vec(mask_1d);
    n_total = numel(dose_sub_raw);
    dose_sub = dose_sub_raw(~isnan(dose_sub_raw));
    n_valid = numel(dose_sub);

    if ~isempty(dose_sub) && n_valid >= min_subvol_voxels
        % D95 = dose covering 95% of the sub-volume = 5th percentile of
        % the dose distribution.  The 5th percentile gives the minimum dose
        % that at least 95% of voxels exceed — the standard definition of
        % D95 in radiation oncology DVH analysis.  A low D95 in the
        % restricted-diffusion sub-volume suggests under-dosing of the most
        % aggressive tumor component.
        d95 = prctile(dose_sub, 5);
        % V50 uses all subvolume voxels (including NaN) as denominator so
        % that unknown-dose voxels conservatively deflate coverage.  Using
        % only valid voxels would inflate V50 when dose data is sparse.
        v50 = sum(dose_sub >= 50) / max(n_total, 1) * 100;
        nan_frac = 1 - n_valid / max(n_total, 1);
        if nan_frac > 0.2
            warning('calculate_subvolume_metrics:highNaNFrac', ...
                '%.0f%% of subvolume voxels have NaN dose — V50 may be deflated.', ...
                nan_frac * 100);
        end
    end
end
