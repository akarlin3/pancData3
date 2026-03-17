function motion = detect_motion_artifacts(dwi_4d, b_values, mask)
% DETECT_MOTION_ARTIFACTS  Identify problematic DWI volumes before model fitting.
%
%   For each b-value volume, computes quality metrics and flags volumes
%   with suspected motion corruption based on normalized mutual information
%   with the b=0 volume and signal dropout percentage.
%
% Inputs:
%   dwi_4d    - 4D DWI array [x, y, z, n_bvalues]
%   b_values  - Vector of b-values (length = size(dwi_4d, 4))
%   mask      - 3D binary mask (same spatial dimensions as dwi_4d)
%
% Outputs:
%   motion    - Struct with fields:
%               .per_volume     - struct array with per-volume quality metrics
%               .flagged        - logical vector of flagged (corrupted) volumes
%               .n_flagged      - number of flagged volumes
%               .noise_floor    - estimated noise floor

    dwi = double(dwi_4d);
    n_vols = size(dwi, 4);

    % Initialize output
    motion = struct();
    motion.per_volume = struct('b_value', {}, 'cv', {}, 'nmi', {}, 'dropout_pct', {});
    motion.flagged = false(n_vols, 1);
    motion.n_flagged = 0;
    motion.noise_floor = NaN;

    if n_vols < 2 || isempty(mask) || ~any(mask(:))
        return;
    end

    % Estimate noise floor from background (outside mask)
    bg_mask = ~mask;
    if ndims(dwi) == 4
        b0_idx = find(b_values == min(b_values), 1);
        b0_vol = dwi(:,:,:,b0_idx);
        bg_vals = b0_vol(bg_mask);
        bg_vals = bg_vals(isfinite(bg_vals) & bg_vals > 0);
        if ~isempty(bg_vals)
            motion.noise_floor = median(bg_vals);
        else
            motion.noise_floor = 0;
        end
    else
        motion.noise_floor = 0;
    end

    dropout_threshold = 2 * motion.noise_floor;

    % Reference: b=0 volume for NMI computation
    b0_idx = find(b_values == min(b_values), 1);
    ref_vol = dwi(:,:,:,b0_idx);
    ref_masked = ref_vol(mask);

    for v = 1:n_vols
        vol = dwi(:,:,:,v);
        vol_masked = vol(mask);
        valid = isfinite(vol_masked);

        pv = struct();
        pv.b_value = b_values(v);

        % Coefficient of variation within mask
        if sum(valid) > 1 && mean(vol_masked(valid)) > 0
            pv.cv = std(vol_masked(valid)) / mean(vol_masked(valid));
        else
            pv.cv = NaN;
        end

        % Normalized mutual information with b=0
        if v == b0_idx
            pv.nmi = 1.0;
        else
            pv.nmi = compute_nmi(ref_masked, vol_masked);
        end

        % Signal dropout percentage
        if dropout_threshold > 0 && sum(valid) > 0
            pv.dropout_pct = 100 * sum(vol_masked(valid) < dropout_threshold) / sum(valid);
        else
            pv.dropout_pct = 0;
        end

        motion.per_volume(v) = pv;

        % Flag criteria: NMI < 0.6 (or NaN, indicating too few valid
        % voxels for reliable NMI) OR dropout > 15%
        if ((pv.nmi < 0.6 || isnan(pv.nmi)) && v ~= b0_idx) || pv.dropout_pct > 15
            motion.flagged(v) = true;
        end
    end

    motion.n_flagged = sum(motion.flagged);
end


function nmi = compute_nmi(x, y)
% Compute normalized mutual information between two vectors using histcounts2.
    n_bins = 32;
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);

    if length(x) < 10
        nmi = NaN;
        return;
    end

    % Guard: constant-signal volumes collapse to a single bin, yielding
    % zero entropy and NaN NMI.  Return NaN early.
    if max(x) == min(x) || max(y) == min(y)
        nmi = NaN;
        return;
    end

    % Joint histogram via histcounts2 (robust edge handling)
    x_edges = linspace(min(x), max(x) + eps(max(abs(x))), n_bins + 1);
    y_edges = linspace(min(y), max(y) + eps(max(abs(y))), n_bins + 1);
    joint_hist = histcounts2(x, y, x_edges, y_edges);

    % Marginal histograms from joint (consistent binning)
    mx = sum(joint_hist, 2);
    my = sum(joint_hist, 1);

    n_total = sum(joint_hist(:));

    % Marginal entropies
    hx = hist_entropy_from_counts(mx, n_total);
    hy = hist_entropy_from_counts(my, n_total);

    if hx == 0 || hy == 0
        nmi = NaN;
        return;
    end

    % Joint entropy
    hxy = hist_entropy_from_counts(joint_hist(:), n_total);

    % NMI = (H(X) + H(Y)) / H(X,Y)
    nmi = (hx + hy) / hxy;
end


function h = hist_entropy_from_counts(counts, n_total)
% Compute Shannon entropy (bits) from a count vector.
    prob = counts(:) / n_total;
    prob = prob(prob > 0);
    h = -sum(prob .* log2(prob));
end
