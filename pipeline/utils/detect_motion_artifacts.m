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

        % Flag criteria: NMI < 0.6 OR dropout > 15%
        if (pv.nmi < 0.6 && v ~= b0_idx) || pv.dropout_pct > 15
            motion.flagged(v) = true;
        end
    end

    motion.n_flagged = sum(motion.flagged);
end


function nmi = compute_nmi(x, y)
% Compute normalized mutual information between two vectors
    n_bins = 32;
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);

    if length(x) < 10
        nmi = NaN;
        return;
    end

    % Marginal entropies
    hx = hist_entropy(x, n_bins);
    hy = hist_entropy(y, n_bins);

    if hx == 0 || hy == 0
        nmi = NaN;
        return;
    end

    % Joint entropy
    x_edges = linspace(min(x), max(x) + eps, n_bins + 1);
    y_edges = linspace(min(y), max(y) + eps, n_bins + 1);

    joint_hist = zeros(n_bins, n_bins);
    for i = 1:length(x)
        xi = min(n_bins, max(1, sum(x(i) >= x_edges(1:end-1))));
        yi = min(n_bins, max(1, sum(y(i) >= y_edges(1:end-1))));
        joint_hist(xi, yi) = joint_hist(xi, yi) + 1;
    end
    joint_prob = joint_hist / sum(joint_hist(:));
    joint_prob = joint_prob(joint_prob > 0);
    hxy = -sum(joint_prob .* log2(joint_prob));

    % NMI = (H(X) + H(Y)) / H(X,Y)
    nmi = (hx + hy) / hxy;
end


function h = hist_entropy(x, n_bins)
    [counts, ~] = histcounts(x, n_bins);
    prob = counts / sum(counts);
    prob = prob(prob > 0);
    h = -sum(prob .* log2(prob));
end
