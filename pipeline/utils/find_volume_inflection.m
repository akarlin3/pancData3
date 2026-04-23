function [knee_thresh, knee_idx, knee_curvature, curvature] = find_volume_inflection(thresholds, vol_frac)
%FIND_VOLUME_INFLECTION  Knee of a monotonic curve via discrete 2nd derivative.
%
%   [knee_thresh, knee_idx, knee_curvature, curvature] = ...
%       find_volume_inflection(thresholds, vol_frac)
%
%   Finds the saturation knee of vol_frac(thresholds) — the point of
%   most concave curvature in a 3-point smoothed copy of the curve.
%   The algorithm:
%     1. Require at least 3 finite samples in vol_frac (below this
%        floor the second derivative is undefined); otherwise return
%        NaN everywhere.
%     2. Apply NaN-aware 3-point moving-average smoothing.  Boundary
%        points use a 2-point average so the smoothed vector has the
%        same length as the input.
%     3. Compute the discrete second derivative at interior points
%        (indices 2..n-1); endpoints remain NaN in `curvature`.
%        Uniform threshold spacing is assumed — the finite difference
%        is unscaled, but this is fine because argmax of the unscaled
%        difference equals argmax of the true d^2V/dt^2.
%     4. The knee index is argmin(curvature): the most concave point
%        (most negative second derivative = saturation "knee").  NaN
%        entries in curvature are skipped.
%
%   Inputs:
%       thresholds - [1xN] monotonically-increasing threshold values.
%       vol_frac   - [1xN] monotonic (ish) cumulative fraction curve
%                    evaluated at those thresholds.  NaN entries are
%                    tolerated; endpoints are automatically excluded
%                    from the second-derivative computation.
%
%   Outputs:
%       knee_thresh    - threshold value at the knee (NaN if undefined)
%       knee_idx       - index of the knee in `thresholds` (NaN if undefined)
%       knee_curvature - second-derivative value at the knee (NaN if
%                        undefined); negative => concave down = saturation.
%       curvature      - [1xN] discrete second derivative, NaN at endpoints
%                        and at any interior point whose 3-point smoothed
%                        window is incomplete.

    n = numel(vol_frac);
    curvature = nan(1, n);
    knee_thresh    = NaN;
    knee_idx       = NaN;
    knee_curvature = NaN;

    if n < 3
        return;
    end

    finite_mask = ~isnan(vol_frac);
    if sum(finite_mask) < 3
        return;
    end

    % 3-point moving-average smoothing (NaN-aware).  Boundary points
    % use a 2-point average so the smoothed vector has the same length
    % as the input.
    sm = nan(1, n);
    v = vol_frac;
    for i = 1:n
        lo = max(1, i - 1);
        hi = min(n, i + 1);
        win = v(lo:hi);
        win = win(~isnan(win));
        if ~isempty(win)
            sm(i) = mean(win);
        end
    end

    % Discrete second derivative at interior points.
    for i = 2:n-1
        if ~isnan(sm(i-1)) && ~isnan(sm(i)) && ~isnan(sm(i+1))
            curvature(i) = sm(i+1) - 2*sm(i) + sm(i-1);
        end
    end

    if all(isnan(curvature))
        return;
    end
    [~, knee_idx] = min(curvature);
    knee_thresh    = thresholds(knee_idx);
    knee_curvature = curvature(knee_idx);
end
