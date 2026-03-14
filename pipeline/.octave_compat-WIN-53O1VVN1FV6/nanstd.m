% NANSTD  Octave-compatible shim for MATLAB's nanstd (Statistics Toolbox).
%
%   y = nanstd(x) computes the standard deviation of x, ignoring NaN values.
%   y = nanstd(x, flag) uses flag=0 for sample std (N-1, default) or
%       flag=1 for population std (N).
%   y = nanstd(x, flag, dim) operates along dimension dim.
%
%   MATLAB deprecated nanstd in R2015a in favor of std(..., 'omitnan'),
%   but older code and Octave still need this function.
%
%   Implementation approach:
%   1. Compute nanmean manually (zero out NaN, sum, divide by count).
%   2. Broadcast the mean back to original size via repmat.
%   3. Compute squared deviations (zeroing NaN positions).
%   4. Divide by (n-1) or (n) depending on the flag.
function y = nanstd(x, flag, dim)
    % Default: sample standard deviation (Bessel's correction, N-1).
    if nargin < 2 || isempty(flag)
        flag = 0;
    end
    % Default dimension: first non-singleton, matching MATLAB convention.
    if nargin < 3
        dim = find(size(x) ~= 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end

    % Count non-NaN elements per slice along dim.
    mask = ~isnan(x);
    n = sum(mask, dim);

    % Compute the mean, treating NaN as zero in the sum.
    x(isnan(x)) = 0;
    mean_x = sum(x, dim) ./ n;

    % Broadcast mean_x back to size of x for element-wise subtraction.
    % repmat is used because Octave's implicit expansion (broadcasting) was
    % not available in older versions.
    sz = size(x);
    rep_dims = ones(1, ndims(x));
    rep_dims(dim) = sz(dim);
    mean_x_rep = repmat(mean_x, rep_dims);

    % Compute deviations; zero out positions that were originally NaN so
    % they do not contribute to the sum of squared deviations.
    diff = x - mean_x_rep;
    diff(~mask) = 0;

    % Choose denominator: (n-1) for sample std, n for population std.
    % max(..., 1) prevents division by zero when n=1 with flag=0.
    if flag == 0
        denom = max(n - 1, 1);
    else
        denom = n;
    end

    y = sqrt(sum(diff.^2, dim) ./ denom);
    % Return NaN for all-NaN slices.
    y(n == 0) = NaN;
end
