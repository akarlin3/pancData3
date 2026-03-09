% NANMEAN  Octave-compatible shim for MATLAB's nanmean (Statistics Toolbox).
%
%   y = nanmean(x) computes the mean of x, ignoring NaN values.
%   y = nanmean(x, dim) operates along dimension dim.
%
%   MATLAB deprecated nanmean in R2015a in favor of mean(..., 'omitnan'),
%   but older code and Octave (which lacks the 'omitnan' flag in some
%   versions) still need this function.
%
%   Implementation: replaces NaN with 0, sums, and divides by the count
%   of non-NaN elements. Returns NaN for slices that are entirely NaN.
function y = nanmean(x, dim)
    % Default dimension: first non-singleton, matching MATLAB's mean() behavior.
    if nargin == 1
        dim = find(size(x) ~= 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end
    % Build a mask of valid (non-NaN) elements and count them per slice.
    mask = ~isnan(x);
    n = sum(mask, dim);
    % Zero out NaN entries so they don't contaminate the sum.
    x(isnan(x)) = 0;
    y = sum(x, dim) ./ n;
    % Slices with no valid data should return NaN, not 0/0 = NaN (already NaN)
    % or Inf. Explicitly set to NaN for clarity.
    y(n == 0) = NaN;
end
