function p1 = compute_histogram_laplace(vec, bin_edges)
% COMPUTE_HISTOGRAM_LAPLACE — Laplace-smoothed probability distribution.
% Adds 1 pseudo-count per bin (Laplace smoothing / additive smoothing) to
% prevent zero-probability bins that would cause log(0) = -Inf in KL
% divergence or other information-theoretic comparisons.  The smoothing
% has negligible effect when n_binned >> nbins (typical for >100 voxels).
%
% Inputs:
%   vec       - numeric vector of voxel values (NaNs are excluded)
%   bin_edges - vector of histogram bin edges
%
% Output:
%   p1 - Laplace-smoothed probability distribution (1 x nbins)
if exist('OCTAVE_VERSION', 'builtin')
    vec_f = vec(~isnan(vec));
    c1 = histc(vec_f, bin_edges);
    % histc includes a count for exact matches of the last edge; merge it
    % into the penultimate bin to match histcounts behavior.
    c1(end-1) = c1(end-1) + c1(end);
    c1 = c1(1:end-1);
else
    [c1, ~] = histcounts(vec, bin_edges);
end
n_binned = sum(c1);
nbins = length(c1);
if n_binned > 0
    % Laplace smoothing: P(bin) = (count + 1) / (total + n_bins)
    p1 = (c1 + 1) / (n_binned + nbins);
else
    p1 = zeros(size(c1));
end
end
