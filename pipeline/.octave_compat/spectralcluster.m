% SPECTRALCLUSTER  Octave-compatible fallback for MATLAB's spectralcluster (R2019b+).
%
%   idx = spectralcluster(X, k) clusters the rows of X into k groups.
%
%   MATLAB's spectralcluster constructs a similarity graph from the data,
%   computes the graph Laplacian's eigenvectors, and then applies k-means
%   in the spectral embedding space. Octave does not provide this function.
%
%   This shim falls back to plain k-means clustering, which is a reasonable
%   approximation when the feature matrix has already been z-score normalized
%   (as extract_tumor_core.m does before calling this function). The 'Replicates'
%   parameter runs k-means 5 times with different random initializations to
%   reduce sensitivity to the starting centroids.
%
%   Behavioral differences from MATLAB's spectralcluster:
%   - No similarity graph or Laplacian eigenvector computation; clusters are
%     found directly in the input feature space via k-means.
%   - May produce different cluster assignments than true spectral clustering,
%     especially for non-convex or manifold-structured data.
%   - Additional name-value arguments (e.g., 'SimilarityGraph', 'KNNGraphK')
%     are accepted via varargin but silently ignored.
%   - Emits a warning on each call so users know the fallback is in effect.
function idx = spectralcluster(X, k, varargin)
    % Warn that this is a k-means approximation, not true spectral clustering.
    warning('spectralcluster:octaveShim', ...
        'Using k-means fallback for spectralcluster (Octave compatibility).');
    % Run k-means with 5 replicates for robustness against random initialization.
    [idx, ~] = kmeans(X, k, 'Replicates', 5);
end
