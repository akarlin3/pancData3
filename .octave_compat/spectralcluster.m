function idx = spectralcluster(X, k, varargin)
% SPECTRALCLUSTER  Octave-compatible fallback using k-means.
%
%   idx = spectralcluster(X, k) clusters the rows of X into k groups.
%
%   MATLAB's spectralcluster (introduced in R2019b) is not available in
%   Octave.  This shim provides k-means as a reasonable approximation.
%   The caller (extract_tumor_core.m) already z-score-normalizes the
%   feature matrix, so k-means in the normalized space is a reasonable
%   substitute for most tumor core segmentation tasks.

    warning('spectralcluster:octaveShim', ...
        'Using k-means fallback for spectralcluster (Octave compatibility).');
    [idx, ~] = kmeans(X, k, 'Replicates', 5);
end
