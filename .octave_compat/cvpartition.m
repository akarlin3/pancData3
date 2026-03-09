% CVPARTITION  Octave-compatible stub for MATLAB's cvpartition.
%
%   MATLAB's cvpartition (Statistics Toolbox) creates cross-validation
%   partition objects for K-Fold, Leave-One-Out, Holdout, etc. This stub
%   returns a minimal struct with a NumTestSets field so that code which
%   queries the number of folds does not error.
%
%   Behavioral differences from MATLAB's cvpartition:
%   - Does not actually partition the data; training() and test() methods
%     are not implemented.
%   - For 'KFold' style, always reports 5 folds regardless of the 'k'
%     argument (sufficient for code paths that only read NumTestSets).
%   - For group-based partitioning, sets NumTestSets to the number of
%     unique groups as a rough approximation.
%   - This is a stub, not a full implementation; callers that need actual
%     fold indices will need a more complete shim.
function c = cvpartition(n, varargin)
    % Return a struct mimicking cvpartition's public properties.
    c = struct();
    c.NumTestSets = 1;
    if ischar(n) || isstring(n)
        % Group-based partitioning, e.g. cvpartition(group, 'Leaveout').
        % Set NumTestSets to the number of unique groups.
        groups = varargin{1};
        c.NumTestSets = length(unique(groups));
    else
        % Numeric N with 'KFold' or similar -- default to 5 folds.
        c.NumTestSets = 5;
    end
end
