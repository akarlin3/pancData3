function [X_clean, keep] = remove_constant_columns(X)
% REMOVE_CONSTANT_COLUMNS  Remove columns that have zero variance.
%
%   [X_clean, keep] = remove_constant_columns(X)
%
%   A column is removed if all its finite values are identical (or if it
%   contains only NaN).  This prevents rank-deficiency warnings in coxphfit
%   and other regression routines.
%
%   Inputs:
%     X     - Numeric matrix (n x p).
%
%   Outputs:
%     X_clean - Matrix with constant columns removed.
%     keep    - Logical vector (1 x p) indicating retained columns.

    col_range = max(X, [], 1, 'omitnan') - min(X, [], 1, 'omitnan');
    keep = (col_range > 0) & any(isfinite(X), 1);
    X_clean = X(:, keep);
end
