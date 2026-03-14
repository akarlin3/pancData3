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
%     X     - Numeric matrix (n x p), where rows are patients/observations
%             and columns are DWI-derived features (e.g., median ADC,
%             perfusion fraction, dose metrics).
%
%   Outputs:
%     X_clean - Matrix with constant columns removed.
%     keep    - Logical vector (1 x p) indicating retained columns. Used by
%               the caller to map surviving features back to their original
%               names for coefficient interpretation.
%
% --- Analytical Rationale ---
% Constant columns (zero variance) arise in the DWI feature matrix for
% several reasons:
%   - A dose metric (e.g., V50Gy) may be identical across all patients if
%     the treatment protocol prescribes uniform dose coverage.
%   - An IVIM parameter may be clipped to a boundary value (e.g., f=0) for
%     all patients after outlier cleaning.
%   - After KNN imputation, a feature with mostly-missing data may be
%     filled with a single nearest-neighbor value.
%   - All-NaN columns occur when a metric is unavailable for the entire
%     cohort (e.g., GTVn metrics when no patients have nodal disease).
%
% Including constant columns in Cox proportional hazards regression
% (coxphfit) or logistic regression causes rank deficiency: the design
% matrix X'X becomes singular, and coefficient estimation either fails or
% produces arbitrarily large estimates with infinite confidence intervals.
% Removing these columns before model fitting is a standard preprocessing
% step in survival analysis pipelines.

    % --- Early exit for empty input ---
    if size(X, 1) == 0
        keep = false(1, size(X, 2));
        X_clean = X(:, keep);
        return;
    end

    % --- Variance Detection via Range ---
    % Using range (max - min) rather than var() is computationally cheaper
    % and avoids numerical precision issues: var() on nearly-constant columns
    % can produce tiny nonzero values due to floating-point arithmetic,
    % whereas range == 0 is an exact test for constancy.
    % 'omitnan' ensures that NaN-heavy columns are evaluated based on their
    % finite values only; a column with 29 NaNs and one finite value is
    % still constant (range = 0) and should be removed.
    col_range = max(X, [], 1, 'omitnan') - min(X, [], 1, 'omitnan');

    % A column is kept only if it has nonzero range AND at least one finite
    % value. The isfinite check catches all-NaN columns where max and min
    % both return NaN, making col_range = NaN (which would pass the > 0
    % test as false, but the explicit check provides clarity).
    keep = (col_range > 0) & any(isfinite(X), 1);
    X_clean = X(:, keep);
end
