function p = perform_statistical_test(data, groups, test_type)
% PERFORM_STATISTICAL_TEST - Handles testing and p-value extraction safely
%
% Syntax:
%   p = perform_statistical_test(data, groups, test_type)
%
% Inputs:
%   data      - Numeric array of data points
%   groups    - Numeric or categorical array of group labels
%   test_type - String specifying the test type (default: 'ranksum')
%
% Outputs:
%   p         - p-value of the test (NaN if insufficient data)
%
%   Analytical Rationale — Why Wilcoxon Rank-Sum is the Default:
%   -------------------------------------------------------------
%   DWI-derived parameters (ADC, D, f, D*) in pancreatic tumors rarely
%   follow normal distributions.  Tumor heterogeneity creates skewed
%   distributions (especially f and D*), and the small sample sizes
%   typical of single-institution pancreatic cancer studies (N ~ 30-80
%   per group) make it difficult to verify normality assumptions.
%
%   The Wilcoxon rank-sum (Mann-Whitney U) test is a non-parametric
%   alternative to the two-sample t-test that:
%   - Makes no distributional assumptions (valid for skewed, heavy-tailed,
%     or bimodal distributions common in heterogeneous tumors).
%   - Is robust to outliers (uses ranks, not raw values), which is
%     important because IVIM fitting can produce extreme outlier values
%     for D* in poorly-perfused regions.
%   - Has only slightly less statistical power than the t-test when data
%     IS normally distributed (asymptotic relative efficiency = 0.955),
%     so there is minimal cost to using it as the default.
%   - Is valid for small samples (n >= 5 per group), unlike the t-test
%     which requires n >= ~30 for the Central Limit Theorem to compensate
%     for non-normality.

    if nargin < 3
        test_type = 'ranksum';
    end

    % Default to NaN: returned when there is insufficient data or the
    % test cannot be computed.  NaN (rather than 1.0 or 0) is the safest
    % sentinel because downstream code that checks p < alpha will
    % correctly treat NaN as "not significant" (NaN < 0.05 is false),
    % preventing false discoveries from missing data.
    p = NaN;

    % Ensure data and groups are column vectors for consistent indexing
    data = data(:);
    groups = groups(:);

    % Remove NaN/missing entries.  NaN values in DWI data arise from
    % failed IVIM fits (non-convergence), voxels outside the FOV, or
    % missing clinical data.  These must be excluded pairwise — a
    % patient with a valid group label but NaN parameter value cannot
    % contribute to the test.
    % Note: isnan() is undefined for categorical arrays in MATLAB, so
    % we use ismissing() which handles both categorical and numeric types.
    if iscategorical(groups)
        valid_idx = ~isnan(data) & ~ismissing(groups);
    else
        valid_idx = ~isnan(data) & ~isnan(groups);
    end
    data = data(valid_idx);
    groups = groups(valid_idx);

    unique_groups = unique(groups);

    switch lower(test_type)
        case 'ranksum'
            % Ranksum (Wilcoxon rank-sum / Mann-Whitney U) expects exactly
            % two groups.  In the pancreatic DWI context, typical group
            % comparisons include:
            %   - Local failure vs. no local failure (treatment outcome)
            %   - High vs. low baseline ADC (median split)
            %   - Responders vs. non-responders (based on ADC change > 20%)
            if length(unique_groups) == 2
                g1 = unique_groups(1);
                g2 = unique_groups(2);

                y1 = data(groups == g1);
                y2 = data(groups == g2);

                % Require at least 5 observations per group.  With n < 5,
                % the rank-sum test has very few possible permutations,
                % making exact p-values coarse (only a few discrete values
                % are possible) and the test clinically uninterpretable.
                % The n >= 5 threshold aligns with common recommendations
                % for the minimum sample size for non-parametric tests in
                % biomedical research.
                if length(y1) >= 5 && length(y2) >= 5
                    p = ranksum(y1, y2);
                end
            end
        otherwise
            error('Unsupported test_type: %s', test_type);
    end
end
