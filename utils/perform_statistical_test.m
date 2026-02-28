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

    if nargin < 3
        test_type = 'ranksum';
    end

    p = NaN;

    % Ensure data and groups are column vectors
    data = data(:);
    groups = groups(:);

    % Remove NaNs
    valid_idx = ~isnan(data) & ~isnan(groups);
    data = data(valid_idx);
    groups = groups(valid_idx);

    unique_groups = unique(groups);

    switch lower(test_type)
        case 'ranksum'
            % Ranksum expects exactly two groups
            if length(unique_groups) == 2
                g1 = unique_groups(1);
                g2 = unique_groups(2);

                y1 = data(groups == g1);
                y2 = data(groups == g2);

                % Check if we have enough data in both groups to run Wilcoxon rank-sum test
                if length(y1) > 0 && length(y2) > 0 && (length(y1) + length(y2) > 2)
                    p = ranksum(y1, y2);
                end
            end
        otherwise
            error('Unsupported test_type: %s', test_type);
    end
end
