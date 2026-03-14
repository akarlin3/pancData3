function p_str = format_p_value(p)
% FORMAT_P_VALUE Formats a p-value for display in plots and console output.
%
%   p_str = format_p_value(p)
%
%   Inputs:
%       p - Numeric p-value from a statistical test (e.g., Wilcoxon rank-sum).
%
%   Outputs:
%       p_str - Formatted string suitable for plot annotation or logging.
%
% --- Analytical Rationale ---
% In medical physics publications and clinical reports, p-values are
% conventionally displayed to 3 decimal places. However, displaying
% 'p = 0.000' for highly significant results is misleading -- it implies
% exact zero probability, which is mathematically impossible for continuous
% test statistics. The convention 'p < 0.001' is used instead, following
% AMA (American Medical Association) style guidelines for biomedical
% reporting.
%
% NaN p-values arise when the statistical test cannot be computed (e.g.,
% all values in one group are identical, or a group has insufficient
% sample size). Displaying 'p = NaN' explicitly flags these cases rather
% than silently omitting the annotation, which could be misinterpreted as
% an oversight.
    if isnan(p)
        p_str = 'p = NaN';
    elseif p < 0.001
        % Report as inequality per biomedical reporting convention.
        p_str = 'p < 0.001';
    else
        % Three decimal places provides sufficient precision for clinical
        % significance assessment (alpha = 0.05 threshold) without
        % implying spurious precision.
        p_str = sprintf('p = %.3f', p);
    end
end
