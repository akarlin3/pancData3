function p_str = format_p_value(p)
% FORMAT_P_VALUE Formats a p-value for display in plots.
% Prevents displaying 'p = 0.000' for very small values.
    if isnan(p)
        p_str = 'p = NaN';
    elseif p < 0.001
        p_str = 'p < 0.001';
    else
        p_str = sprintf('p = %.3f', p);
    end
end
