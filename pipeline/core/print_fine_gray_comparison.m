function print_fine_gray_comparison(cox_results, fg_results, feat_names)
% PRINT_FINE_GRAY_COMPARISON  CSH HR vs Fine-Gray sHR comparison table.
fprintf('\n  --- CSH vs Fine-Gray Subdistribution Hazard Comparison ---\n');
fprintf('  %-10s %8s %8s %8s %8s  |  %8s %8s %8s %8s\n', ...
    'Feature', 'CSH_HR', 'CI_lo', 'CI_hi', 'p', 'sHR', 'CI_lo', 'CI_hi', 'p');
fprintf('  %s\n', repmat('-', 1, 86));

for fi = 1:numel(feat_names)
    % CSH results
    if cox_results.success
        b_csh = cox_results.coefficients(fi);
        se_csh = cox_results.se(fi);
        hr_csh = cox_results.hazard_ratios(fi);
        p_csh = cox_results.p_values(fi);
        ci_lo_csh = exp(b_csh - 1.96 * se_csh);
        ci_hi_csh = exp(b_csh + 1.96 * se_csh);
    else
        hr_csh = NaN; ci_lo_csh = NaN; ci_hi_csh = NaN; p_csh = NaN;
    end

    % Fine-Gray results
    b_fg = fg_results.coefficients(fi);
    se_fg = fg_results.se(fi);
    hr_fg = fg_results.hazard_ratios(fi);
    p_fg = fg_results.p_values(fi);
    ci_lo_fg = exp(b_fg - 1.96 * se_fg);
    ci_hi_fg = exp(b_fg + 1.96 * se_fg);

    fprintf('  %-10s %8.3f %8.3f %8.3f %8.4f  |  %8.3f %8.3f %8.3f %8.4f\n', ...
        feat_names{fi}, hr_csh, ci_lo_csh, ci_hi_csh, p_csh, ...
        hr_fg, ci_lo_fg, ci_hi_fg, p_fg);
end
end
