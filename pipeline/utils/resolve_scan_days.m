function td_scan_days = resolve_scan_days(summary_metrics, config_struct)
% RESOLVE_SCAN_DAYS  Three-level scan day resolution for survival analysis.
%
%   td_scan_days = resolve_scan_days(summary_metrics, config_struct)
%
%   Resolves the actual calendar days (relative to treatment start) when
%   each fraction's MRI was acquired.  Uses a three-level fallback strategy
%   in decreasing order of data quality:
%
%     1. DICOM StudyDate headers (most accurate — actual scanner timestamps)
%     2. config.json td_scan_days (user-specified from clinical records)
%     3. Empty [] — lets metrics_survival use built-in defaults with an
%        immortal-time-bias warning
%
%   Correct scan day assignment is critical for time-dependent Cox models:
%   using assumed scan days when actual dates differ introduces immortal
%   time bias, systematically distorting hazard ratio estimates.
%
%   Inputs:
%     summary_metrics - Pipeline summary metrics struct.  If it has an
%                       fx_dates field (cell array of DICOM date strings),
%                       these are converted to scan days via
%                       compute_scan_days_from_dates.
%     config_struct   - Pipeline configuration struct.  May contain
%                       .td_scan_days (numeric vector) as a fallback.
%
%   Outputs:
%     td_scan_days - Numeric vector of scan days relative to treatment
%                    start, or [] if no source could provide valid days.

    td_scan_days = [];

    % Level 1: DICOM StudyDate headers
    if isfield(summary_metrics, 'fx_dates') && ~isempty(summary_metrics.fx_dates)
        n_dates = sum(~cellfun('isempty', summary_metrics.fx_dates(:)));
        fprintf('      💡 Found %d DICOM StudyDate entries across cohort.\n', n_dates);
        td_scan_days = compute_scan_days_from_dates(summary_metrics.fx_dates);
        if isempty(td_scan_days)
            fprintf('      ⚠️  Could not derive scan days from DICOM dates (insufficient valid dates or non-monotonic).\n');
        end
    else
        fprintf('      💡 No DICOM StudyDate data available (fx_dates empty).\n');
    end

    % Level 2: config.json td_scan_days
    if isempty(td_scan_days) && isfield(config_struct, 'td_scan_days') && ~isempty(config_struct.td_scan_days)
        td_scan_days = config_struct.td_scan_days;
        fprintf('      💡 Using td_scan_days from config.json.\n');
    end

    % Level 3: empty — metrics_survival will use its own defaults
    if isempty(td_scan_days)
        fprintf('      💡 Set "td_scan_days" in config.json to override defaults.\n');
    end
end
