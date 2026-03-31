function stats = compute_median_followup(baseline)
% COMPUTE_MEDIAN_FOLLOWUP  Median follow-up time for the cohort.
%
%   stats = compute_median_followup(baseline)
%
%   Computes median follow-up time using two methods:
%     1. Simple median of observed follow-up times across all patients
%     2. Reverse Kaplan-Meier (standard reporting method): treats censoring
%        as the "event" and actual events as censored, giving an unbiased
%        estimate of potential follow-up time.
%
%   Inputs:
%     baseline - Struct from metrics_baseline / load_baseline_from_disk
%                with fields: m_lf, m_total_time, m_total_follow_up_time,
%                valid_pts
%
%   Outputs:
%     stats - Struct with fields:
%       .median_followup_days    - Simple median (days)
%       .median_followup_months  - Simple median (months, /30.44)
%       .rkm_median_days         - Reverse KM median (days)
%       .rkm_median_months       - Reverse KM median (months)
%       .n_patients              - Number of patients included
%       .n_events                - Number of LF events
%       .n_censored              - Number of censored patients
%       .n_competing             - Number of competing risk events
%       .range_months            - [min, max] follow-up in months
%
%   Example:
%     baseline = load_baseline_from_disk('metrics_baseline_results_Standard.mat');
%     stats = compute_median_followup(baseline);
%     fprintf('Median follow-up: %.1f months (reverse KM: %.1f months)\n', ...
%         stats.median_followup_months, stats.rkm_median_months);

    % =====================================================================
    % ASSEMBLE PER-PATIENT FOLLOW-UP TIME
    % =====================================================================
    % For each patient, the observed follow-up is:
    %   LF patients (lf==1): time to event (m_total_time)
    %   Censored (lf==0) or competing risk (lf==2): time to last contact
    %                                                (m_total_follow_up_time)
    lf  = baseline.m_lf;
    t_event   = baseline.m_total_time;
    t_followup = baseline.m_total_follow_up_time;

    % Restrict to valid patients (those with baseline imaging)
    if isfield(baseline, 'valid_pts') && ~isempty(baseline.valid_pts)
        idx = baseline.valid_pts;
        lf = lf(idx);
        t_event = t_event(idx);
        t_followup = t_followup(idx);
    end

    % Build unified time vector: event time for LF, follow-up time for censored
    obs_time = t_event;
    cens_mask = (lf == 0 | lf == 2) & ~isnan(t_followup);
    obs_time(cens_mask) = t_followup(cens_mask);

    % Remove patients with missing time data
    valid = ~isnan(obs_time) & obs_time > 0;
    obs_time = obs_time(valid);
    lf = lf(valid);

    n_patients = numel(obs_time);
    n_events   = sum(lf == 1);
    n_censored = sum(lf == 0);
    n_competing = sum(lf == 2);

    % =====================================================================
    % METHOD 1: SIMPLE MEDIAN
    % =====================================================================
    median_days = median(obs_time);

    % =====================================================================
    % METHOD 2: REVERSE KAPLAN-MEIER
    % =====================================================================
    % In the reverse KM, we flip the event indicator: censored patients
    % become "events" (they completed follow-up) and LF/competing-risk
    % patients become "censored" (follow-up was interrupted by the event).
    % This gives an unbiased estimate of potential follow-up.
    rkm_event = double(lf == 0 | lf == 2);  % "event" = reached end of follow-up
    rkm_median = reverse_km_median(obs_time, rkm_event);

    % =====================================================================
    % PACK OUTPUT
    % =====================================================================
    days_per_month = 30.44;

    stats.median_followup_days   = median_days;
    stats.median_followup_months = median_days / days_per_month;
    stats.rkm_median_days        = rkm_median;
    stats.rkm_median_months      = rkm_median / days_per_month;
    stats.n_patients             = n_patients;
    stats.n_events               = n_events;
    stats.n_censored             = n_censored;
    stats.n_competing            = n_competing;
    stats.range_months           = [min(obs_time), max(obs_time)] / days_per_month;

    % =====================================================================
    % DISPLAY SUMMARY
    % =====================================================================
    fprintf('\n=== MEDIAN FOLLOW-UP TIME ===\n');
    fprintf('  Cohort size:       %d patients\n', n_patients);
    fprintf('    LF events:       %d\n', n_events);
    fprintf('    Censored:        %d\n', n_censored);
    fprintf('    Competing risk:  %d\n', n_competing);
    fprintf('  Simple median:     %.0f days  (%.1f months)\n', ...
        median_days, stats.median_followup_months);
    fprintf('  Reverse KM median: %.0f days  (%.1f months)\n', ...
        rkm_median, stats.rkm_median_months);
    fprintf('  Range:             %.1f - %.1f months\n', ...
        stats.range_months(1), stats.range_months(2));
    fprintf('============================\n\n');
end

function med = reverse_km_median(times, events)
% REVERSE_KM_MEDIAN  Compute median from Kaplan-Meier survival curve.
%   times  - observation times (n x 1)
%   events - event indicator (1 = event, 0 = censored) (n x 1)

    [sorted_times, order] = sort(times);
    sorted_events = events(order);

    n = numel(times);
    at_risk = n;
    surv = 1.0;
    med = max(times);  % fallback if survival never drops below 0.5

    for i = 1:n
        if sorted_events(i) == 1
            surv = surv * (1 - 1 / at_risk);
            if surv <= 0.5
                med = sorted_times(i);
                return;
            end
        end
        at_risk = at_risk - 1;
    end
end
