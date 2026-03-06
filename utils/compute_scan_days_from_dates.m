function scan_days = compute_scan_days_from_dates(fx_dates)
% COMPUTE_SCAN_DAYS_FROM_DATES  Convert DICOM StudyDate strings to scan-day offsets.
%
%   scan_days = compute_scan_days_from_dates(fx_dates)
%
%   Inputs:
%     fx_dates - Cell matrix (patients x fractions) of DICOM StudyDate
%                strings in 'YYYYMMDD' format.  Empty cells are treated as
%                missing fractions.
%
%   Outputs:
%     scan_days - 1 x nFx vector of median days-since-first-scan across
%                 patients.  Returns [] if fewer than 2 fractions have valid
%                 dates or if parsing fails.
%
% --- Analytical Rationale ---
% In adaptive pancreatic radiotherapy, DWI scans are acquired at specific
% treatment fractions (e.g., simulation, fraction 5, fraction 10, fraction 25).
% However, treatment schedules vary across patients due to holidays, machine
% downtime, or clinical holds. The actual calendar dates of "fraction 5" can
% differ by days or weeks between patients.
%
% For survival analysis and time-dependent modeling (Cox regression,
% build_td_panel), the time axis must be in continuous days rather than
% discrete fraction indices. Using fraction numbers (1, 2, 3, ...) as the
% time axis would incorrectly assume equal spacing, distorting hazard rate
% estimates and dose-response temporal correlations.
%
% This function computes a consensus (median) scan-day timeline across the
% cohort, which serves as the shared time axis for longitudinal feature
% panels. The median is preferred over the mean because it is robust to
% outlier patients with unusually delayed or accelerated treatment schedules.

scan_days = [];
if isempty(fx_dates)
    return;
end

[nPat, nFx] = size(fx_dates);

% --- DICOM Date Parsing ---
% DICOM StudyDate uses the 'YYYYMMDD' format (e.g., '20250115' for
% January 15, 2025) per the DICOM PS3.5 standard. Convert to MATLAB's
% datenum for arithmetic (days between scans).
% Invalid or missing dates are left as NaN to propagate gracefully.
date_nums = nan(nPat, nFx);
for j = 1:nPat
    for k = 1:nFx
        ds = fx_dates{j, k};
        if ~isempty(ds) && ischar(ds) && length(ds) == 8
            try
                date_nums(j, k) = datenum(ds, 'yyyymmdd');
            catch
                % Malformed date string (e.g., '00000000' from anonymized
                % DICOM headers): leave as NaN rather than crashing.
            end
        end
    end
end

% --- Per-Patient Baseline Normalization ---
% Normalize each patient's dates relative to a common reference fraction
% (the fraction column with the most data). Only patients who have the
% reference fraction contribute to the consensus timeline; patients
% missing it are left as NaN rather than using a misaligned fallback
% (normalizing to their own earliest scan puts them on a different time
% scale, which contaminates the per-column median and can break
% monotonicity).
fx_counts = sum(~isnan(date_nums), 1);
[~, ref_col] = max(fx_counts);
days_from_baseline = nan(nPat, nFx);
for j = 1:nPat
    if ~isnan(date_nums(j, ref_col))
        days_from_baseline(j, :) = date_nums(j, :) - date_nums(j, ref_col);
    end
end

% --- Cohort-Level Consensus Timeline ---
% Take the median across patients for each fraction column. The median
% provides a robust central estimate of when each fraction's scan
% typically occurs, tolerating patients whose treatment was delayed
% (e.g., toxicity-related breaks) without distorting the timeline.
median_days = nan(1, nFx);
for k = 1:nFx
    col = days_from_baseline(:, k);
    col = col(~isnan(col));
    if ~isempty(col)
        median_days(k) = median(col);
    end
end

% --- Positional Correspondence Preservation ---
% Return the full vector preserving positional correspondence with
% feature array columns. Filtering out NaN entries here would shift
% the mapping (e.g., fraction 4's features paired with fraction 3's
% scan day), so we keep NaN for missing fractions and let the caller
% handle them. This 1:1 mapping between scan_days indices and feature
% columns is a critical invariant for build_td_panel.
valid_fx = ~isnan(median_days);
if sum(valid_fx) < 2
    % Need at least 2 valid timepoints to define a longitudinal trajectory.
    % A single timepoint provides no temporal information.
    return;
end
scan_days = median_days;

% --- Monotonicity Validation ---
% Ensure the valid (non-NaN) entries are strictly increasing. This is a
% physical constraint: later fractions must occur on later calendar dates.
% Non-monotonic scan days indicate data corruption (e.g., mislabeled
% fraction folders, incorrect DICOM dates). build_td_panel requires
% strictly increasing time values for interpolation and panel construction,
% so we fall back to default spacing rather than propagating corrupt timing.
valid_days = median_days(valid_fx);
if any(diff(valid_days) <= 0)
    warning('compute_scan_days_from_dates:notIncreasing', ...
        'Computed scan days [%s] are not strictly increasing. Falling back to defaults.', ...
        num2str(median_days, '%.1f '));
    scan_days = [];
end
end
