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

scan_days = [];
if isempty(fx_dates)
    return;
end

[nPat, nFx] = size(fx_dates);

% Parse DICOM date strings ('YYYYMMDD') into MATLAB datenums
date_nums = nan(nPat, nFx);
for j = 1:nPat
    for k = 1:nFx
        ds = fx_dates{j, k};
        if ~isempty(ds) && ischar(ds) && length(ds) == 8
            try
                date_nums(j, k) = datenum(ds, 'yyyymmdd');
            catch
                % leave as NaN
            end
        end
    end
end

% Compute per-patient days since that patient's first scan
days_from_baseline = nan(nPat, nFx);
for j = 1:nPat
    valid = ~isnan(date_nums(j, :));
    if any(valid)
        days_from_baseline(j, :) = date_nums(j, :) - min(date_nums(j, valid));
    end
end

% Take the median across patients for each fraction
median_days = nan(1, nFx);
for k = 1:nFx
    col = days_from_baseline(:, k);
    col = col(~isnan(col));
    if ~isempty(col)
        median_days(k) = median(col);
    end
end

% Keep only fractions with valid medians
valid_fx = ~isnan(median_days);
if sum(valid_fx) < 2
    return;
end
scan_days = median_days(valid_fx);

% Ensure strictly increasing (required by build_td_panel)
if any(diff(scan_days) <= 0)
    warning('compute_scan_days_from_dates:notIncreasing', ...
        'Computed scan days are not strictly increasing. Falling back to defaults.');
    scan_days = [];
end
end
