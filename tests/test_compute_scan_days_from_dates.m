function tests = test_compute_scan_days_from_dates
% TEST_COMPUTE_SCAN_DAYS_FROM_DATES  Unit tests for compute_scan_days_from_dates.
%
% Tests the utility that converts a cell array of DICOM acquisition date
% strings (patients x fractions) into a consensus scan-day vector by
% computing per-patient day offsets relative to a reference fraction and
% taking the median across patients.
%
% Coverage:
%   - Basic conversion with uniform patient schedules
%   - Missing fractions (empty date strings) handled via partial medians
%   - Empty input returns empty output
%   - Single-fraction input returns empty (need >= 2 fractions)
%   - Varying patient schedules resolved via median
%   - Patients missing the reference fraction excluded from median
    tests = functiontests(localfunctions);
end

function test_basic_conversion(testCase)
    % Verifies the simplest case: 3 patients with identical schedules.
    % All patients scanned on the same dates, so the median offset at each
    % fraction should exactly match the known day offsets [0, 7, 90].
    base = datenum('20240101', 'yyyymmdd');
    fx_dates = {datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd')};
    scan_days = compute_scan_days_from_dates(fx_dates);
    testCase.verifyEqual(scan_days, [0, 7, 90]);
end

function test_missing_fractions(testCase)
    % Verifies that a missing date (empty string) for one patient at one
    % fraction does not crash the function.  The missing entry is excluded
    % from the median; with only patient 1 contributing, the median for
    % fraction 2 equals patient 1's offset of 5 days.
    base = datenum('20240101', 'yyyymmdd');
    fx_dates = {datestr(base, 'yyyymmdd'),   datestr(base+5, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   '',                           datestr(base+90, 'yyyymmdd')};
    scan_days = compute_scan_days_from_dates(fx_dates);
    % Fraction 2 has only 1 patient contributing median = 5
    testCase.verifyEqual(scan_days, [0, 5, 90]);
end

function test_empty_input(testCase)
    % An empty cell array should return an empty result without error.
    scan_days = compute_scan_days_from_dates({});
    testCase.verifyTrue(isempty(scan_days));
end

function test_single_fraction(testCase)
    % A single-column date cell (one fraction per patient) cannot produce
    % meaningful day offsets because there is no second timepoint to
    % compute a delta against.  The function should return empty.
    fx_dates = {'20240101'; '20240101'};
    scan_days = compute_scan_days_from_dates(fx_dates);
    testCase.verifyTrue(isempty(scan_days));
end

function test_varying_patient_dates_median(testCase)
    % Two patients with slightly different schedules; median should be used
    fx_dates = {'20240101', '20240106'; ...   % offsets [0, 5]
                '20240101', '20240108'};      % offsets [0, 7]
    scan_days = compute_scan_days_from_dates(fx_dates);
    testCase.verifyEqual(scan_days, [0, 6]);  % median of [5, 7] = 6
end

function test_missing_ref_col_excluded(testCase)
    % Patients missing the reference fraction should be excluded from the
    % median rather than using a misaligned fallback that can break
    % monotonicity.  Here ref_col=2 (3 of 4 patients have Fx2).
    % Patient 4 has only Fx3.  Old fallback would normalize patient 4 to
    % its earliest scan (Fx3=day0), contaminating the Fx3 column median.
    base = datenum('20240101', 'yyyymmdd');
    fx_dates = {datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                '',                          '',                           datestr(base+90, 'yyyymmdd')};
    scan_days = compute_scan_days_from_dates(fx_dates);
    % Patient 4 excluded (missing ref_col=1). Median from patients 1-3.
    testCase.verifyEqual(scan_days, [0, 7, 90]);
end
