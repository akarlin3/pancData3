function tests = test_compute_scan_days_from_dates
% TEST_COMPUTE_SCAN_DAYS_FROM_DATES  Unit tests for compute_scan_days_from_dates
    tests = functiontests(localfunctions);
end

function test_basic_conversion(testCase)
    % 3 patients x 3 fractions with known offsets [0, 7, 90]
    base = datenum('20240101', 'yyyymmdd');
    fx_dates = {datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   datestr(base+7, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd')};
    scan_days = compute_scan_days_from_dates(fx_dates);
    testCase.verifyEqual(scan_days, [0, 7, 90]);
end

function test_missing_fractions(testCase)
    % Patient 2 is missing fraction 2
    base = datenum('20240101', 'yyyymmdd');
    fx_dates = {datestr(base, 'yyyymmdd'),   datestr(base+5, 'yyyymmdd'),  datestr(base+90, 'yyyymmdd'); ...
                datestr(base, 'yyyymmdd'),   '',                           datestr(base+90, 'yyyymmdd')};
    scan_days = compute_scan_days_from_dates(fx_dates);
    % Fraction 2 has only 1 patient contributing median = 5
    testCase.verifyEqual(scan_days, [0, 5, 90]);
end

function test_empty_input(testCase)
    scan_days = compute_scan_days_from_dates({});
    testCase.verifyTrue(isempty(scan_days));
end

function test_single_fraction(testCase)
    % Only one fraction: should return empty (need >= 2)
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
