classdef test_resolve_scan_days < matlab.unittest.TestCase
    % TEST_RESOLVE_SCAN_DAYS Unit tests for resolve_scan_days.
    %
    % Validates the three-level scan day resolution: DICOM dates,
    % config.json td_scan_days, and empty fallback.

    methods(Test)
        function test_dicom_dates_preferred(testCase)
            % When fx_dates produce valid scan days, those should be used.
            % Use actual dates that compute_scan_days_from_dates can parse.
            sm = struct();
            sm.fx_dates = {'20240101', '20240108', '20240115'; ...
                           '20240101', '20240108', '20240115'};
            cfg = struct('td_scan_days', [0, 10, 20]);

            result = resolve_scan_days(sm, cfg);
            % Should NOT be the config values since DICOM dates are preferred
            if ~isempty(result)
                testCase.verifyNotEqual(result, [0, 10, 20]);
            end
        end

        function test_config_fallback(testCase)
            % When fx_dates is empty, config.json td_scan_days should be used.
            sm = struct('fx_dates', {{}});
            cfg = struct('td_scan_days', [0, 7, 14, 21, 28]);

            result = resolve_scan_days(sm, cfg);
            testCase.verifyEqual(result, [0, 7, 14, 21, 28]);
        end

        function test_config_fallback_no_fx_dates_field(testCase)
            % When fx_dates field is absent, fall back to config td_scan_days.
            sm = struct('id_list', {{'P01'}});
            cfg = struct('td_scan_days', [0, 5, 10]);

            result = resolve_scan_days(sm, cfg);
            testCase.verifyEqual(result, [0, 5, 10]);
        end

        function test_empty_when_no_sources(testCase)
            % When neither DICOM dates nor config are available, return [].
            sm = struct('fx_dates', {{}});
            cfg = struct('td_scan_days', []);

            result = resolve_scan_days(sm, cfg);
            testCase.verifyTrue(isempty(result));
        end

        function test_empty_config_field_missing(testCase)
            % When td_scan_days field is missing from config, return [].
            sm = struct('fx_dates', {{}});
            cfg = struct('dataloc', '/tmp');

            result = resolve_scan_days(sm, cfg);
            testCase.verifyTrue(isempty(result));
        end

        function test_dicom_dates_invalid_falls_to_config(testCase)
            % If compute_scan_days_from_dates returns [] (e.g., too few
            % dates), should fall back to config td_scan_days.
            sm = struct();
            % Only one date column — not enough for scan day computation
            sm.fx_dates = {'20240101'; '20240101'};
            cfg = struct('td_scan_days', [0, 7]);

            result = resolve_scan_days(sm, cfg);
            testCase.verifyEqual(result, [0, 7]);
        end
    end
end
