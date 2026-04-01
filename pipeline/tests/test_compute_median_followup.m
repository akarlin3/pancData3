classdef test_compute_median_followup < matlab.unittest.TestCase
    % TEST_COMPUTE_MEDIAN_FOLLOWUP  Tests for compute_median_followup.
    %
    % Covers:
    %   - All-censored cohort (simple and reverse KM agree)
    %   - Mixed events and censored patients
    %   - Competing risk patients counted correctly
    %   - valid_pts subsetting
    %   - Output struct completeness

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testAllCensored(testCase)
            % When all patients are censored, simple median = reverse KM median
            baseline.m_lf = [0; 0; 0; 0; 0];
            baseline.m_total_time = [100; 200; 300; 400; 500];
            baseline.m_total_follow_up_time = [100; 200; 300; 400; 500];
            baseline.valid_pts = (1:5)';

            stats = compute_median_followup(baseline);

            testCase.verifyEqual(stats.n_patients, 5);
            testCase.verifyEqual(stats.n_events, 0);
            testCase.verifyEqual(stats.n_censored, 5);
            testCase.verifyEqual(stats.median_followup_days, 300);
            % Reverse KM: all are "events", so median = simple median
            testCase.verifyEqual(stats.rkm_median_days, 300);
        end

        function testMixedCohort(testCase)
            % 3 censored at 200, 400, 600 days; 2 LF events at 150, 350
            baseline.m_lf = [1; 0; 1; 0; 0];
            baseline.m_total_time = [150; NaN; 350; NaN; NaN];
            baseline.m_total_follow_up_time = [NaN; 200; NaN; 400; 600];
            baseline.valid_pts = (1:5)';

            stats = compute_median_followup(baseline);

            testCase.verifyEqual(stats.n_patients, 5);
            testCase.verifyEqual(stats.n_events, 2);
            testCase.verifyEqual(stats.n_censored, 3);
            % Observed times: [150, 200, 350, 400, 600] -> median = 350
            testCase.verifyEqual(stats.median_followup_days, 350);
        end

        function testCompetingRisk(testCase)
            % Competing risk patients should be counted separately
            baseline.m_lf = [0; 1; 2; 0; 2];
            baseline.m_total_time = [NaN; 180; 250; NaN; 300];
            baseline.m_total_follow_up_time = [500; NaN; 250; 400; 300];
            baseline.valid_pts = (1:5)';

            stats = compute_median_followup(baseline);

            testCase.verifyEqual(stats.n_competing, 2);
            testCase.verifyEqual(stats.n_events, 1);
            testCase.verifyEqual(stats.n_censored, 2);
        end

        function testValidPtsSubsetting(testCase)
            % Only patients in valid_pts should be included
            baseline.m_lf = [0; 0; 0; 0; 0];
            baseline.m_total_time = [100; 200; 300; 400; 500];
            baseline.m_total_follow_up_time = [100; 200; 300; 400; 500];
            baseline.valid_pts = [2; 4];  % Only patients 2 and 4

            stats = compute_median_followup(baseline);

            testCase.verifyEqual(stats.n_patients, 2);
            testCase.verifyEqual(stats.median_followup_days, 300);  % median(200, 400)
        end

        function testOutputFields(testCase)
            baseline.m_lf = [0; 1];
            baseline.m_total_time = [NaN; 90];
            baseline.m_total_follow_up_time = [365; NaN];
            baseline.valid_pts = [1; 2];

            stats = compute_median_followup(baseline);

            testCase.verifyTrue(isfield(stats, 'median_followup_days'));
            testCase.verifyTrue(isfield(stats, 'median_followup_months'));
            testCase.verifyTrue(isfield(stats, 'rkm_median_days'));
            testCase.verifyTrue(isfield(stats, 'rkm_median_months'));
            testCase.verifyTrue(isfield(stats, 'n_patients'));
            testCase.verifyTrue(isfield(stats, 'range_months'));
            testCase.verifyEqual(numel(stats.range_months), 2);
        end

        function testNoValidPtsField(testCase)
            % Should work without valid_pts (use all patients)
            baseline.m_lf = [0; 0; 0];
            baseline.m_total_time = [100; 200; 300];
            baseline.m_total_follow_up_time = [100; 200; 300];

            stats = compute_median_followup(baseline);

            testCase.verifyEqual(stats.n_patients, 3);
            testCase.verifyEqual(stats.median_followup_days, 200);
        end

        function testMonthConversion(testCase)
            baseline.m_lf = [0];
            baseline.m_total_time = [NaN];
            baseline.m_total_follow_up_time = [365];
            baseline.valid_pts = 1;

            stats = compute_median_followup(baseline);

            testCase.verifyEqual(stats.median_followup_months, 365 / 30.44, 'AbsTol', 0.01);
        end

    end
end
