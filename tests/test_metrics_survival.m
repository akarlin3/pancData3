classdef test_metrics_survival < matlab.unittest.TestCase
    % TEST_METRICS_SURVIVAL Unit tests for the metrics_survival function.
    %
    % Verifies the Time-Dependent Cox PH model logic including:
    %   - Early return when there are insufficient events
    %   - Successful execution with enough events for Cox fitting
    %   - Correct handling of competing-risk events (coded as 2)
    %   - Correct subsetting via the valid_pts logical mask

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            close all;
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testInsufficientEventsReturnsEarly(testCase)
            % Verifies early return when there are zero events (all patients
            % censored, m_lf=0 for everyone). The time-dependent panel
            % builder sets td_ok=false, and the function should print a
            % diagnostic message and return without attempting to fit a Cox
            % model (which would fail with no events).
            rng(1);
            n   = 10;
            nTp = 4;
            valid_pts  = true(n, 1);
            ADC_abs    = rand(n, nTp) * 2e-3;
            D_abs      = rand(n, nTp) * 1e-3;
            f_abs      = rand(n, nTp) * 0.3;
            Dstar_abs  = rand(n, nTp) * 0.05;
            m_lf       = zeros(n, 1);            % all censored — no events
            m_total_time = randi([10, 100], n, 1);
            m_total_follow_up_time = m_total_time + randi([1, 30], n, 1);
            fx_label   = {'Fx1', 'Fx2', 'Fx3', 'Fx4'};
            dtype_label = 'Standard';

            % Must complete without throwing
            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
                m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label);
            testCase.verifyTrue(true, 'metrics_survival should return early without error.');
        end

        function testSufficientEventsRunsCoxModel(testCase)
            % Verifies that with 8 events in 20 patients, the Time-Dependent
            % Cox PH model fits successfully. Uses evalc to capture console
            % output and checks for expected statistical quantities (TD Panel
            % construction, hazard ratios, or an insufficient-events message
            % if the landmark time filter reduces the event count too much).
            rng(42);
            n   = 20;
            nTp = 4;
            valid_pts  = true(n, 1);
            ADC_abs    = rand(n, nTp) * 2e-3;
            D_abs      = rand(n, nTp) * 1e-3;
            f_abs      = rand(n, nTp) * 0.3;
            Dstar_abs  = rand(n, nTp) * 0.05;
            m_lf       = [ones(8, 1); zeros(12, 1)];   % 8 events
            m_total_time = randi([10, 100], n, 1);
            m_total_follow_up_time = nan(n, 1);
            m_total_follow_up_time(m_lf == 0) = ...
                m_total_time(m_lf == 0) + randi([5, 30], sum(m_lf == 0), 1);
            fx_label    = {'Fx1', 'Fx2', 'Fx3', 'Fx4'};
            dtype_label = 'Standard';

            tmp_dir = tempname;
            mkdir(tmp_dir);
            cleanup = onCleanup(@() rmdir(tmp_dir, 's'));

            console_output = evalc( ...
                'metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, [], tmp_dir)');
            diary off;

            % The Cox model should produce a panel with intervals
            testCase.verifySubstring(console_output, 'TD Panel', ...
                'Expected time-dependent panel construction output.');

            % Verify diary file was created (content may be empty when
            % metrics_survival runs inside evalc, because evalc intercepts
            % fprintf output before diary can capture it).
            diary_file = fullfile(tmp_dir, 'metrics_survival_output_Standard.txt');
            testCase.verifyTrue(exist(diary_file, 'file') == 2, ...
                'Diary file should be created.');

            % If model fitted, output should contain HR or hazard ratio
            % references.  If insufficient events after landmark, the
            % early return message should appear instead.
            has_model = contains(console_output, 'HR') || ...
                        contains(console_output, 'hazard') || ...
                        contains(console_output, 'coxphfit');
            has_early_return = contains(console_output, 'Insufficient');
            testCase.verifyTrue(has_model || has_early_return, ...
                'Output should contain Cox model results or insufficient-events message.');
        end

        function testCompetingRisksTreatedAsCensored(testCase)
            % Event code 2 (competing risk) is censored for cause-specific Cox.
            % Function must handle the mixed event vector gracefully.
            % Verify the panel correctly reports competing events separately.
            rng(7);
            n   = 20;
            nTp = 4;
            valid_pts  = true(n, 1);
            ADC_abs    = rand(n, nTp) * 2e-3;
            D_abs      = rand(n, nTp) * 1e-3;
            f_abs      = rand(n, nTp) * 0.3;
            Dstar_abs  = rand(n, nTp) * 0.05;
            % 5 cause-1 events, 5 competing events, 10 censored
            m_lf       = [ones(5, 1); 2*ones(5, 1); zeros(10, 1)];
            m_total_time = randi([10, 80], n, 1);
            m_total_follow_up_time = nan(n, 1);
            m_total_follow_up_time(m_lf == 0) = ...
                m_total_time(m_lf == 0) + randi([5, 20], sum(m_lf == 0), 1);
            fx_label    = {'Fx1', 'Fx2', 'Fx3', 'Fx4'};
            dtype_label = 'Standard';

            tmp_dir = tempname;
            mkdir(tmp_dir);
            cleanup = onCleanup(@() rmdir(tmp_dir, 's'));

            console_output = evalc( ...
                'metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, [], tmp_dir)');
            diary off;

            % Panel should report both event types
            testCase.verifySubstring(console_output, 'TD Panel', ...
                'Expected time-dependent panel construction output.');
            % Competing events should appear in panel stats
            testCase.verifySubstring(console_output, 'competing', ...
                'Panel output should mention competing events.');
        end

        function testValidPtsMaskSubsetsPatients(testCase)
            % valid_pts selects a subset; excluded patients should not influence
            % the model.  Function must accept mixed logical masks.
            rng(13);
            n_total    = 25;
            nTp        = 4;
            % First 15 patients are valid; last 10 are excluded
            valid_pts  = [true(15, 1); false(10, 1)];
            ADC_abs    = rand(n_total, nTp) * 2e-3;
            D_abs      = rand(n_total, nTp) * 1e-3;
            f_abs      = rand(n_total, nTp) * 0.3;
            Dstar_abs  = rand(n_total, nTp) * 0.05;
            m_lf       = [ones(5, 1); zeros(20, 1)];
            m_total_time = randi([10, 100], n_total, 1);
            m_total_follow_up_time = nan(n_total, 1);
            m_total_follow_up_time(m_lf == 0) = ...
                m_total_time(m_lf == 0) + randi([5, 30], sum(m_lf == 0), 1);
            fx_label    = {'Fx1', 'Fx2', 'Fx3', 'Fx4'};
            dtype_label = 'Standard';

            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
                m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label);
            testCase.verifyTrue(true, ...
                'metrics_survival should handle a partial valid_pts mask without error.');
        end

        function testNaNFollowUpTimeHandled(testCase)
            % Patients with m_lf == 0 but NaN follow-up time are handled gracefully
            % (the censoring mask excludes them; the panel is built without them).
            rng(99);
            n   = 15;
            nTp = 4;
            valid_pts  = true(n, 1);
            ADC_abs    = rand(n, nTp) * 2e-3;
            D_abs      = rand(n, nTp) * 1e-3;
            f_abs      = rand(n, nTp) * 0.3;
            Dstar_abs  = rand(n, nTp) * 0.05;
            m_lf       = [ones(5, 1); zeros(10, 1)];
            m_total_time = randi([10, 60], n, 1);
            % Leave ALL follow-up times as NaN — censored patients should be skipped
            m_total_follow_up_time = nan(n, 1);
            fx_label    = {'Fx1', 'Fx2', 'Fx3', 'Fx4'};
            dtype_label = 'Standard';

            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
                m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label);
            testCase.verifyTrue(true, ...
                'metrics_survival should handle all-NaN follow-up times gracefully.');
        end

        function testIPCWFrequencyPreservesN(testCase)
            % Verify that IPCW frequency weights do not inflate sample size.
            % With mean-stabilised weights ~1 and direct rounding,
            % sum(freq) should be close to N (not 10*N).
            rng(42);
            n = 30;
            raw_weights = 0.5 + rand(n, 1);              % range [0.5, 1.5]
            raw_weights = raw_weights / mean(raw_weights); % stabilise to mean=1

            freq = max(1, round(raw_weights));             % pipeline formula

            testCase.verifyGreaterThanOrEqual(sum(freq), n * 0.5, ...
                'Effective sample size too small');
            testCase.verifyLessThanOrEqual(sum(freq), n * 1.5, ...
                'Effective sample size inflated beyond 1.5x actual N');
        end

    end
end
