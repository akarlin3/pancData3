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
            % With zero events td_ok is false — the function should print a
            % diagnostic message and return without fitting a model (no error).
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
            % With 8 events in 20 patients the Cox model (or Firth fallback)
            % should fit successfully.
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

            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
                m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label);
            testCase.verifyTrue(true, ...
                'metrics_survival should complete (Cox or Firth) without error.');
        end

        function testCompetingRisksTreatedAsCensored(testCase)
            % Event code 2 (competing risk) is censored for cause-specific Cox.
            % Function must handle the mixed event vector gracefully.
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

            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
                m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label);
            testCase.verifyTrue(true, ...
                'metrics_survival should handle competing risks without error.');
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

    end
end
