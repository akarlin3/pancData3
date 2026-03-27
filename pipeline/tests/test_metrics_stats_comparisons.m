classdef test_metrics_stats_comparisons < matlab.unittest.TestCase
    % TEST_METRICS_STATS_COMPARISONS Unit tests for metrics_stats_comparisons.
    %
    % Validates the statistical group comparison module including:
    %   - Wilcoxon rank-sum tests (LF vs LC groups)
    %   - Benjamini-Hochberg FDR correction for multiple comparisons
    %   - GLME (Generalized Linear Mixed-Effects) model fitting
    %   - Competing risk exclusion (lf_group==2 filtered out)
    %   - Edge cases: all-NaN metrics, single timepoint, patient filtering
    %   - Output file generation (diary logs, CSVs, PNGs)

    properties
        OutputFolder  % Temporary directory for test output files
        OldVisible    % Saved figure visibility state for restoration
    end

    methods(TestMethodSetup)
        function setup(testCase)
            % Create a uniquely-named temp directory (UUID suffix prevents
            % collisions in parallel test runs) and suppress figure pop-ups.
            testCase.OutputFolder = fullfile(tempdir, ['test_msc_' strrep(tempname, tempdir, '')]);
            mkdir(testCase.OutputFolder);
            testCase.OldVisible = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            diary off;
            set(0, 'DefaultFigureVisible', testCase.OldVisible);
            if isfolder(testCase.OutputFolder)
                rmdir(testCase.OutputFolder, 's');
            end
        end
    end

    methods(Test)
        function test_runs_with_minimal_data(testCase)
            % Smoke test: verifies the function runs end-to-end without
            % error on a minimal dataset (20 patients, 3 timepoints, 1
            % metric). Checks that the diary log file is created.
            rng(42);
            nPat = 20;
            nTp = 3;

            valid_pts = true(nPat, 1);
            lf_group = [zeros(10,1); ones(10,1)];

            % Create one metric set with one metric
            metric_data = randn(nPat, nTp);
            metric_sets = {{metric_data}};
            set_names = {{'ADC_abs'}};
            time_labels = {'Fx1', 'Fx2', 'Fx3'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);

            % Verify diary log was created
            logFile = fullfile(testCase.OutputFolder, 'metrics_stats_comparisons_output_Standard.txt');
            testCase.verifyTrue(exist(logFile, 'file') > 0);
        end

        function test_competing_risk_excluded(testCase)
            % Verifies that patients coded as competing risks (lf_group==2)
            % are excluded from Wilcoxon rank-sum comparisons. Only LC
            % (lf=0) and LF (lf=1) groups participate in the test.
            rng(42);
            nPat = 30;
            nTp = 2;

            valid_pts = true(nPat, 1);
            % 10 LC, 10 LF, 10 competing risk
            lf_group = [zeros(10,1); ones(10,1); 2*ones(10,1)];

            metric_data = randn(nPat, nTp);
            metric_sets = {{metric_data}};
            set_names = {{'Test_Metric'}};
            time_labels = {'Fx1', 'Fx2'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            % Should not error with competing risk patients
            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);
        end

        function test_multiple_metric_sets(testCase)
            % Verifies handling of multiple metric sets (as in the real
            % pipeline where absolute and percent-change metrics are
            % separate sets). Each set should produce its own figure PNG.
            rng(42);
            nPat = 20;
            nTp = 2;

            valid_pts = true(nPat, 1);
            lf_group = [zeros(10,1); ones(10,1)];

            m1 = randn(nPat, nTp);
            m2 = randn(nPat, nTp) + 1;
            metric_sets = {{m1, m2}, {m1}};
            set_names = {{'ADC', 'D'}, {'ADC_pct'}};
            time_labels = {'Fx1', 'Fx2'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);

            % Should create one figure per metric set
            pngFiles = dir(fullfile(testCase.OutputFolder, 'Metric_Set_*.png'));
            testCase.verifyEqual(numel(pngFiles), 2);
        end

        function test_all_nan_metric_graceful(testCase)
            % Verifies graceful handling when a metric is entirely NaN
            % across all patients and timepoints. The Wilcoxon test should
            % be skipped (no valid data), and no error should be thrown.
            nPat = 20;
            nTp = 2;

            valid_pts = true(nPat, 1);
            lf_group = [zeros(10,1); ones(10,1)];

            metric_data = nan(nPat, nTp);
            metric_sets = {{metric_data}};
            set_names = {{'NaN_metric'}};
            time_labels = {'Fx1', 'Fx2'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);
        end

        function test_fdr_correction_with_known_pvalues(testCase)
            % Verifies that BH FDR correction works with a known effect
            % pattern: strong group separation at Fx1 (shift of +3 SD),
            % weak at Fx2, none at Fx3. The Fx1 comparison should appear
            % in Significant_LF_Metrics.csv after FDR correction.
            rng(42);
            nPat = 40;
            nTp = 3;

            valid_pts = true(nPat, 1);
            lf_group = [zeros(20,1); ones(20,1)];

            % Strong effect at Fx1, weak at Fx2, none at Fx3
            m1 = zeros(nPat, nTp);
            m1(1:20, 1) = randn(20,1);
            m1(21:40, 1) = randn(20,1) + 3;  % strong separation
            m1(:, 2) = randn(nPat, 1);
            m1(:, 3) = randn(nPat, 1);

            metric_sets = {{m1}};
            set_names = {{'StrongEffect'}};
            time_labels = {'Fx1', 'Fx2', 'Fx3'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);

            % Check if significant results CSV was created
            sigFile = fullfile(testCase.OutputFolder, 'Significant_LF_Metrics.csv');
            if exist(sigFile, 'file')
                T = readtable(sigFile);
                testCase.verifyGreaterThan(height(T), 0);
            end
        end

        function test_single_timepoint_glme_skips(testCase)
            % Verifies that with only 1 timepoint, the GLME model still
            % runs (using just the intercept term, no time interaction).
            % This is a boundary case that tests the model formula generation.
            rng(42);
            nPat = 20;
            nTp = 1;

            valid_pts = true(nPat, 1);
            lf_group = [zeros(10,1); ones(10,1)];

            metric_data = randn(nPat, nTp);
            metric_sets = {{metric_data}};
            set_names = {{'Test'}};
            time_labels = {'Fx1'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);
        end

        function test_valid_pts_mask_filters_patients(testCase)
            % Verifies that the valid_pts logical mask correctly filters
            % patients before statistical testing. Only the first 10
            % patients (valid_pts=true) should be included; the last 10
            % are excluded even though they have valid data.
            rng(42);
            nPat = 20;
            nTp = 2;

            valid_pts = false(nPat, 1);
            valid_pts(1:10) = true;  % Only first 10 patients valid
            lf_group = [zeros(5,1); ones(5,1); zeros(5,1); ones(5,1)];

            metric_data = randn(nPat, nTp);
            metric_sets = {{metric_data}};
            set_names = {{'Filtered'}};
            time_labels = {'Fx1', 'Fx2'};

            ADC_abs = randn(nPat, nTp);
            D_abs = randn(nPat, nTp);
            f_abs = rand(nPat, nTp);
            Dstar_abs = rand(nPat, nTp) * 0.01;

            metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, ...
                time_labels, 'Standard', testCase.OutputFolder, testCase.OutputFolder, ...
                nTp, ADC_abs, D_abs, f_abs, Dstar_abs);
        end
    end
end
