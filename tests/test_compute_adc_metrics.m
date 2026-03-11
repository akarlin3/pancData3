classdef test_compute_adc_metrics < matlab.unittest.TestCase
    % TEST_COMPUTE_ADC_METRICS Unit tests for compute_adc_metrics.
    %
    % Validates ADC metric computation: GTV volume, mean/SD, sub-volume
    % delineation, high-ADC fraction, motion corruption flag, histogram,
    % and KS test logic.

    properties
        Config
        BinEdges
    end

    methods(TestMethodSetup)
        function setupConfig(testCase)
            testCase.Config = struct( ...
                'adc_thresh', 1e-3, ...
                'high_adc_thresh', 1.15e-3, ...
                'd_thresh', 1e-3, ...
                'f_thresh', 0.1, ...
                'dstar_thresh', 0.01, ...
                'core_method', 'adc_threshold', ...
                'fdm_parameter', 'adc', ...
                'fdm_thresh', 0.0004);
            testCase.BinEdges = linspace(0, 3e-3, 31);
        end
    end

    methods(Test)
        function test_empty_adc_returns_nan(testCase)
            % Empty ADC vector should return all-NaN struct.
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, [], [], [], [], ...
                [], 1, 100, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyTrue(isnan(result.adc_mean_val));
            testCase.verifyTrue(isnan(result.gtv_vol_val));
        end

        function test_gtv_volume_correct(testCase)
            % GTV volume = numel(adc_vec) * vox_vol, including NaN voxels.
            adc = [0.001; 0.002; NaN; 0.0015];
            vox_vol = 0.05;  % cm^3
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], vox_vol, 100, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyEqual(result.gtv_vol_val, 4 * 0.05, 'AbsTol', 1e-10);
        end

        function test_mean_adc_ignores_nan(testCase)
            % Mean should be computed over finite values only.
            adc = [0.001; NaN; 0.002];
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], 1, 100, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyEqual(result.adc_mean_val, 0.0015, 'AbsTol', 1e-10);
        end

        function test_sub_volume_with_threshold(testCase)
            % Voxels below adc_thresh should be in the core sub-volume.
            adc = [0.0005; 0.0008; 0.002; 0.003];
            vox_vol = 1;
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], vox_vol, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            % adc_threshold: voxels < 1e-3 are [0.0005, 0.0008] = 2 voxels
            testCase.verifyEqual(result.adc_sub_vol_val, 2, 'AbsTol', 1e-10);
        end

        function test_high_adc_sub_volume(testCase)
            % High-ADC sub-volume counts voxels above high_adc_thresh.
            adc = [0.0005; 0.0012; 0.002; 0.003];
            vox_vol = 0.1;
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], vox_vol, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            % Voxels > 1.15e-3: [0.0012, 0.002, 0.003] = 3 voxels
            testCase.verifyEqual(result.high_adc_sub_vol_val, 3 * 0.1, 'AbsTol', 1e-10);
        end

        function test_motion_corruption_flag(testCase)
            % Fraction of voxels above adc_max should be computed.
            adc = [0.001; 0.002; 0.004; 0.005];  % adc_max = 3e-3
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], 1, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            % 2 out of 4 above 3e-3
            testCase.verifyEqual(result.fx_corrupted_val, 0.5, 'AbsTol', 1e-10);
        end

        function test_ks_test_skipped_at_baseline(testCase)
            % At k=1 (baseline), KS test should not be performed.
            adc = [0.001; 0.002; 0.0015];
            baseline = [0.001; 0.002; 0.0015];
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                baseline, 1, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyTrue(isnan(result.ks_stat_adc));
            testCase.verifyTrue(isnan(result.ks_pval_adc));
        end

        function test_histogram_laplace_smoothed(testCase)
            % Histogram should have Laplace smoothing (no zero bins when data present).
            adc = [0.0005; 0.001; 0.0015; 0.002];
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], 1, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyEqual(numel(result.adc_histogram), numel(testCase.BinEdges) - 1);
            % All bins should be > 0 due to Laplace smoothing
            testCase.verifyTrue(all(result.adc_histogram > 0));
        end

        function test_sd_computed(testCase)
            % Standard deviation should be computed for non-empty vectors.
            adc = [0.001; 0.002; 0.003];
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(testCase.Config, adc, [], [], [], ...
                [], 1, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyFalse(isnan(result.adc_sd_val));
            testCase.verifyGreaterThan(result.adc_sd_val, 0);
        end

        function test_fdm_fractions_at_baseline_nan(testCase)
            % fDM fractions should remain NaN at baseline (k=1).
            adc = [0.001; 0.002];
            cfg = testCase.Config;
            cfg.core_method = 'fdm';
            core_opts = struct('timepoint_index', 1);
            result = compute_adc_metrics(cfg, adc, [], [], [], ...
                [], 1, 1, testCase.BinEdges, 1.15e-3, 3e-3, ...
                false, [], core_opts, 1, 1, [], 1);
            testCase.verifyTrue(isnan(result.fdm_responding_pc));
            testCase.verifyTrue(isnan(result.fdm_progressing_pc));
            testCase.verifyTrue(isnan(result.fdm_stable_pc));
        end
    end
end
