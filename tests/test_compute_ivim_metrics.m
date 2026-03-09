classdef test_compute_ivim_metrics < matlab.unittest.TestCase
    % TEST_COMPUTE_IVIM_METRICS Unit tests for compute_ivim_metrics.
    %
    % Validates IVIM metric computation: failed-fit filtering, mean/SD,
    % sub-volume with independent vs unified methods, KS test, and
    % histogram generation.

    properties
        Config
        BinEdges
    end

    methods(TestMethodSetup)
        function setupConfig(testCase)
            testCase.Config = struct( ...
                'core_method', 'adc_threshold', ...
                'fdm_parameter', 'adc', ...
                'fdm_thresh', 0.0004);
            testCase.BinEdges = linspace(0, 3e-3, 31);
        end
    end

    methods(Test)
        function test_empty_d_returns_nan(testCase)
            % Empty D vector should return all-NaN struct.
            result = compute_ivim_metrics(testCase.Config, [], [], [], ...
                [], false(0,1), 1, 100, testCase.BinEdges, 1e-3, 0.1, 1);
            testCase.verifyTrue(isnan(result.d_mean_val));
            testCase.verifyTrue(isnan(result.f_mean_val));
            testCase.verifyTrue(isnan(result.dstar_mean_val));
        end

        function test_failed_fit_filtering(testCase)
            % Voxels with f=0 and D*=0 or D*=NaN should be NaN'd out.
            d_vec = [0.001; 0.002; 0.003];
            f_vec = [0.1; 0; 0.15];       % voxel 2: f=0
            dstar_vec = [0.01; 0; 0.02];  % voxel 2: D*=0 (failed fit)
            adc_mask = true(3,1);

            result = compute_ivim_metrics(testCase.Config, d_vec, f_vec, dstar_vec, ...
                [], adc_mask, 1, 1, testCase.BinEdges, 1e-3, 0.1, 1);

            % Voxel 2 should be NaN in the filtered vectors
            testCase.verifyTrue(isnan(result.f_vec(2)));
            testCase.verifyTrue(isnan(result.dstar_vec(2)));
            % Voxels 1 and 3 should be unchanged
            testCase.verifyEqual(result.f_vec(1), 0.1);
            testCase.verifyEqual(result.f_vec(3), 0.15);
        end

        function test_d_mean_correct(testCase)
            % D mean should be correctly computed.
            d_vec = [0.001; 0.002; 0.003];
            f_vec = [0.1; 0.2; 0.15];
            dstar_vec = [0.01; 0.02; 0.015];
            adc_mask = true(3,1);

            result = compute_ivim_metrics(testCase.Config, d_vec, f_vec, dstar_vec, ...
                [], adc_mask, 1, 1, testCase.BinEdges, 1e-3, 0.1, 1);
            testCase.verifyEqual(result.d_mean_val, 0.002, 'AbsTol', 1e-10);
        end

        function test_independent_threshold_f_sub_vol(testCase)
            % With independent thresholds (non-unified method), f sub-volume
            % should count voxels where f < f_thresh.
            d_vec = [0.001; 0.002; 0.003];
            f_vec = [0.05; 0.08; 0.15];  % 2 voxels below f_thresh=0.1
            dstar_vec = [0.01; 0.02; 0.015];
            adc_mask = true(3,1);
            vox_vol = 0.5;

            result = compute_ivim_metrics(testCase.Config, d_vec, f_vec, dstar_vec, ...
                [], adc_mask, vox_vol, 1, testCase.BinEdges, 1e-3, 0.1, 1);
            testCase.verifyEqual(result.f_sub_vol_val, 2 * 0.5, 'AbsTol', 1e-10);
        end

        function test_unified_method_uses_adc_mask(testCase)
            % Unified methods (percentile, spectral, fdm) should use the
            % ADC-derived mask for the f sub-volume instead of f_thresh.
            cfg = testCase.Config;
            cfg.core_method = 'percentile';

            d_vec = [0.001; 0.002; 0.003];
            f_vec = [0.05; 0.2; 0.15];
            dstar_vec = [0.01; 0.02; 0.015];
            % Only voxels 1 and 3 are in the core
            adc_mask = [true; false; true];
            vox_vol = 1;

            result = compute_ivim_metrics(cfg, d_vec, f_vec, dstar_vec, ...
                [], adc_mask, vox_vol, 1, testCase.BinEdges, 1e-3, 0.1, 1);
            % f sub-volume should be 2 voxels (mask-based, not threshold-based)
            testCase.verifyEqual(result.f_sub_vol_val, 2, 'AbsTol', 1e-10);
        end

        function test_ks_test_skipped_at_baseline(testCase)
            % At k=1 (baseline), KS test should not run.
            d_vec = [0.001; 0.002];
            f_vec = [0.1; 0.2];
            dstar_vec = [0.01; 0.02];
            d_baseline = [0.001; 0.002];
            adc_mask = true(2,1);

            result = compute_ivim_metrics(testCase.Config, d_vec, f_vec, dstar_vec, ...
                d_baseline, adc_mask, 1, 1, testCase.BinEdges, 1e-3, 0.1, 1);
            testCase.verifyTrue(isnan(result.ks_stat_d));
            testCase.verifyTrue(isnan(result.ks_pval_d));
        end

        function test_histogram_length(testCase)
            % Histogram output should match bin_edges - 1 in length.
            d_vec = [0.001; 0.002; 0.003];
            f_vec = [0.1; 0.2; 0.15];
            dstar_vec = [0.01; 0.02; 0.015];
            adc_mask = true(3,1);

            result = compute_ivim_metrics(testCase.Config, d_vec, f_vec, dstar_vec, ...
                [], adc_mask, 1, 1, testCase.BinEdges, 1e-3, 0.1, 1);
            testCase.verifyEqual(numel(result.d_histogram), numel(testCase.BinEdges) - 1);
        end

        function test_dstar_mean_correct(testCase)
            % D* mean should ignore NaN from failed fits.
            d_vec = [0.001; 0.002; 0.003];
            f_vec = [0.1; 0; 0.15];       % voxel 2 failed
            dstar_vec = [0.01; 0; 0.02];  % voxel 2 failed
            adc_mask = true(3,1);

            result = compute_ivim_metrics(testCase.Config, d_vec, f_vec, dstar_vec, ...
                [], adc_mask, 1, 1, testCase.BinEdges, 1e-3, 0.1, 1);
            % After filtering, D* = [0.01, NaN, 0.02], mean = 0.015
            testCase.verifyEqual(result.dstar_mean_val, 0.015, 'AbsTol', 1e-10);
        end
    end
end
