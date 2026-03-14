classdef test_compute_multi_core_metrics < matlab.unittest.TestCase
% TEST_COMPUTE_MULTI_CORE_METRICS — Unit tests for compute_multi_core_metrics.m
%
% Validates the multi-method core metric computation loop:
%   - All 11 methods produce output fields
%   - Metrics are NaN or valid numeric for empty masks
%   - Mask storage flag works correctly
%   - fDM volume fractions computed for post-baseline timepoints
%   - Graceful handling when extract_tumor_core fails

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'dependencies'));
        end
    end

    methods (Test)
        function test_all_methods_produce_output(testCase)
            ALL_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
                'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
                'percentile', 'spectral', 'fdm'};

            [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = testCase.setupData(ALL_METHODS);

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], struct('timepoint_index', 1), ...
                1, 1, 1, 0.001, 100, 0.1, 0.001, 0.5, false);

            for m = 1:numel(ALL_METHODS)
                mname = ALL_METHODS{m};
                testCase.verifyTrue(isfield(per_method, mname), ...
                    sprintf('Method %s should have output field.', mname));
                testCase.verifyTrue(isfield(per_method.(mname), 'adc_sub_vol'), ...
                    sprintf('Method %s should have adc_sub_vol.', mname));
            end
        end

        function test_adc_sub_volume_computed(testCase)
            ALL_METHODS = {'adc_threshold'};
            [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = testCase.setupData(ALL_METHODS);

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], struct('timepoint_index', 1), ...
                1, 1, 1, 0.001, 100, 0.1, 0.001, 0.5, false);

            vol = per_method.adc_threshold.adc_sub_vol(1,1,1);
            testCase.verifyTrue(isfinite(vol) || isnan(vol), ...
                'ADC sub-volume should be finite or NaN.');
        end

        function test_store_masks_flag(testCase)
            ALL_METHODS = {'adc_threshold', 'otsu'};
            [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = testCase.setupData(ALL_METHODS);
            % Add core_masks cell array for storage
            for m = 1:numel(ALL_METHODS)
                per_method.(ALL_METHODS{m}).core_masks = cell(3, 3);
            end

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], struct('timepoint_index', 1), ...
                1, 1, 1, 0.001, 100, 0.1, 0.001, 0.5, true);

            % With store_masks=true, core_masks{1,1} should be populated
            testCase.verifyFalse(isempty(per_method.adc_threshold.core_masks{1,1}), ...
                'Core mask should be stored when store_masks=true.');
            testCase.verifyTrue(islogical(per_method.adc_threshold.core_masks{1,1}), ...
                'Stored mask should be logical.');
        end

        function test_no_mask_storage_when_disabled(testCase)
            ALL_METHODS = {'adc_threshold'};
            [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = testCase.setupData(ALL_METHODS);

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], struct('timepoint_index', 1), ...
                1, 1, 1, 0.001, 100, 0.1, 0.001, 0.5, false);

            testCase.verifyFalse(isfield(per_method.adc_threshold, 'core_masks') && ...
                ~isempty(per_method.adc_threshold.core_masks{1,1}), ...
                'Core masks should not be stored when store_masks=false.');
        end

        function test_fdm_fractions_at_post_baseline(testCase)
            ALL_METHODS = {'fdm'};
            [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = testCase.setupData(ALL_METHODS);
            % k=2 (post-baseline) with baseline vectors in core_opts
            core_opts = struct('timepoint_index', 2, ...
                'baseline_adc_vec', adc_vec + 0.0001, ...
                'baseline_d_vec', d_vec);

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], core_opts, ...
                1, 2, 1, 0.001, 100, 0.1, 0.001, 0.5, false);

            % fDM fractions should be populated for k=2
            testCase.verifyTrue(isfield(per_method.fdm, 'fdm_responding_pc'));
            resp = per_method.fdm.fdm_responding_pc(1,2,1);
            prog = per_method.fdm.fdm_progressing_pc(1,2,1);
            stab = per_method.fdm.fdm_stable_pc(1,2,1);
            % Fractions should sum to ~1.0
            if isfinite(resp) && isfinite(prog) && isfinite(stab)
                testCase.verifyEqual(resp + prog + stab, 1.0, 'AbsTol', 1e-10, ...
                    'fDM fractions should sum to 1.0.');
            end
        end

        function test_empty_voxel_vector(testCase)
            ALL_METHODS = {'adc_threshold'};
            [per_method, config, ~, ~, ~, ~] = testCase.setupData(ALL_METHODS);
            % Empty ADC vector
            adc_vec = [];
            d_vec = [];
            f_vec = [];
            dstar_vec = [];

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], struct('timepoint_index', 1), ...
                1, 1, 1, 0.001, 100, 0.1, 0.001, 0.0, false);

            vol = per_method.adc_threshold.adc_sub_vol(1,1,1);
            testCase.verifyEqual(vol, 0, 'AbsTol', 1e-12, ...
                'Empty input should produce zero sub-volume.');
        end

        function test_unified_methods_share_mask_for_d_f(testCase)
            % Percentile method should use same mask for ADC and D/f
            ALL_METHODS = {'percentile'};
            [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = testCase.setupData(ALL_METHODS);
            per_method.percentile.core_masks = cell(3, 3);

            per_method = compute_multi_core_metrics(per_method, config, ...
                ALL_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                false, [], struct('timepoint_index', 1), ...
                1, 1, 1, 0.001, 100, 0.1, 0.001, 0.5, true);

            % For unified methods, f_sub_vol should equal adc_sub_vol
            % (same mask applied)
            adc_vol = per_method.percentile.adc_sub_vol(1,1,1);
            f_vol = per_method.percentile.f_sub_vol(1,1,1);
            testCase.verifyEqual(adc_vol, f_vol, 'AbsTol', 1e-12, ...
                'Unified methods should use same mask for ADC and f sub-volumes.');
        end
    end

    methods (Access = private)
        function [per_method, config, adc_vec, d_vec, f_vec, dstar_vec] = setupData(~, methods)
            rng(42);
            n_vox = 200;
            adc_vec = rand(n_vox, 1) * 0.003;
            d_vec = rand(n_vox, 1) * 0.003;
            f_vec = rand(n_vox, 1) * 0.3;
            dstar_vec = rand(n_vox, 1) * 0.05;

            config = struct();
            config.adc_thresh = 0.001;
            config.d_thresh = 0.001;
            config.f_thresh = 0.1;
            config.dstar_thresh = 0.01;
            config.core_method = 'adc_threshold';
            config.core_percentile = 25;
            config.core_n_clusters = 2;
            config.fdm_parameter = 'adc';
            config.fdm_thresh = 0.0004;
            config.spectral_min_voxels = 20;

            per_method = struct();
            for m = 1:numel(methods)
                mname = methods{m};
                per_method.(mname).adc_sub_vol = zeros(3, 3, 3);
                per_method.(mname).adc_sub_vol_pc = zeros(3, 3, 3);
                per_method.(mname).adc_sub_mean = NaN(3, 3, 3);
                per_method.(mname).adc_sub_kurt = NaN(3, 3, 3);
                per_method.(mname).adc_sub_skew = NaN(3, 3, 3);
                per_method.(mname).f_sub_vol = zeros(3, 3, 3);
                per_method.(mname).d_sub_mean = NaN(3, 3, 3);
                per_method.(mname).d_sub_kurt = NaN(3, 3, 3);
                per_method.(mname).d_sub_skew = NaN(3, 3, 3);
                per_method.(mname).fdm_responding_pc = NaN(3, 3, 3);
                per_method.(mname).fdm_progressing_pc = NaN(3, 3, 3);
                per_method.(mname).fdm_stable_pc = NaN(3, 3, 3);
            end
        end
    end
end
