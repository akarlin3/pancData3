classdef test_compute_core_failure_rates < matlab.unittest.TestCase
    % TEST_COMPUTE_CORE_FAILURE_RATES  Tests for the core failure rate aggregation.
    %
    % Validates that compute_core_failure_rates correctly:
    %   - Returns a struct with all expected fields
    %   - Produces rate matrices with correct dimensions [11 x 3]
    %   - All rates are in [0, 1]
    %   - Known all-NaN patients contribute to all_nan_rate
    %   - adc_threshold never reports fallback

    properties
        TempDir
        ConfigStruct
        DataVectors
        IdList
        GtvLocations
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'core'));
            addpath(fullfile(repoRoot, 'utils'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end

            set(0, 'DefaultFigureVisible', 'off');

            [testCase.ConfigStruct, testCase.DataVectors, testCase.IdList, testCase.GtvLocations] = ...
                buildTestData(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            diary off;
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)

        function testOutputStructFields(testCase)
            % Verify all expected fields exist.
            ft = compute_core_failure_rates(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            expected = {'method_names', 'pipeline_names', 'n_patients', ...
                'n_timepoints', 'fallback_rate', 'empty_rate', ...
                'insufficient_rate', 'all_nan_rate', 'any_failure_rate', ...
                'median_core_voxels', 'core_voxel_counts', ...
                'was_fallback', 'was_empty', 'was_insufficient', 'was_all_nan'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(ft, expected{i}), ...
                    sprintf('Output should contain field %s.', expected{i}));
            end
        end

        function testRateDimensions(testCase)
            % All rate matrices should be [11 x 3].
            ft = compute_core_failure_rates(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            testCase.verifyEqual(size(ft.fallback_rate), [11 3]);
            testCase.verifyEqual(size(ft.empty_rate), [11 3]);
            testCase.verifyEqual(size(ft.insufficient_rate), [11 3]);
            testCase.verifyEqual(size(ft.all_nan_rate), [11 3]);
            testCase.verifyEqual(size(ft.any_failure_rate), [11 3]);
        end

        function testRatesInZeroOneRange(testCase)
            % All non-NaN rates should be in [0, 1].
            ft = compute_core_failure_rates(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            rate_fields = {'fallback_rate', 'empty_rate', 'insufficient_rate', ...
                'all_nan_rate', 'any_failure_rate'};
            for i = 1:numel(rate_fields)
                vals = ft.(rate_fields{i});
                valid = vals(~isnan(vals));
                testCase.verifyTrue(all(valid >= 0 & valid <= 1), ...
                    sprintf('%s should be in [0, 1].', rate_fields{i}));
            end
        end

        function testKnownNaNPatient(testCase)
            % Patient 2 has all-NaN Standard vectors; all_nan_rate should be > 0.
            ft = compute_core_failure_rates(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            % Standard pipeline (p=1): patient 2 is all-NaN
            % all_nan_rate for adc_threshold (method 1) should be > 0
            testCase.verifyGreaterThan(ft.all_nan_rate(1, 1), 0, ...
                'all_nan_rate should be > 0 when a patient has all-NaN data.');
        end

        function testAdcThresholdNeverFallback(testCase)
            % adc_threshold (method 1) fallback_rate should always be 0.
            ft = compute_core_failure_rates(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            for p = 1:3
                val = ft.fallback_rate(1, p);
                if ~isnan(val)
                    testCase.verifyEqual(val, 0, 'AbsTol', 1e-10, ...
                        sprintf('adc_threshold should never fall back (pipeline %d).', p));
                end
            end
        end

    end
end


function [cfg, dv, id_list, gtv_locations] = buildTestData(tempDir)
% Build mock test data for core failure rate tests.

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';
    cfg.adc_thresh = 0.001;
    cfg.high_adc_thresh = 0.00115;
    cfg.d_thresh = 0.001;
    cfg.f_thresh = 0.1;
    cfg.dstar_thresh = 0.01;
    cfg.min_vox_hist = 20;
    cfg.core_method = 'adc_threshold';
    cfg.core_percentile = 25;
    cfg.core_n_clusters = 2;
    cfg.fdm_parameter = 'adc';
    cfg.fdm_thresh = 0.0004;
    cfg.adc_max = 0.003;

    rng(42);
    n_vox = 100;

    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
        'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
        'dstar_vector_ivimnet', [], ...
        'vox_vol', [], 'vox_dims', []);
    dv = repmat(empty_entry, 2, 2, 1);

    % Patient 1, Fx1 -- bimodal data with all pipelines
    adc1 = [0.0005 + 0.0001*randn(30,1); 0.0015 + 0.0002*randn(70,1)];
    d1 = [0.0004 + 0.0001*randn(30,1); 0.0014 + 0.0002*randn(70,1)];
    f1 = [0.05 + 0.01*randn(30,1); 0.2 + 0.05*randn(70,1)];
    dstar1 = 0.02 * ones(n_vox, 1);

    dv(1,1,1).adc_vector = adc1;
    dv(1,1,1).d_vector = d1;
    dv(1,1,1).f_vector = f1;
    dv(1,1,1).dstar_vector = dstar1;
    dv(1,1,1).adc_vector_dncnn = adc1 + 0.00005*randn(n_vox, 1);
    dv(1,1,1).d_vector_dncnn = d1 + 0.00005*randn(n_vox, 1);
    dv(1,1,1).f_vector_dncnn = f1 + 0.005*randn(n_vox, 1);
    dv(1,1,1).dstar_vector_dncnn = dstar1;
    dv(1,1,1).d_vector_ivimnet = d1 + 0.0001*randn(n_vox, 1);
    dv(1,1,1).f_vector_ivimnet = f1 + 0.01*randn(n_vox, 1);
    dv(1,1,1).dstar_vector_ivimnet = dstar1;

    % Patient 1, Fx2
    dv(1,2,1).adc_vector = adc1 + 0.0003*randn(n_vox, 1);
    dv(1,2,1).d_vector = d1 + 0.0002*randn(n_vox, 1);
    dv(1,2,1).f_vector = f1 + 0.02*randn(n_vox, 1);
    dv(1,2,1).dstar_vector = dstar1;
    dv(1,2,1).adc_vector_dncnn = adc1 + 0.0004*randn(n_vox, 1);
    dv(1,2,1).d_vector_dncnn = d1 + 0.0003*randn(n_vox, 1);
    dv(1,2,1).f_vector_dncnn = f1 + 0.03*randn(n_vox, 1);
    dv(1,2,1).dstar_vector_dncnn = dstar1;
    dv(1,2,1).d_vector_ivimnet = d1 + 0.0002*randn(n_vox, 1);
    dv(1,2,1).f_vector_ivimnet = f1 + 0.02*randn(n_vox, 1);
    dv(1,2,1).dstar_vector_ivimnet = dstar1;

    % Patient 2, Fx1 -- all NaN for Standard (triggers all_nan failure)
    dv(2,1,1).adc_vector = nan(50, 1);
    dv(2,1,1).d_vector = nan(50, 1);
    dv(2,1,1).f_vector = nan(50, 1);
    dv(2,1,1).dstar_vector = nan(50, 1);

    % Patient 2, Fx2 -- empty
    dv(2,2,1).adc_vector = [];
    dv(2,2,1).d_vector = [];
    dv(2,2,1).f_vector = [];
    dv(2,2,1).dstar_vector = [];

    id_list = {'P01', 'P02'};
    gtv_locations = cell(2, 2, 1);
    gtv_locations{1, 1, 1} = '';
    gtv_locations{1, 2, 1} = '';
    gtv_locations{2, 1, 1} = '';
    gtv_locations{2, 2, 1} = '';
end
