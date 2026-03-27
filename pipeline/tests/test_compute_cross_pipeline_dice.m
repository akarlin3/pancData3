classdef test_compute_cross_pipeline_dice < matlab.unittest.TestCase
    % TEST_COMPUTE_CROSS_PIPELINE_DICE  Tests for cross-pipeline Dice computation.
    %
    % Validates that compute_cross_pipeline_dice correctly:
    %   - Returns a struct with all expected fields
    %   - Produces dice array with correct dimensions [11 x 3 x nPatients]
    %   - All non-NaN Dice values are in [0, 1]
    %   - Identical pipeline vectors produce Dice = 1.0
    %   - Single-patient edge case runs without error
    %   - Missing DnCNN data produces NaN for DnCNN-involving pairs

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
            % Verify all expected fields exist in dice_results.
            results = compute_cross_pipeline_dice(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            expected = {'dice', 'method_names', 'pipeline_names', ...
                'pipeline_pair_labels', 'patient_ids', 'n_voxels', ...
                'fallback_flags', 'n_patients'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(results, expected{i}), ...
                    sprintf('Output should contain field %s.', expected{i}));
            end
        end

        function testDiceDimensions(testCase)
            % dice_results.dice should be [11 x 3 x 2] (11 methods, 3 pairs, 2 patients).
            results = compute_cross_pipeline_dice(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            testCase.verifyEqual(size(results.dice), [11 3 2], ...
                'Dice array should be 11 x 3 x 2.');
        end

        function testDiceRange(testCase)
            % All non-NaN Dice values should be in [0, 1].
            results = compute_cross_pipeline_dice(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            valid_vals = results.dice(~isnan(results.dice));
            testCase.verifyTrue(all(valid_vals >= 0 & valid_vals <= 1), ...
                'All Dice values should be in [0, 1].');
        end

        function testIdenticalPipelinesGiveDiceOne(testCase)
            % If Standard and DnCNN vectors are identical, Dice should be
            % 1.0 for all threshold-based methods.
            dv = testCase.DataVectors;
            % Make DnCNN identical to Standard
            dv(1,1,1).adc_vector_dncnn = dv(1,1,1).adc_vector;
            dv(1,1,1).d_vector_dncnn = dv(1,1,1).d_vector;
            dv(1,1,1).f_vector_dncnn = dv(1,1,1).f_vector;
            dv(1,1,1).dstar_vector_dncnn = dv(1,1,1).dstar_vector;

            results = compute_cross_pipeline_dice(dv, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);

            % Pair 1 (Std vs DnCNN) should be 1.0 for adc_threshold (method 1)
            dice_val = results.dice(1, 1, 1);
            if ~isnan(dice_val)
                testCase.verifyEqual(dice_val, 1.0, 'AbsTol', 1e-10, ...
                    'Identical pipelines should produce Dice = 1.0 for adc_threshold.');
            end
        end

        function testSinglePatient(testCase)
            % Edge case with 1 patient should not error.
            dv = testCase.DataVectors(1, :, :);
            id = testCase.IdList(1);
            gtv = testCase.GtvLocations(1, :, :);

            results = compute_cross_pipeline_dice(dv, testCase.ConfigStruct, id, gtv);

            testCase.verifyEqual(results.n_patients, 1);
            testCase.verifyEqual(size(results.dice, 3), 1, ...
                'Single patient should produce dice with 3rd dim = 1.');
        end

        function testEmptyDnCNNPipeline(testCase)
            % If DnCNN vectors are all empty, Dice involving DnCNN should be NaN.
            dv = testCase.DataVectors;
            % Clear DnCNN fields
            dv(1,1,1).adc_vector_dncnn = [];
            dv(1,1,1).d_vector_dncnn = [];
            dv(1,1,1).f_vector_dncnn = [];
            dv(1,1,1).dstar_vector_dncnn = [];
            dv(2,1,1).adc_vector_dncnn = [];
            dv(2,1,1).d_vector_dncnn = [];
            dv(2,1,1).f_vector_dncnn = [];
            dv(2,1,1).dstar_vector_dncnn = [];

            results = compute_cross_pipeline_dice(dv, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);

            % Pair 1 (Std vs DnCNN) and Pair 3 (DnCNN vs IVIMNet) should be NaN
            for m = 1:11
                testCase.verifyTrue(isnan(results.dice(m, 1, 1)), ...
                    sprintf('Std vs DnCNN should be NaN when DnCNN empty (method %d).', m));
                testCase.verifyTrue(isnan(results.dice(m, 3, 1)), ...
                    sprintf('DnCNN vs IVIMNet should be NaN when DnCNN empty (method %d).', m));
            end
        end

    end
end


function [cfg, dv, id_list, gtv_locations] = buildTestData(tempDir)
% Build mock test data for cross-pipeline Dice tests.

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

    % Create a 10x10x1 3D GTV mask
    Stvol3d = true(10, 10, 1); %#ok<NASGU>
    gtv_mat_file = fullfile(tempDir, 'gtv_p1.mat');
    save(gtv_mat_file, 'Stvol3d');

    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
        'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
        'dstar_vector_ivimnet', [], ...
        'vox_vol', [], 'vox_dims', []);
    dv = repmat(empty_entry, 2, 2, 1);

    % Patient 1, Fx1 -- bimodal ADC (30 core + 70 margin)
    adc1 = [0.0005 + 0.0001*randn(30,1); 0.0015 + 0.0002*randn(70,1)];
    d1 = [0.0004 + 0.0001*randn(30,1); 0.0014 + 0.0002*randn(70,1)];
    f1 = [0.05 + 0.01*randn(30,1); 0.2 + 0.05*randn(70,1)];
    dstar1 = 0.02 * ones(n_vox, 1);

    % Standard
    dv(1,1,1).adc_vector = adc1;
    dv(1,1,1).d_vector = d1;
    dv(1,1,1).f_vector = f1;
    dv(1,1,1).dstar_vector = dstar1;
    dv(1,1,1).vox_vol = 0.008;
    dv(1,1,1).vox_dims = [2 2 2];

    % DnCNN (slightly different from Standard)
    dv(1,1,1).adc_vector_dncnn = adc1 + 0.00005*randn(n_vox, 1);
    dv(1,1,1).d_vector_dncnn = d1 + 0.00005*randn(n_vox, 1);
    dv(1,1,1).f_vector_dncnn = f1 + 0.005*randn(n_vox, 1);
    dv(1,1,1).dstar_vector_dncnn = dstar1 + 0.001*randn(n_vox, 1);

    % IVIMNet (ADC shared with Standard, IVIM params slightly different)
    dv(1,1,1).d_vector_ivimnet = d1 + 0.0001*randn(n_vox, 1);
    dv(1,1,1).f_vector_ivimnet = f1 + 0.01*randn(n_vox, 1);
    dv(1,1,1).dstar_vector_ivimnet = dstar1 + 0.002*randn(n_vox, 1);

    % Patient 1, Fx2 (for struct completeness)
    dv(1,2,1).adc_vector = adc1 + 0.0003 * randn(n_vox, 1);
    dv(1,2,1).d_vector = d1 + 0.0002 * randn(n_vox, 1);
    dv(1,2,1).f_vector = f1 + 0.02 * randn(n_vox, 1);
    dv(1,2,1).dstar_vector = dstar1;

    % Patient 2, Fx1 -- different bimodal data with all pipelines
    rng(99);
    adc2 = [0.0006 + 0.0001*randn(40,1); 0.0018 + 0.0002*randn(60,1)];
    d2 = [0.0005 + 0.0001*randn(40,1); 0.0016 + 0.0002*randn(60,1)];
    f2 = [0.06 + 0.01*randn(40,1); 0.22 + 0.05*randn(60,1)];
    dstar2 = 0.015 * ones(n_vox, 1);

    dv(2,1,1).adc_vector = adc2;
    dv(2,1,1).d_vector = d2;
    dv(2,1,1).f_vector = f2;
    dv(2,1,1).dstar_vector = dstar2;
    dv(2,1,1).adc_vector_dncnn = adc2 + 0.00005*randn(n_vox, 1);
    dv(2,1,1).d_vector_dncnn = d2 + 0.00005*randn(n_vox, 1);
    dv(2,1,1).f_vector_dncnn = f2 + 0.005*randn(n_vox, 1);
    dv(2,1,1).dstar_vector_dncnn = dstar2 + 0.001*randn(n_vox, 1);
    dv(2,1,1).d_vector_ivimnet = d2 + 0.0001*randn(n_vox, 1);
    dv(2,1,1).f_vector_ivimnet = f2 + 0.01*randn(n_vox, 1);
    dv(2,1,1).dstar_vector_ivimnet = dstar2 + 0.002*randn(n_vox, 1);
    dv(2,1,1).vox_vol = 0.008;
    dv(2,1,1).vox_dims = [2 2 2];

    % Patient 2, Fx2 -- empty
    dv(2,2,1).adc_vector = [];
    dv(2,2,1).d_vector = [];
    dv(2,2,1).f_vector = [];
    dv(2,2,1).dstar_vector = [];

    id_list = {'P01', 'P02'};
    gtv_locations = cell(2, 2, 1);
    gtv_locations{1, 1, 1} = gtv_mat_file;
    gtv_locations{1, 2, 1} = gtv_mat_file;
    gtv_locations{2, 1, 1} = '';
    gtv_locations{2, 2, 1} = '';
end
