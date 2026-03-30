classdef test_compute_subvolume_stability < matlab.unittest.TestCase
    % TEST_COMPUTE_SUBVOLUME_STABILITY  Tests for sub-volume stability over fractions.

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
                buildStabilityTestData(testCase.TempDir);
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
            result = compute_subvolume_stability(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            expected = {'dice_vs_baseline', 'method_names', 'n_patients', 'n_timepoints'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(result, expected{i}), ...
                    sprintf('Output should contain field %s.', expected{i}));
            end
        end

        function testFx1DiceIsOne(testCase)
            % Fx1 Dice vs Fx1 should be 1.0 for all methods with non-empty masks.
            result = compute_subvolume_stability(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            dice_fx1 = result.dice_vs_baseline(:, 1, :);
            valid = dice_fx1(~isnan(dice_fx1));
            testCase.verifyTrue(all(abs(valid - 1.0) < 1e-10), ...
                'All Fx1 Dice values should be 1.0.');
        end

        function testOutputDimensions(testCase)
            result = compute_subvolume_stability(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            n_methods = numel(result.method_names);
            n_patients = numel(testCase.IdList);
            nTp = size(testCase.DataVectors, 2);

            testCase.verifyEqual(size(result.dice_vs_baseline), [n_methods, nTp, n_patients]);
        end

        function testDiceRange(testCase)
            % All Dice values should be in [0, 1] or NaN.
            result = compute_subvolume_stability(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            dice_all = result.dice_vs_baseline(:);
            valid = dice_all(~isnan(dice_all));
            testCase.verifyTrue(all(valid >= 0 & valid <= 1), ...
                'All valid Dice values should be in [0, 1].');
        end

    end
end


function [cfg, dv, id_list, gtv_locations] = buildStabilityTestData(tempDir)
% Build mock test data with 3 patients x 3 timepoints for stability tests.

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';
    cfg.adc_thresh = 0.001;
    cfg.high_adc_thresh = 0.00115;
    cfg.d_thresh = 0.001;
    cfg.f_thresh = 0.1;
    cfg.dstar_thresh = 0.01;
    cfg.min_vox_hist = 10;
    cfg.core_method = 'adc_threshold';
    cfg.core_percentile = 25;
    cfg.core_n_clusters = 2;
    cfg.fdm_parameter = 'adc';
    cfg.fdm_thresh = 0.0004;
    cfg.adc_max = 0.003;
    cfg.spectral_min_voxels = 20;

    rng(42);
    n_vox = 80;
    n_patients = 3;
    nTp = 3;

    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
        'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
        'dstar_vector_ivimnet', [], ...
        'vox_vol', [], 'vox_dims', []);
    dv = repmat(empty_entry, n_patients, nTp, 1);

    for j = 1:n_patients
        base_adc = [0.0005 + 0.0001*randn(25,1); 0.0015 + 0.0002*randn(55,1)];
        base_d = [0.0004 + 0.0001*randn(25,1); 0.0014 + 0.0002*randn(55,1)];
        base_f = [0.05 + 0.01*randn(25,1); 0.2 + 0.05*randn(55,1)];
        base_dstar = 0.02 * ones(n_vox, 1);

        for k = 1:nTp
            % Add increasing noise to simulate treatment response
            noise = 0.0003 * (k - 1);
            dv(j, k, 1).adc_vector = base_adc + noise * randn(n_vox, 1);
            dv(j, k, 1).d_vector = base_d + noise * randn(n_vox, 1);
            dv(j, k, 1).f_vector = base_f + noise * 10 * randn(n_vox, 1);
            dv(j, k, 1).dstar_vector = base_dstar;
        end
    end

    id_list = {'P01', 'P02', 'P03'};
    gtv_locations = cell(n_patients, nTp, 1);
end
