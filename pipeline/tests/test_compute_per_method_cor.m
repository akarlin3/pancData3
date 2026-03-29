classdef test_compute_per_method_cor < matlab.unittest.TestCase
    % TEST_COMPUTE_PER_METHOD_COR  Tests for per-method Coefficient of Reproducibility.

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
                buildCorTestData(testCase.TempDir);
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
            result = compute_per_method_cor(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            expected = {'method_names', 'median_wcv', 'cor', ...
                'n_patients_with_repeats', 'per_patient_wcv', 'subvol_per_repeat'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(result, expected{i}), ...
                    sprintf('Output should contain field %s.', expected{i}));
            end
        end

        function testCorPositive(testCase)
            % CoR should be > 0 for methods with varied repeat data.
            result = compute_per_method_cor(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            % At least adc_threshold should have computable CoR
            valid_cor = result.cor(~isnan(result.cor));
            testCase.verifyGreaterThan(numel(valid_cor), 0, ...
                'At least one method should have computable CoR.');
            testCase.verifyTrue(all(valid_cor >= 0), ...
                'All valid CoR values should be >= 0.');
        end

        function testIdenticalRepeatsNearZeroCor(testCase)
            % Patients with identical repeat data should have CoR ≈ 0.
            [cfg, dv, ids, gtv] = buildIdenticalRepeatData(testCase.TempDir);
            result = compute_per_method_cor(dv, cfg, ids, gtv);

            % adc_threshold (deterministic) with identical data should give CoR near 0
            adc_idx = find(strcmp(result.method_names, 'adc_threshold'));
            if ~isempty(adc_idx) && ~isnan(result.cor(adc_idx))
                testCase.verifyLessThan(result.cor(adc_idx), 5, ...
                    'Identical repeat data should give CoR near 0 for deterministic methods.');
            end
        end

        function testOutputDimensions(testCase)
            % Verify output dimensions match expected sizes.
            result = compute_per_method_cor(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            n_methods = numel(result.method_names);
            n_patients = numel(testCase.IdList);
            n_repeats = size(testCase.DataVectors, 3);

            testCase.verifyEqual(size(result.median_wcv), [n_methods, 1]);
            testCase.verifyEqual(size(result.cor), [n_methods, 1]);
            testCase.verifyEqual(size(result.per_patient_wcv, 1), n_methods);
            testCase.verifyEqual(size(result.per_patient_wcv, 2), n_patients);
            testCase.verifyEqual(size(result.subvol_per_repeat, 3), n_repeats);
        end

        function testPatientsWithRepeatsCount(testCase)
            % Should detect patients with repeat scans.
            result = compute_per_method_cor(testCase.DataVectors, ...
                testCase.ConfigStruct, testCase.IdList, testCase.GtvLocations);

            testCase.verifyGreaterThan(result.n_patients_with_repeats, 0, ...
                'Should find at least one patient with repeat scans.');
        end

    end
end


function [cfg, dv, id_list, gtv_locations] = buildCorTestData(tempDir)
% Build mock test data with repeat scans for CoR tests.

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
    n_patients = 5;
    n_repeats = 3;

    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
        'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
        'dstar_vector_ivimnet', [], ...
        'vox_vol', [], 'vox_dims', []);

    % [nPatients x nTp x nRepeats] - Fx1 only (1 timepoint), 3 repeats
    dv = repmat(empty_entry, n_patients, 1, n_repeats);

    for j = 1:n_patients
        base_adc = [0.0005 + 0.0001*randn(25,1); 0.0015 + 0.0002*randn(55,1)];
        base_d = [0.0004 + 0.0001*randn(25,1); 0.0014 + 0.0002*randn(55,1)];
        base_f = [0.05 + 0.01*randn(25,1); 0.2 + 0.05*randn(55,1)];
        base_dstar = 0.02 * ones(n_vox, 1);

        for rpi = 1:n_repeats
            noise_scale = 0.0002 * (rpi - 1);
            dv(j, 1, rpi).adc_vector = base_adc + noise_scale * randn(n_vox, 1);
            dv(j, 1, rpi).d_vector = base_d + noise_scale * randn(n_vox, 1);
            dv(j, 1, rpi).f_vector = base_f + noise_scale * 10 * randn(n_vox, 1);
            dv(j, 1, rpi).dstar_vector = base_dstar;
        end
    end

    id_list = arrayfun(@(x) sprintf('P%02d', x), 1:n_patients, 'UniformOutput', false);
    gtv_locations = cell(n_patients, 1, n_repeats);
end


function [cfg, dv, id_list, gtv_locations] = buildIdenticalRepeatData(tempDir)
% Build data where all repeats are identical (should give CoR ≈ 0).

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

    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
        'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
        'dstar_vector_ivimnet', [], ...
        'vox_vol', [], 'vox_dims', []);

    dv = repmat(empty_entry, 3, 1, 2);

    for j = 1:3
        adc = [0.0005 * ones(30,1); 0.0015 * ones(50,1)];
        d = [0.0004 * ones(30,1); 0.0014 * ones(50,1)];
        f = [0.05 * ones(30,1); 0.2 * ones(50,1)];
        dstar = 0.02 * ones(n_vox, 1);

        for rpi = 1:2
            dv(j, 1, rpi).adc_vector = adc;
            dv(j, 1, rpi).d_vector = d;
            dv(j, 1, rpi).f_vector = f;
            dv(j, 1, rpi).dstar_vector = dstar;
        end
    end

    id_list = {'P01', 'P02', 'P03'};
    gtv_locations = cell(3, 1, 2);
end
