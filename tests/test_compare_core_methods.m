classdef test_compare_core_methods < matlab.unittest.TestCase
    % TEST_COMPARE_CORE_METHODS  Tests for the compare_core_methods module.

    properties
        TempDir
        ConfigStruct
        DataVectors
        SummaryMetrics
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

            % Suppress figures
            set(0, 'DefaultFigureVisible', 'off');

            % Build mock config
            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.output_folder = testCase.TempDir;
            testCase.ConfigStruct.dwi_types_to_run = 1;
            testCase.ConfigStruct.dwi_type_name = 'Standard';
            testCase.ConfigStruct.adc_thresh = 0.001;
            testCase.ConfigStruct.high_adc_thresh = 0.00115;
            testCase.ConfigStruct.d_thresh = 0.001;
            testCase.ConfigStruct.f_thresh = 0.1;
            testCase.ConfigStruct.dstar_thresh = 0.01;
            testCase.ConfigStruct.min_vox_hist = 20;  % low for small test data
            testCase.ConfigStruct.core_method = 'adc_threshold';
            testCase.ConfigStruct.core_percentile = 25;
            testCase.ConfigStruct.core_n_clusters = 2;
            testCase.ConfigStruct.fdm_parameter = 'adc';
            testCase.ConfigStruct.fdm_thresh = 0.0004;
            testCase.ConfigStruct.adc_max = 0.003;

            % Build mock data_vectors_gtvp (2 patients, 2 timepoints)
            rng(42);
            n_vox = 100;

            % Create 3D GTV mask and save to temp .mat
            Stvol3d = true(10, 10, 1); %#ok<NASGU>
            gtv_mat_file = fullfile(testCase.TempDir, 'gtv_p1.mat');
            save(gtv_mat_file, 'Stvol3d');

            % Init struct with all required fields
            empty_entry = struct( ...
                'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
                'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
                'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
                'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
                'dstar_vector_ivimnet', [], ...
                'vox_vol', [], 'vox_dims', []);
            testCase.DataVectors = repmat(empty_entry, 2, 2, 1);

            % Patient 1, Fx1 -- bimodal ADC (30 core + 70 margin)
            adc1 = [0.0005 + 0.0001*randn(30,1); 0.0015 + 0.0002*randn(70,1)];
            d1 = [0.0004 + 0.0001*randn(30,1); 0.0014 + 0.0002*randn(70,1)];
            f1 = [0.05 + 0.01*randn(30,1); 0.2 + 0.05*randn(70,1)];
            dstar1 = 0.02 * ones(n_vox, 1);
            testCase.DataVectors(1,1,1).adc_vector = adc1;
            testCase.DataVectors(1,1,1).d_vector = d1;
            testCase.DataVectors(1,1,1).f_vector = f1;
            testCase.DataVectors(1,1,1).dstar_vector = dstar1;
            testCase.DataVectors(1,1,1).vox_vol = 0.008;
            testCase.DataVectors(1,1,1).vox_dims = [2 2 2];

            % Patient 1, Fx2 (shifted ADC for fDM testing)
            testCase.DataVectors(1,2,1).adc_vector = adc1 + 0.0003 * randn(n_vox, 1);
            testCase.DataVectors(1,2,1).d_vector = d1 + 0.0002 * randn(n_vox, 1);
            testCase.DataVectors(1,2,1).f_vector = f1 + 0.02 * randn(n_vox, 1);
            testCase.DataVectors(1,2,1).dstar_vector = dstar1;
            testCase.DataVectors(1,2,1).vox_vol = 0.008;
            testCase.DataVectors(1,2,1).vox_dims = [2 2 2];

            % Patient 2, Fx1 -- all NaN (edge case)
            testCase.DataVectors(2,1,1).adc_vector = nan(50, 1);
            testCase.DataVectors(2,1,1).d_vector = nan(50, 1);
            testCase.DataVectors(2,1,1).f_vector = nan(50, 1);
            testCase.DataVectors(2,1,1).dstar_vector = nan(50, 1);
            testCase.DataVectors(2,1,1).vox_vol = 0.008;
            testCase.DataVectors(2,1,1).vox_dims = [2 2 2];

            % Patient 2, Fx2 -- empty
            testCase.DataVectors(2,2,1).adc_vector = [];
            testCase.DataVectors(2,2,1).d_vector = [];
            testCase.DataVectors(2,2,1).f_vector = [];
            testCase.DataVectors(2,2,1).dstar_vector = [];

            % Summary metrics
            testCase.SummaryMetrics = struct();
            testCase.SummaryMetrics.id_list = {'P01', 'P02'};
            gtv_locs = cell(2, 2, 1);
            gtv_locs{1, 1, 1} = gtv_mat_file;
            gtv_locs{1, 2, 1} = gtv_mat_file;
            gtv_locs{2, 1, 1} = '';
            gtv_locs{2, 2, 1} = '';
            testCase.SummaryMetrics.gtv_locations = gtv_locs;
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            diary off;  % close diary before rmdir (Windows file locking)
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)

        function testBasicExecution(testCase)
            % Should run without error and return a struct with 11 methods
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            testCase.verifyTrue(isstruct(results), 'Output should be a struct.');
            testCase.verifyEqual(numel(results.method_names), 11, ...
                'Should have 11 methods.');
        end

        function testDiceMatrixProperties(testCase)
            % Dice matrix should be symmetric, in [0,1], with diagonal = 1
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            D = results.mean_dice_matrix;
            % Diagonal should be 1
            for i = 1:11
                if ~isnan(D(i,i))
                    testCase.verifyEqual(D(i,i), 1, 'AbsTol', 1e-10, ...
                        sprintf('Self-Dice for method %d should be 1.', i));
                end
            end
            % Symmetry
            for i = 1:11
                for j2 = (i+1):11
                    testCase.verifyEqual(D(i,j2), D(j2,i), 'AbsTol', 1e-12, ...
                        'Dice matrix should be symmetric.');
                end
            end
            % Range
            valid_vals = D(~isnan(D));
            testCase.verifyTrue(all(valid_vals >= 0 & valid_vals <= 1), ...
                'All Dice values should be in [0, 1].');
        end

        function testVolumeFractionsRange(testCase)
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            vf = results.volume_fractions;
            valid_vf = vf(~isnan(vf));
            testCase.verifyTrue(all(valid_vf >= 0 & valid_vf <= 1), ...
                'Volume fractions should be in [0, 1].');
        end

        function testFigureGeneration(testCase)
            compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            dice_fig = fullfile(testCase.TempDir, 'core_method_dice_heatmap_Standard.png');
            testCase.verifyTrue(exist(dice_fig, 'file') > 0, ...
                'Dice heatmap PNG should be created.');

            vol_fig = fullfile(testCase.TempDir, 'core_method_volume_comparison_Standard.png');
            testCase.verifyTrue(exist(vol_fig, 'file') > 0, ...
                'Volume comparison PNG should be created.');
        end

        function testMATFileSaved(testCase)
            compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            mat_file = fullfile(testCase.TempDir, 'compare_core_results_Standard.mat');
            testCase.verifyTrue(exist(mat_file, 'file') > 0, ...
                'Results MAT file should be created.');

            loaded = load(mat_file, 'compare_results');
            testCase.verifyTrue(isfield(loaded, 'compare_results'), ...
                'MAT file should contain compare_results struct.');
        end

        function testAllNaNPatientSkipped(testCase)
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            % Patient 2 has all-NaN data; volume fractions should be NaN
            for m = 1:11
                testCase.verifyTrue(isnan(results.volume_fractions(2, 1, m)), ...
                    sprintf('All-NaN patient volume fraction should be NaN (method %d).', m));
            end
        end

        function testFallbackDetection(testCase)
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            % fDM at baseline (k=1) should be flagged as fallback
            fdm_idx = find(strcmp(results.method_names, 'fdm'));
            testCase.verifyTrue(results.fallback_flags(1, 1, fdm_idx), ...
                'fDM at baseline should be flagged as fallback.');
        end

    end
end
