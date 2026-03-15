classdef test_compare_core_methods < matlab.unittest.TestCase
    % TEST_COMPARE_CORE_METHODS  Tests for the compare_core_methods module.
    %
    % Validates that compare_core_methods correctly:
    %   - Runs all 11 tumor core delineation methods and returns results
    %   - Produces a symmetric Dice matrix with values in [0,1] and
    %     diagonal entries equal to 1
    %   - Computes volume fractions in the valid range [0,1]
    %   - Generates expected output figures (Dice heatmap, volume comparison)
    %   - Saves results to a .mat file
    %   - Handles edge cases: all-NaN patient data, fDM fallback at baseline,
    %     mask size mismatches with Fx1 fallback logic

    properties
        TempDir
        ConfigStruct
        DataVectors
        SummaryMetrics
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Create isolated temp directory and add project paths
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'core'));
            addpath(fullfile(repoRoot, 'utils'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end

            % Suppress figures during testing
            set(0, 'DefaultFigureVisible', 'off');

            % Build mock config with all thresholds needed by extract_tumor_core
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

            % Build mock data_vectors_gtvp (2 patients, 2 timepoints).
            % Patient 1 has bimodal ADC (simulating core + margin).
            % Patient 2 has all-NaN data (edge case for graceful skipping).
            rng(42);
            n_vox = 100;

            % Create a 10x10x1 3D GTV mask and save to temp .mat file
            % (required by 3D-aware methods like active_contours)
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
            % Close diary first to release file locks (Windows),
            % then close figures and remove temp directory.
            diary off;
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)

        function testBasicExecution(testCase)
            % Verify that compare_core_methods runs without error and
            % returns a struct containing results for all 11 core methods.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            testCase.verifyTrue(isstruct(results), 'Output should be a struct.');
            testCase.verifyEqual(numel(results.method_names), 11, ...
                'Should have 11 methods.');
        end

        function testDiceMatrixProperties(testCase)
            % Verify mathematical properties of the mean Dice coefficient
            % matrix: symmetric, values in [0,1], self-comparison = 1.
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
            % Volume fractions (core volume / total GTV volume) must be
            % in [0, 1] for all valid (non-NaN) entries.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            vf = results.volume_fractions;
            valid_vf = vf(~isnan(vf));
            testCase.verifyTrue(all(valid_vf >= 0 & valid_vf <= 1), ...
                'Volume fractions should be in [0, 1].');
        end

        function testFigureGeneration(testCase)
            % Verify that compare_core_methods produces the expected
            % output PNG figures in the output folder.
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
            % Verify that results are persisted to a .mat file and the
            % saved struct contains the expected compare_results variable.
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
            % Patient 2 has all-NaN voxel data (simulating missing scan).
            % All 11 methods should produce NaN volume fractions for this
            % patient rather than erroring.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            % Patient 2 has all-NaN data; volume fractions should be NaN
            for m = 1:11
                testCase.verifyTrue(isnan(results.volume_fractions(2, 1, m)), ...
                    sprintf('All-NaN patient volume fraction should be NaN (method %d).', m));
            end
        end

        function testFallbackDetection(testCase)
            % The fDM (functional Diffusion Map) method requires two
            % timepoints to compute voxel-wise changes. At baseline (k=1),
            % there is no prior timepoint, so it should fall back to
            % adc_threshold and be flagged accordingly.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            % fDM at baseline (k=1) should be flagged as fallback
            fdm_idx = find(strcmp(results.method_names, 'fdm'));
            testCase.verifyTrue(results.fallback_flags(1, 1, fdm_idx), ...
                'fDM at baseline should be flagged as fallback.');
        end

        function testFx1MaskFallback(testCase)
            % When DIR warps Fx2 vectors into Fx1 space, the native Fx2
            % GTV mask has a different voxel count than the warped vectors.
            % The code should detect this mismatch and fall back to the
            % Fx1 mask (which matches the vector length), enabling 3D
            % methods like active_contours to run successfully.
            n_vox = 100;
            rng(42);

            % Fx1 mask: 100 voxels (matches the warped Fx2 vectors)
            Stvol3d = true(10, 10, 1); %#ok<NASGU>
            fx1_file = fullfile(testCase.TempDir, 'gtv_fx1_fallback.mat');
            save(fx1_file, 'Stvol3d');

            % Fx2 mask: 90 voxels (mismatch — simulates native fraction)
            Stvol3d = true(9, 10, 1); %#ok<NASGU>
            fx2_file = fullfile(testCase.TempDir, 'gtv_fx2_fallback.mat');
            save(fx2_file, 'Stvol3d');

            % Set gtv_locations: Fx1 → 100 voxels, Fx2 → 90 voxels
            sm = testCase.SummaryMetrics;
            sm.gtv_locations{1, 1, 1} = fx1_file;
            sm.gtv_locations{1, 2, 1} = fx2_file;

            % Vectors for both timepoints have 100 voxels (Fx1 space)
            dv = testCase.DataVectors;
            adc1 = [0.0005 + 0.0001*randn(30,1); 0.0015 + 0.0002*randn(70,1)];
            d1 = [0.0004 + 0.0001*randn(30,1); 0.0014 + 0.0002*randn(70,1)];
            f1 = [0.05 + 0.01*randn(30,1); 0.2 + 0.05*randn(70,1)];
            dstar1 = 0.02 * ones(n_vox, 1);
            dv(1,2,1).adc_vector = adc1;
            dv(1,2,1).d_vector = d1;
            dv(1,2,1).f_vector = f1;
            dv(1,2,1).dstar_vector = dstar1;

            results = compare_core_methods(dv, sm, testCase.ConfigStruct);

            % Active contours at Fx2 should NOT be flagged as fallback
            % because the Fx1 mask fallback provides a valid 3D mask.
            ac_idx = find(strcmp(results.method_names, 'active_contours'));
            testCase.verifyFalse(results.fallback_flags(1, 2, ac_idx), ...
                'Active contours at Fx2 should use Fx1 fallback mask, not fall back to ADC threshold.');
        end

        function testFx1MaskFallbackBothMismatch(testCase)
            % Edge case: when BOTH the native Fx2 mask AND the Fx1 mask
            % have voxel counts that differ from the vector length,
            % 3D methods (e.g., active_contours) cannot reconstruct
            % the volume and should gracefully fall back to a 1D method.
            n_vox = 80;  % Neither 100-voxel nor 90-voxel mask will match
            rng(42);

            % Fx1 mask: 100 voxels
            Stvol3d = true(10, 10, 1); %#ok<NASGU>
            fx1_file = fullfile(testCase.TempDir, 'gtv_fx1_nomatch.mat');
            save(fx1_file, 'Stvol3d');

            % Fx2 mask: 90 voxels
            Stvol3d = true(9, 10, 1); %#ok<NASGU>
            fx2_file = fullfile(testCase.TempDir, 'gtv_fx2_nomatch.mat');
            save(fx2_file, 'Stvol3d');

            sm = testCase.SummaryMetrics;
            sm.gtv_locations{1, 1, 1} = fx1_file;
            sm.gtv_locations{1, 2, 1} = fx2_file;

            % Vectors have 80 voxels — matches neither mask
            dv = testCase.DataVectors;
            dv(1,1,1).adc_vector = 0.001 * ones(n_vox, 1);
            dv(1,1,1).d_vector = 0.0008 * ones(n_vox, 1);
            dv(1,1,1).f_vector = 0.15 * ones(n_vox, 1);
            dv(1,1,1).dstar_vector = 0.01 * ones(n_vox, 1);
            dv(1,1,1).vox_vol = 0.008;
            dv(1,1,1).vox_dims = [2 2 2];
            dv(1,2,1).adc_vector = 0.001 * ones(n_vox, 1);
            dv(1,2,1).d_vector = 0.0008 * ones(n_vox, 1);
            dv(1,2,1).f_vector = 0.15 * ones(n_vox, 1);
            dv(1,2,1).dstar_vector = 0.01 * ones(n_vox, 1);
            dv(1,2,1).vox_vol = 0.008;
            dv(1,2,1).vox_dims = [2 2 2];

            % Should not error — 3D methods fall back gracefully
            results = compare_core_methods(dv, sm, testCase.ConfigStruct);

            % Active contours should be flagged as fallback at Fx2
            ac_idx = find(strcmp(results.method_names, 'active_contours'));
            testCase.verifyTrue(results.fallback_flags(1, 2, ac_idx), ...
                'Active contours should fall back when both masks mismatch.');
        end

        function testHD95MatrixProperties(testCase)
            % HD95 matrix should be symmetric with zero diagonal (self-
            % distance is always 0) and non-negative off-diagonal values.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            H = results.mean_hd95_matrix;
            % Check symmetry for non-NaN entries
            for i = 1:11
                for j2 = (i+1):11
                    if ~isnan(H(i,j2)) && ~isnan(H(j2,i))
                        testCase.verifyEqual(H(i,j2), H(j2,i), 'AbsTol', 1e-10, ...
                            'HD95 matrix should be symmetric.');
                    end
                end
            end
            % Diagonal should be 0 (self-comparison) for valid entries
            for i = 1:11
                if ~isnan(H(i,i))
                    testCase.verifyEqual(H(i,i), 0, 'AbsTol', 1e-10, ...
                        sprintf('HD95 self-distance for method %d should be 0.', i));
                end
            end
            % Non-negative
            valid_hd = H(~isnan(H));
            testCase.verifyTrue(all(valid_hd >= 0), ...
                'All HD95 values should be non-negative.');
        end

        function testDiceCountTracked(testCase)
            % dice_count should track how many valid observations were
            % used for each pairwise comparison.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            % Patient 1 Fx1 has valid data, so dice_count should be >= 1
            testCase.verifyTrue(all(results.dice_count(:) >= 0), ...
                'Dice count should be non-negative.');
            % At least one observation (Patient 1 Fx1)
            testCase.verifyTrue(any(results.dice_count(:) > 0), ...
                'At least one Dice observation should exist.');
        end

        function testResultsStructFields(testCase)
            % Verify all expected fields are present in output struct.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            expected = {'method_names', 'mean_dice_matrix', 'mean_hd95_matrix', ...
                'dice_count', 'hd95_count', 'volume_fractions', ...
                'fallback_flags', 'all_dice', 'all_hd95', ...
                'n_patients', 'nTp', 'dwi_type_name'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(results, expected{i}), ...
                    sprintf('Output should contain field %s.', expected{i}));
            end
        end

        function testDiceMatrixDimensions(testCase)
            % Mean Dice matrix should be 11x11 (one per method).
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            testCase.verifyEqual(size(results.mean_dice_matrix), [11 11], ...
                'Mean Dice matrix should be 11x11.');
        end

        function testVolumeFractionsDimensions(testCase)
            % Volume fractions should be nPatients x nTp x 11.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            expected_size = [2, 2, 11];  % 2 patients, 2 timepoints, 11 methods
            testCase.verifyEqual(size(results.volume_fractions), expected_size, ...
                'Volume fractions should be nPatients x nTp x nMethods.');
        end

        function testFDMatFx2NotFallback(testCase)
            % fDM at Fx2 (k=2) should NOT be flagged as fallback because
            % baseline data is available for comparison.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            fdm_idx = find(strcmp(results.method_names, 'fdm'));
            testCase.verifyFalse(results.fallback_flags(1, 2, fdm_idx), ...
                'fDM at Fx2 should NOT be a fallback (baseline exists).');
        end

        function testAllDiceCellArraySize(testCase)
            % all_dice cell array should match nPatients x nTp.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            testCase.verifyEqual(size(results.all_dice), [2, 2], ...
                'all_dice should be 2x2 cell array.');
        end

        function testEmptyPatientDiceIsEmpty(testCase)
            % Patient 2 has empty Fx2 data; its dice matrix should be empty.
            results = compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            testCase.verifyTrue(isempty(results.all_dice{2, 2}), ...
                'Empty patient Fx2 should have empty Dice matrix.');
        end

        function testDiaryFileCreated(testCase)
            % compare_core_methods should create a diary log file.
            compare_core_methods(testCase.DataVectors, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);
            diary off;

            diary_file = fullfile(testCase.TempDir, 'compare_core_methods_output_Standard.txt');
            testCase.verifyTrue(exist(diary_file, 'file') > 0, ...
                'Diary log file should be created.');
        end

    end
end
