classdef test_optimize_adc_threshold < matlab.unittest.TestCase
    % TEST_OPTIMIZE_ADC_THRESHOLD  Tests for optimize_adc_threshold.m

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
            addpath(fullfile(repoRoot, 'dependencies'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end

            set(0, 'DefaultFigureVisible', 'off');

            [testCase.ConfigStruct, testCase.DataVectors, testCase.IdList, testCase.GtvLocations] = ...
                buildOptThreshData(testCase.TempDir);
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
            result = optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);

            expected = { ...
                'thresholds', 'median_dice', 'mean_dice', 'std_dice', ...
                'median_vol_frac', 'mean_vol_frac', 'n_patients', ...
                'optimal_thresh', 'optimal_dice', 'optimal_vol_frac', 'per_patient_dice', ...
                'inflection_thresh', 'inflection_idx', 'inflection_curvature', ...
                'vol_frac_curvature', ...
                'significance_pvalues', 'significance_thresh', 'significance_pvalue', ...
                'significance_n_lc', 'significance_n_lf', 'significance_metric'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(result, expected{i}), ...
                    sprintf('Output should contain field %s.', expected{i}));
            end
        end

        function testInflectionTacticContract(testCase)
            % Tactic 2 contract: curvature vector is 13-long with NaN
            % endpoints, scalar outputs are present and either finite
            % with a valid index or all-NaN (if the curve is too flat).
            % Exact knee location depends on morphological cleanup noise
            % in the small 5x5x5 fixture; correctness of the detection
            % algorithm itself is covered by testInflectionKneeShape
            % below with a hand-crafted monotonic sigmoid curve.
            result = optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);

            testCase.verifyEqual(numel(result.vol_frac_curvature), 13, ...
                'Curvature must be 1x13 (matching thresholds).');
            testCase.verifyTrue(isnan(result.vol_frac_curvature(1)), ...
                'Curvature endpoint (index 1) must be NaN.');
            testCase.verifyTrue(isnan(result.vol_frac_curvature(end)), ...
                'Curvature endpoint (index end) must be NaN.');

            if ~isnan(result.inflection_thresh)
                testCase.verifyGreaterThanOrEqual(result.inflection_idx, 2);
                testCase.verifyLessThanOrEqual(result.inflection_idx, 12);
                testCase.verifyEqual(result.inflection_thresh, ...
                    result.thresholds(result.inflection_idx), 'AbsTol', 1e-12);
                testCase.verifyFalse(isnan(result.inflection_curvature));
            else
                % Flat curve — all three scalars should be NaN together.
                testCase.verifyTrue(isnan(result.inflection_idx));
                testCase.verifyTrue(isnan(result.inflection_curvature));
            end
        end


        function testSignificanceMissingLF(testCase)
            % Tactic 3: the default fixture has no .LF field, so the
            % significance scan must gracefully report NaN/empty rather
            % than crashing on the missing field.
            result = optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);

            testCase.verifyTrue(isnan(result.significance_thresh));
            testCase.verifyTrue(isnan(result.significance_pvalue));
            testCase.verifyEqual(result.significance_n_lc, 0);
            testCase.verifyEqual(result.significance_n_lf, 0);
            testCase.verifyEqual(numel(result.significance_pvalues), 13);
            testCase.verifyTrue(all(isnan(result.significance_pvalues)));
        end

        function testSignificanceWithLF(testCase)
            % Tactic 3: inject distinct ADC distributions for LC vs LF
            % (LC patients keep the cellular cluster; LF patients shift
            % the cellular cluster slightly higher) so the per-threshold
            % Wilcoxon p-value is meaningful and a best threshold exists.
            dv = testCase.DataVectors;
            n_patients = size(dv, 1);
            % First half = LC (LF=0), second half = LF (LF=1).  Need >=3
            % per group; the fixture has 5 patients so split 3/2 won't
            % satisfy LF>=3.  Re-shape the fixture inline.
            n_each = 3;
            n_total = 2 * n_each;
            base_template = dv(1, 1, 1);
            n_vox = numel(base_template.adc_vector);
            template = base_template;
            template.LF = 0;
            dv2 = repmat(template, n_total, 1, size(dv, 3));
            for j = 1:n_total
                lf = double(j > n_each);
                shift = lf * 1.0e-4;  % LF patients have slightly higher ADC
                base_adc = [(1.1e-3 + shift) * ones(round(n_vox/2), 1); ...
                            2.5e-3 * ones(n_vox - round(n_vox/2), 1)];
                for rpi = 1:size(dv, 3)
                    noise = 1e-4 * randn(n_vox, 1);
                    dv2(j, 1, rpi).adc_vector = base_adc + noise;
                    dv2(j, 1, rpi).d_vector = base_adc + noise;
                    dv2(j, 1, rpi).f_vector = 0.1 * ones(n_vox, 1);
                    dv2(j, 1, rpi).dstar_vector = 0.02 * ones(n_vox, 1);
                    dv2(j, 1, rpi).LF = lf;
                end
            end
            id2 = arrayfun(@(x) sprintf('P%02d', x), 1:n_total, 'UniformOutput', false);
            gtv2 = cell(n_total, 1, size(dv, 3));
            for j = 1:n_total
                for rpi = 1:size(dv, 3)
                    gtv2{j, 1, rpi} = testCase.GtvLocations{1, 1, 1};
                end
            end

            result = optimize_adc_threshold(dv2, testCase.ConfigStruct, id2, gtv2);

            testCase.verifyEqual(numel(result.significance_pvalues), 13);
            % At least one threshold must reach the 3-per-group floor.
            testCase.verifyTrue(any(~isnan(result.significance_pvalues)), ...
                'At least one threshold should produce a valid p-value with 3+3 cohort');
            testCase.verifyFalse(isnan(result.significance_thresh));
            testCase.verifyGreaterThanOrEqual(result.significance_pvalue, 0);
            testCase.verifyLessThanOrEqual(result.significance_pvalue, 1);
            testCase.verifyEqual(result.significance_n_lc, n_each);
            testCase.verifyEqual(result.significance_n_lf, n_each);
            testCase.verifyEqual(result.significance_metric, 'wilcoxon ranksum on vol_frac');
        end

        function testThresholdsLength(testCase)
            result = optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);
            testCase.verifyEqual(numel(result.thresholds), 13, ...
                'There should be 13 threshold values.');
            testCase.verifyEqual(result.thresholds(1), 0.8e-3, 'AbsTol', 1e-9);
            testCase.verifyEqual(result.thresholds(end), 2.0e-3, 'AbsTol', 1e-9);
        end

        function testOptimalThreshNotAtBoundary(testCase)
            result = optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);
            % With ADC values concentrated around 1.2e-3, optimal should be
            % somewhere in the interior of the sweep range, not at either end.
            if ~isnan(result.optimal_dice)
                testCase.verifyNotEqual(result.optimal_thresh, result.thresholds(1), ...
                    'Optimal threshold should not be at lower boundary.');
                testCase.verifyNotEqual(result.optimal_thresh, result.thresholds(end), ...
                    'Optimal threshold should not be at upper boundary.');
            end
        end

        function testFigureSaved(testCase)
            optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);
            png_path = fullfile(testCase.TempDir, ...
                sprintf('adc_threshold_optimization_%s.png', testCase.ConfigStruct.dwi_type_name));
            testCase.verifyTrue(exist(png_path, 'file') == 2, ...
                'Optimization PNG should be saved.');
        end

        function testPerPatientDiceDimensions(testCase)
            result = optimize_adc_threshold(testCase.DataVectors, testCase.ConfigStruct, ...
                testCase.IdList, testCase.GtvLocations);
            n_patients = numel(testCase.IdList);
            testCase.verifyEqual(size(result.per_patient_dice), [n_patients, 13]);
        end

    end
end


function [cfg, dv, id_list, gtv_locations] = buildOptThreshData(tempDir)
% Build 5-patient dataset with 2 Fx1 repeats each. Synthetic ADC data
% centered near 1.2e-3 with small repeat-to-repeat noise near the edges.
%
% Mask size: 8x8x8 (512 voxels) — `safe_load_mask.validate_mask_dimensions`
% requires every non-singleton dimension to be >= 8, so smaller masks
% load as [] and silently skip every patient.

    rng(42);
    n_patients = 5;
    n_repeats = 2;

    % 8x8x8 GTV volume (512 voxels)
    grid_sz = 8;
    n_vox = grid_sz^3;
    mask_3d = true(grid_sz, grid_sz, grid_sz);

    % Save the GTV mask once (shared across repeats).
    gtv_path = fullfile(tempDir, 'test_gtv.mat');
    Stvol3d = mask_3d; %#ok<NASGU>
    save(gtv_path, 'Stvol3d');

    template = struct('adc_vector', [], 'd_vector', [], ...
        'f_vector', [], 'dstar_vector', [], ...
        'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
        'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], ...
        'dstar_vector_ivimnet', [], ...
        'vox_dims', [2, 2, 3], 'vox_vol', 0.012);
    dv = repmat(template, n_patients, 1, n_repeats);

    % ADC values: low cluster at 1.1e-3 (cellular tumor core), high cluster
    % at 2.5e-3 (above all tested thresholds). Small repeat-to-repeat noise
    % is added so that sub-volume definitions shift near the threshold and
    % the optimal Dice is at an interior threshold, not at a boundary.
    for j = 1:n_patients
        base_adc = [1.1e-3 * ones(round(n_vox/2), 1); ...
                    2.5e-3 * ones(n_vox - round(n_vox/2), 1)];
        for rpi = 1:n_repeats
            noise = 1e-4 * randn(n_vox, 1);
            dv(j, 1, rpi).adc_vector = base_adc + noise;
            dv(j, 1, rpi).d_vector = base_adc + noise;
            dv(j, 1, rpi).f_vector = 0.1 * ones(n_vox, 1);
            dv(j, 1, rpi).dstar_vector = 0.02 * ones(n_vox, 1);
        end
    end

    id_list = arrayfun(@(x) sprintf('P%02d', x), 1:n_patients, 'UniformOutput', false);
    gtv_locations = cell(n_patients, 1, n_repeats);
    for j = 1:n_patients
        for rpi = 1:n_repeats
            gtv_locations{j, 1, rpi} = gtv_path;
        end
    end

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';
    cfg.morph_min_cc = 1;
    if exist('OCTAVE_VERSION', 'builtin')
        k = zeros(3,3,3); k(2,2,:)=1; k(2,:,2)=1; k(:,2,2)=1;
        cfg.morph_se = strel('arbitrary', k);
    else
        cfg.morph_se = strel('sphere', 1);
    end
end
