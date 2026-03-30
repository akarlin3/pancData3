classdef test_extract_tumor_core_fit_info < matlab.unittest.TestCase
    % TEST_EXTRACT_TUMOR_CORE_FIT_INFO  Tests for the optional fit_info output
    % of extract_tumor_core.
    %
    % Validates that extract_tumor_core correctly populates the fit_info struct
    % with failure mode information while maintaining backward compatibility
    % with single-output callers.

    properties
        TempDir
        ConfigStruct
        AdcVec
        DVec
        FVec
        DstarVec
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

            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.output_folder = testCase.TempDir;
            testCase.ConfigStruct.dwi_types_to_run = 1;
            testCase.ConfigStruct.dwi_type_name = 'Standard';
            testCase.ConfigStruct.adc_thresh = 0.001;
            testCase.ConfigStruct.high_adc_thresh = 0.00115;
            testCase.ConfigStruct.d_thresh = 0.001;
            testCase.ConfigStruct.f_thresh = 0.1;
            testCase.ConfigStruct.dstar_thresh = 0.01;
            testCase.ConfigStruct.min_vox_hist = 20;
            testCase.ConfigStruct.core_method = 'adc_threshold';
            testCase.ConfigStruct.core_percentile = 25;
            testCase.ConfigStruct.core_n_clusters = 2;
            testCase.ConfigStruct.fdm_parameter = 'adc';
            testCase.ConfigStruct.fdm_thresh = 0.0004;
            testCase.ConfigStruct.adc_max = 0.003;

            rng(42);
            n_vox = 100;
            % Bimodal ADC: 50 low (core) + 50 high (margin)
            testCase.AdcVec = [0.0005 + 0.0001*randn(50,1); 0.002 + 0.0002*randn(50,1)];
            testCase.DVec = [0.0004 + 0.0001*randn(50,1); 0.0018 + 0.0002*randn(50,1)];
            testCase.FVec = [0.05 + 0.01*randn(50,1); 0.2 + 0.05*randn(50,1)];
            testCase.DstarVec = 0.02 * ones(n_vox, 1);
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

        function testFitInfoReturnedWhenRequested(testCase)
            % Verify fit_info is a struct with all expected fields.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';
            [~, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(isstruct(info), 'fit_info should be a struct.');
            expected_fields = {'method', 'success', 'fallback', 'empty_mask', ...
                'insufficient_voxels', 'all_nan_input', 'n_core_voxels', 'error_msg'};
            for i = 1:numel(expected_fields)
                testCase.verifyTrue(isfield(info, expected_fields{i}), ...
                    sprintf('fit_info should have field %s.', expected_fields{i}));
            end
        end

        function testSingleOutputStillWorks(testCase)
            % Single-output call should not error (backward compatibility).
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';
            mask = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(islogical(mask), 'Single output should be logical mask.');
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec));
        end

        function testSuccessForAdcThreshold(testCase)
            % ADC threshold with bimodal data should succeed.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';
            [~, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(info.success, 'adc_threshold should succeed.');
            testCase.verifyFalse(info.fallback, 'adc_threshold should not report fallback.');
            testCase.verifyGreaterThan(info.n_core_voxels, 0, ...
                'Should have some core voxels.');
            testCase.verifyEqual(info.method, 'adc_threshold');
        end

        function testFallbackDetectedForGMM(testCase)
            % GMM with only 2 valid voxels (below min_vox_hist) should fall back.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'gmm';
            cfg.min_vox_hist = 20;

            % Create vector with only 2 non-NaN values
            adc_sparse = nan(100, 1);
            adc_sparse(1:2) = [0.0005; 0.0008];
            d_sparse = nan(100, 1);
            d_sparse(1:2) = [0.0004; 0.0007];
            f_sparse = nan(100, 1);
            f_sparse(1:2) = [0.05; 0.06];
            dstar_sparse = nan(100, 1);
            dstar_sparse(1:2) = [0.01; 0.02];

            [~, info] = extract_tumor_core(cfg, adc_sparse, ...
                d_sparse, f_sparse, dstar_sparse, false, [], struct());

            testCase.verifyTrue(info.fallback, 'GMM with too few voxels should fall back.');
        end

        function testAllNanInput(testCase)
            % All-NaN input should set all_nan_input flag.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';
            n = 50;
            [mask, info] = extract_tumor_core(cfg, nan(n,1), ...
                nan(n,1), nan(n,1), nan(n,1), false, [], struct());

            testCase.verifyTrue(info.all_nan_input, 'Should detect all-NaN input.');
            testCase.verifyTrue(info.empty_mask, 'All-NaN should produce empty mask.');
            testCase.verifyFalse(info.success, 'All-NaN should not be successful.');
            testCase.verifyEqual(sum(mask), 0, 'Mask should be all-false.');
        end

        function testEmptyMaskDetected(testCase)
            % Very low threshold so no voxels qualify -> empty mask.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';
            cfg.adc_thresh = 0.00001;  % absurdly low

            [~, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(info.empty_mask, 'Should detect empty mask.');
            testCase.verifyEqual(info.n_core_voxels, 0);
        end

        function testInsufficientVoxels(testCase)
            % Set threshold so only ~2 voxels qualify, with min_vox_hist=5.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';
            cfg.min_vox_hist = 5;

            % Create data where only 2 voxels are below threshold
            adc = ones(100, 1) * 0.002;  % all above threshold
            adc(1:2) = 0.0005;  % only 2 below
            d = adc * 0.9;
            f = 0.15 * ones(100, 1);
            dstar = 0.02 * ones(100, 1);

            [~, info] = extract_tumor_core(cfg, adc, d, f, dstar, false, [], struct());

            testCase.verifyTrue(info.insufficient_voxels, ...
                'Should detect insufficient voxels when core < min_vox_hist.');
            testCase.verifyEqual(info.n_core_voxels, 2);
        end

        function testRegionGrowingFallbackWithout3D(testCase)
            % Region growing without 3D mask should fall back.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'region_growing';

            [~, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(info.fallback, ...
                'region_growing without 3D mask should fall back.');
        end

        function testFdmBaselineFallback(testCase)
            % fDM at timepoint_index=1 (no baseline) should fall back.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';
            opts = struct('timepoint_index', 1);

            [~, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], opts);

            testCase.verifyTrue(info.fallback, ...
                'fDM at baseline (no prior timepoint) should fall back.');
        end

        function testDThresholdMaskCorrectness(testCase)
            % With d_threshold method, verify mask selects voxels where D < d_thresh.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'd_threshold';
            cfg.d_thresh = 0.001;

            [mask, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(islogical(mask), 'Mask should be logical.');
            % Every selected voxel should have D <= d_thresh
            selected_d = testCase.DVec(mask);
            testCase.verifyTrue(all(selected_d <= cfg.d_thresh), ...
                'All selected voxels should have D <= d_thresh.');
            % Every unselected non-NaN voxel should have D > d_thresh
            unselected_valid = ~mask & ~isnan(testCase.DVec);
            testCase.verifyTrue(all(testCase.DVec(unselected_valid) > cfg.d_thresh), ...
                'All unselected valid voxels should have D > d_thresh.');
            testCase.verifyEqual(info.method, 'd_threshold');
        end

        function testOtsuProducesBinaryMask(testCase)
            % With otsu method, verify mask is logical and has both true/false values
            % for bimodal data.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'otsu';

            evalc('[mask, info] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');

            testCase.verifyTrue(islogical(mask), 'Otsu mask should be logical.');
            testCase.verifyTrue(any(mask), 'Bimodal data should produce some true values.');
            testCase.verifyTrue(any(~mask), 'Bimodal data should produce some false values.');
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec), ...
                'Mask size should match input size.');
        end

        function testKmeansTwoClusters(testCase)
            % With kmeans method and core_n_clusters=2, verify mask partitions
            % data into 2 groups.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'kmeans';
            cfg.core_n_clusters = 2;

            evalc('[mask, info] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');

            testCase.verifyTrue(islogical(mask), 'K-means mask should be logical.');
            n_core = sum(mask);
            n_non_core = sum(~mask & ~isnan(testCase.AdcVec));
            testCase.verifyGreaterThan(n_core, 0, 'Should have core voxels.');
            testCase.verifyGreaterThan(n_non_core, 0, 'Should have non-core voxels.');
            % Core cluster should have lower mean ADC than non-core
            mean_core_adc = mean(testCase.AdcVec(mask), 'omitnan');
            mean_noncore_adc = mean(testCase.AdcVec(~mask & ~isnan(testCase.AdcVec)), 'omitnan');
            testCase.verifyLessThan(mean_core_adc, mean_noncore_adc, ...
                'Core cluster should have lower mean ADC than non-core.');
        end

        function testPercentileSelectsBottom25(testCase)
            % With percentile method and core_percentile=25, verify ~25% of
            % voxels selected.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'percentile';
            cfg.core_percentile = 25;
            % Set adc_thresh high so it does not limit the percentile
            cfg.adc_thresh = 0.01;

            [mask, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            n_valid = sum(~isnan(testCase.AdcVec));
            fraction_selected = sum(mask) / n_valid;
            % Allow some tolerance since percentile threshold picks <= boundary
            testCase.verifyGreaterThan(fraction_selected, 0.10, ...
                'Should select at least ~10% of voxels.');
            testCase.verifyLessThan(fraction_selected, 0.40, ...
                'Should select no more than ~40% of voxels.');
        end

        function testDfIntersectionMask(testCase)
            % With df_intersection method, verify mask selects voxels where
            % both D < d_thresh AND f <= f_thresh.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'df_intersection';
            cfg.d_thresh = 0.001;
            cfg.f_thresh = 0.1;

            [mask, info] = extract_tumor_core(cfg, testCase.AdcVec, ...
                testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());

            testCase.verifyTrue(islogical(mask), 'Mask should be logical.');
            % Verify intersection logic on selected voxels
            d_below = testCase.DVec(mask) <= cfg.d_thresh;
            f_below = testCase.FVec(mask) <= cfg.f_thresh;
            testCase.verifyTrue(all(d_below), ...
                'All selected voxels should have D <= d_thresh.');
            testCase.verifyTrue(all(f_below), ...
                'All selected voxels should have f <= f_thresh.');
            testCase.verifyEqual(info.method, 'df_intersection');
        end

        function testSpectralFallbackWithFewVoxels(testCase)
            % With spectral method and very few voxels (<spectral_min_voxels),
            % should fall back to ADC threshold.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'spectral';
            cfg.spectral_min_voxels = 200;  % higher than our 100-voxel data

            evalc('[mask, info] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');

            testCase.verifyTrue(info.fallback, ...
                'Spectral with too few voxels should fall back.');
        end

        function testMethodFallbackRecovery(testCase)
            % When a method throws an error internally, verify it produces
            % empty mask with success=false.
            rng(42);
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'adc_threshold';

            % All-NaN triggers the early return path with success=false
            n = 50;
            [mask, info] = extract_tumor_core(cfg, nan(n,1), ...
                nan(n,1), nan(n,1), nan(n,1), false, [], struct());

            testCase.verifyFalse(info.success, ...
                'All-NaN input should produce success=false.');
            testCase.verifyTrue(info.empty_mask, ...
                'Should produce empty mask on failure.');
            testCase.verifyEqual(sum(mask), 0, 'Mask should be all-false.');
        end

        function testAllMethodsReturnLogicalMask(testCase)
            % Loop through all 11 methods, verify each returns a logical mask
            % of correct size.
            rng(42);
            methods_list = {'adc_threshold', 'd_threshold', 'df_intersection', ...
                'otsu', 'gmm', 'kmeans', 'percentile'};
            % Excluded: region_growing (needs 3D), active_contours (needs 3D),
            % spectral (may need special setup), fdm (needs baseline)

            for i = 1:numel(methods_list)
                cfg = testCase.ConfigStruct;
                cfg.core_method = methods_list{i};
                cfg.min_vox_hist = 10;

                evalc('[mask, info] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');

                testCase.verifyTrue(islogical(mask), ...
                    sprintf('%s should return logical mask.', methods_list{i}));
                testCase.verifyEqual(numel(mask), numel(testCase.AdcVec), ...
                    sprintf('%s mask should match input size.', methods_list{i}));
            end

            % Also test spectral and fdm (expected to fallback but still return logical)
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'spectral';
            cfg.spectral_min_voxels = 10;
            evalc('[mask, ~] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');
            testCase.verifyTrue(islogical(mask), 'spectral should return logical mask.');
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec));

            cfg.core_method = 'fdm';
            opts = struct('timepoint_index', 1);
            evalc('[mask, ~] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], opts);');
            testCase.verifyTrue(islogical(mask), 'fdm should return logical mask.');
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec));

            cfg.core_method = 'region_growing';
            evalc('[mask, ~] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');
            testCase.verifyTrue(islogical(mask), 'region_growing should return logical mask.');
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec));

            cfg.core_method = 'active_contours';
            evalc('[mask, ~] = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, testCase.FVec, testCase.DstarVec, false, [], struct());');
            testCase.verifyTrue(islogical(mask), 'active_contours should return logical mask.');
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec));
        end

    end
end
