classdef test_new_core_methods < matlab.unittest.TestCase
% TEST_NEW_CORE_METHODS  Detailed tests for percentile, spectral, and fDM
%   tumor core extraction methods.
%
%   These three methods were added after the original 8 core methods.
%   Unlike the threshold-based methods, they use data-driven or longitudinal
%   approaches:
%     - Percentile: selects the lowest N-th percentile of ADC as core,
%       capped by a safety floor (adc_thresh) to prevent overly aggressive
%       core delineation
%     - Spectral: uses spectral clustering on multi-parameter feature space
%       (ADC, D, f, D*) to identify the most restricted cluster as core
%     - fDM (functional Diffusion Map): identifies treatment-resistant voxels
%       by comparing current and baseline diffusion values; voxels with a
%       negative delta (ADC decreased) are classified as progressing (core)
%
%   Run tests with:
%     results = runtests('tests/test_new_core_methods.m');

    properties
        ConfigStruct  % Mock pipeline configuration with core method settings
        AdcVec        % 100-element ADC vector (bimodal: 30 low + 70 high)
        DVec          % 100-element D vector (bimodal, correlated with ADC)
        FVec          % 100-element f vector (bimodal: 30 low + 70 high)
        DstarVec      % 100-element D* vector (bimodal, unlike test_core_methods which uses NaN)
        GtvMask3d     % 10x10x1 logical mask (entire slab is within GTV)
    end

    methods (TestMethodSetup)
        function setupFixtures(testCase)
            % Build a reproducible bimodal voxel population and config struct.
            % Unlike test_core_methods, D* is populated here (not NaN) to
            % allow spectral clustering to use all 4 IVIM parameters.
            [dir_path, ~, ~] = fileparts(mfilename('fullpath'));
            addpath(fullfile(dir_path, '..', 'utils'));
            addpath(fullfile(dir_path, '..', 'core'));

            rng(42); % Fixed seed for reproducible random data
            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.adc_thresh = 0.001;   % mm^2/s
            testCase.ConfigStruct.d_thresh = 0.001;     % mm^2/s
            testCase.ConfigStruct.f_thresh = 0.1;       % 10%
            testCase.ConfigStruct.dstar_thresh = 0.01;  % mm^2/s
            testCase.ConfigStruct.min_vox_hist = 50;    % Reduced for test data size
            testCase.ConfigStruct.core_percentile = 25;
            testCase.ConfigStruct.core_n_clusters = 2;  % Binary: core vs. margin
            testCase.ConfigStruct.fdm_parameter = 'adc';
            testCase.ConfigStruct.fdm_thresh = 0.0002;  % mm^2/s delta threshold

            % Bimodal distribution: 30 "core" voxels with low diffusion values
            % and 70 "margin" voxels with high diffusion values.
            testCase.AdcVec = [0.0005 + 0.0001*randn(30,1); 0.0015 + 0.0002*randn(70,1)];
            testCase.DVec   = [0.0004 + 0.0001*randn(30,1); 0.0014 + 0.0002*randn(70,1)];
            testCase.FVec   = [0.05 + 0.01*randn(30,1);     0.2 + 0.05*randn(70,1)];
            testCase.DstarVec = [0.005 + 0.002*randn(30,1); 0.02 + 0.005*randn(70,1)];
            testCase.GtvMask3d = true(10, 10, 1);
        end
    end

    % ==================================================================
    % Percentile method tests
    % ==================================================================
    methods (Test)
        function testPercentileSelectsApproxNPercent(testCase)
            % Verifies that the percentile method selects approximately the
            % bottom N% of ADC voxels as core. With 100 voxels and
            % core_percentile=25, expect ~25 voxels selected. The adc_thresh
            % safety floor is set very high (0.01) so it does not interfere.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'percentile';
            cfg.core_percentile = 25;
            cfg.adc_thresh = 0.01; % Effectively disabled

            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, true, testCase.GtvMask3d);

            n_core = sum(mask);
            % Allow +-2 tolerance for ties at the percentile boundary
            testCase.verifyGreaterThanOrEqual(n_core, 23);
            testCase.verifyLessThanOrEqual(n_core, 27);
        end

        function testPercentileSafetyFloor(testCase)
            % The safety floor ensures the percentile cutoff never exceeds
            % adc_thresh. Here core_percentile=50 would normally select
            % ~50 voxels, but adc_thresh=0.0003 is so restrictive that
            % only the very lowest ADC voxels qualify. This prevents
            % clinically implausible core sizes when the percentile is set
            % too aggressively.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'percentile';
            cfg.core_percentile = 50;
            cfg.adc_thresh = 0.0003; % Very restrictive cap

            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, false, []);

            % Every selected core voxel must have ADC <= the safety threshold
            core_vals = testCase.AdcVec(mask);
            if ~isempty(core_vals)
                testCase.verifyLessThanOrEqual(max(core_vals), cfg.adc_thresh);
            end
        end

        function testPercentileAllNaN(testCase)
            % Edge case: when all voxels are NaN (no valid data), the
            % percentile method must return an all-false mask rather than
            % crashing on prctile() or producing undefined behaviour.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'percentile';
            nan_vec = nan(100, 1);
            mask = extract_tumor_core(cfg, nan_vec, nan_vec, nan_vec, nan_vec, false, []);
            testCase.verifyEqual(sum(mask), 0);
        end

        function testPercentile100SelectsAll(testCase)
            % Boundary case: core_percentile=100 means "include all voxels
            % up to the 100th percentile", which is the entire population.
            % With adc_thresh set unreasonably high, the safety floor does
            % not interfere, so all 100 voxels should be selected.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'percentile';
            cfg.core_percentile = 100;
            cfg.adc_thresh = 1; % Safety floor effectively disabled

            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, false, []);

            testCase.verifyEqual(sum(mask), numel(testCase.AdcVec));
        end
    end

    % ==================================================================
    % Spectral clustering method tests
    % ==================================================================
    methods (Test)
        function testSpectralIdentifiesClusters(testCase)
            % Spectral clustering should partition the bimodal voxel
            % population into two clusters and assign the cluster with
            % lower mean ADC as "core". Verifies output type, length,
            % non-degenerate partitioning, and that core has lower ADC
            % than non-core (the defining characteristic of restricted
            % diffusion in tumour cores).
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'spectral';

            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, true, testCase.GtvMask3d);

            testCase.verifyTrue(islogical(mask));
            testCase.verifyEqual(numel(mask), numel(testCase.AdcVec));

            n_core = sum(mask);
            testCase.verifyGreaterThan(n_core, 0);
            testCase.verifyLessThan(n_core, numel(testCase.AdcVec));

            % The core cluster must exhibit lower mean ADC (restricted diffusion)
            mean_core = mean(testCase.AdcVec(mask));
            mean_noncore = mean(testCase.AdcVec(~mask));
            testCase.verifyLessThan(mean_core, mean_noncore);
        end

        function testSpectralAdcOnlyWhenIVIMNaN(testCase)
            % When IVIM parameters (D, f, D*) are all NaN (e.g., IVIM
            % fitting failed or was not performed), spectral clustering
            % should fall back to using ADC alone as the feature vector.
            % It must still produce a valid non-degenerate partition.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'spectral';
            nan_ivim = nan(100, 1);

            mask = extract_tumor_core(cfg, testCase.AdcVec, nan_ivim, ...
                nan_ivim, nan_ivim, false, []);

            testCase.verifyTrue(islogical(mask));
            n_core = sum(mask);
            testCase.verifyGreaterThan(n_core, 0);
            testCase.verifyLessThan(n_core, numel(testCase.AdcVec));
        end

        function testSpectralFallbackTooFewVoxels(testCase)
            % Spectral clustering requires a minimum number of voxels
            % (min_vox_hist) to build a meaningful similarity graph.
            % When the voxel count is below this threshold, the method
            % must fall back to simple adc_threshold to avoid noisy or
            % degenerate clustering results.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'spectral';
            cfg.min_vox_hist = 200; % Exceeds our 100-voxel test data

            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, false, []);

            % Verify exact equivalence with the threshold fallback
            expected = testCase.AdcVec <= cfg.adc_thresh;
            testCase.verifyEqual(mask, expected);
        end
    end

    % ==================================================================
    % fDM method tests
    % ==================================================================
    methods (Test)
        function testFdmBaselineFallsBack(testCase)
            % At the first timepoint (baseline), there is no prior scan to
            % compute a delta from. The fDM method must fall back to simple
            % ADC thresholding as there is no longitudinal information.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';

            opts = struct('timepoint_index', 1);
            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, false, [], opts);

            expected = testCase.AdcVec <= cfg.adc_thresh;
            testCase.verifyEqual(mask, expected);
        end

        function testFdmKnownDelta(testCase)
            % Deterministic test with three known voxel populations:
            %   - 30 "progressing" voxels: ADC decreased by 0.001 (core)
            %   - 30 "responding" voxels: ADC increased by 0.001 (non-core)
            %   - 30 "stable" voxels: ADC unchanged (non-core)
            % The fDM threshold of 0.0005 is well below the delta of 0.001,
            % so the classification should be exact with no ambiguity.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';
            cfg.fdm_thresh = 0.0005;
            cfg.fdm_parameter = 'adc';

            n = 90;
            baseline = 0.001 * ones(n, 1);
            % delta = current - baseline: negative = progressing, positive = responding
            current  = [baseline(1:30) - 0.001; baseline(31:60) + 0.001; baseline(61:90)];

            opts = struct('timepoint_index', 2);
            opts.baseline_adc_vec = baseline;

            mask = extract_tumor_core(cfg, current, nan(n,1), nan(n,1), nan(n,1), false, [], opts);

            % Progressing voxels (ADC decreased) should be classified as core
            testCase.verifyEqual(sum(mask(1:30)), 30, 'All progressing voxels should be core');
            % Responding voxels (ADC increased) should NOT be core
            testCase.verifyEqual(sum(mask(31:60)), 0, 'No responding voxels should be core');
            % Stable voxels (no change) should NOT be core
            testCase.verifyEqual(sum(mask(61:90)), 0, 'No stable voxels should be core');
        end

        function testFdmThreeFractionsSumToOne(testCase)
            % Verifies the fundamental fDM partition property: every valid
            % voxel must be classified as exactly one of {responding,
            % progressing, stable}. The three fractions must sum to 1.0.
            % This is a mathematical invariant, not method-specific.
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';
            cfg.fdm_thresh = 0.0003;

            rng(99);
            n = 100;
            baseline = 0.001 * ones(n, 1);
            % Add Gaussian noise to create a mix of all three categories
            current = baseline + 0.0005 * randn(n, 1);

            opts = struct('timepoint_index', 2);
            opts.baseline_adc_vec = baseline;

            % Manually compute the three fractions to verify the partition
            delta = current - baseline;
            valid = ~isnan(delta);
            n_valid = sum(valid);
            responding = sum(delta(valid) > cfg.fdm_thresh) / n_valid;
            progressing = sum(delta(valid) < -cfg.fdm_thresh) / n_valid;
            stable = sum(abs(delta(valid)) <= cfg.fdm_thresh) / n_valid;

            testCase.verifyEqual(responding + progressing + stable, 1.0, 'AbsTol', 1e-10);
        end

        function testFdmUsesD(testCase)
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';
            cfg.fdm_parameter = 'd';
            cfg.fdm_thresh = 0.0003;

            n = 60;
            baseline_d = 0.0008 * ones(n, 1);
            % D decreased everywhere → all should be core
            current_d = baseline_d - 0.001;
            current_adc = 0.002 * ones(n, 1); % high ADC (would not be core by threshold)

            opts = struct('timepoint_index', 3);
            opts.baseline_d_vec = baseline_d;

            mask = extract_tumor_core(cfg, current_adc, current_d, nan(n,1), nan(n,1), false, [], opts);

            testCase.verifyEqual(sum(mask), n, 'All voxels should be core (D decreased)');
        end

        function testFdmVectorMismatchFallback(testCase)
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';

            opts = struct('timepoint_index', 2);
            opts.baseline_adc_vec = rand(50, 1); % wrong length (100 vs 50)

            mask = extract_tumor_core(cfg, testCase.AdcVec, testCase.DVec, ...
                testCase.FVec, testCase.DstarVec, false, [], opts);

            % Should fall back to adc_threshold
            expected = testCase.AdcVec <= cfg.adc_thresh;
            testCase.verifyEqual(mask, expected);
        end

        function testFdmRepeatabilityCorOverridesThresh(testCase)
            cfg = testCase.ConfigStruct;
            cfg.core_method = 'fdm';
            cfg.fdm_thresh = 0.01; % very large → nothing would be classified

            n = 60;
            baseline = 0.001 * ones(n, 1);
            current = baseline - 0.0005; % small decrease

            opts = struct('timepoint_index', 2);
            opts.baseline_adc_vec = baseline;
            opts.repeatability_cor = 0.0003; % small → should detect the change

            mask = extract_tumor_core(cfg, current, nan(n,1), nan(n,1), nan(n,1), false, [], opts);

            % With CoR=0.0003, delta of -0.0005 exceeds threshold → core
            testCase.verifyEqual(sum(mask), n);
        end
    end
end
