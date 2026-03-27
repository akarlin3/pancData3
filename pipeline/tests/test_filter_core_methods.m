classdef test_filter_core_methods < matlab.unittest.TestCase
    % TEST_FILTER_CORE_METHODS  Tests for the core method pruning function.
    %
    % Validates that filter_core_methods correctly:
    %   - Retains all methods when threshold is 1.0
    %   - Prunes methods exceeding failure rate threshold
    %   - Handles manual exclusions
    %   - Never prunes adc_threshold (safety guard)
    %   - Produces correct pruned_info structure
    %   - Handles min_core_voxels filter
    %   - Applies combined filters correctly

    properties
        AllMethods
        FailureTable
        ConfigStruct
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end

            testCase.AllMethods = {'adc_threshold', 'd_threshold', 'df_intersection', ...
                'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
                'percentile', 'spectral', 'fdm'};

            testCase.FailureTable = buildFailureTable();

            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.max_core_failure_rate = 1.0;
            testCase.ConfigStruct.excluded_core_methods = {};
            testCase.ConfigStruct.min_core_voxels = 0;
        end
    end

    methods(Test)

        function testNoPruningWhenThresholdIsOne(testCase)
            % All 11 methods should be retained when threshold is 1.0.
            cfg = testCase.ConfigStruct;
            cfg.max_core_failure_rate = 1.0;

            [active, pruned] = filter_core_methods(testCase.AllMethods, ...
                testCase.FailureTable, cfg);

            testCase.verifyEqual(numel(active), 11, ...
                'All 11 methods should be retained with threshold 1.0.');
            testCase.verifyEmpty(pruned, 'No methods should be pruned.');
        end

        function testPrunesHighFailureMethods(testCase)
            % Methods with >25% failure should be pruned.
            cfg = testCase.ConfigStruct;
            cfg.max_core_failure_rate = 0.25;

            [active, pruned] = filter_core_methods(testCase.AllMethods, ...
                testCase.FailureTable, cfg);

            % gmm (30%), active_contours (50%), spectral (40%) should be pruned
            testCase.verifyFalse(ismember('gmm', active), 'gmm should be pruned.');
            testCase.verifyFalse(ismember('active_contours', active), 'active_contours should be pruned.');
            testCase.verifyFalse(ismember('spectral', active), 'spectral should be pruned.');
            testCase.verifyTrue(ismember('adc_threshold', active), 'adc_threshold should be retained.');
            testCase.verifyTrue(ismember('d_threshold', active), 'd_threshold should be retained.');
            testCase.verifyEqual(numel(pruned), 3, 'Should prune exactly 3 methods.');
        end

        function testManualExclusion(testCase)
            % fdm should be pruned via manual exclusion even if rate is low.
            cfg = testCase.ConfigStruct;
            cfg.excluded_core_methods = {'fdm'};

            [active, pruned] = filter_core_methods(testCase.AllMethods, ...
                testCase.FailureTable, cfg);

            testCase.verifyFalse(ismember('fdm', active), 'fdm should be manually excluded.');
            testCase.verifyEqual(numel(pruned), 1);
            testCase.verifyEqual(pruned(1).name, 'fdm');
            testCase.verifyEqual(pruned(1).reason, 'manual_exclusion');
        end

        function testAdcThresholdNeverPruned(testCase)
            % adc_threshold should never be pruned even with 100% failure rate.
            ft = testCase.FailureTable;
            % Artificially set adc_threshold to 100% failure
            ft.any_failure_rate(1, :) = [1.0 1.0 1.0];

            cfg = testCase.ConfigStruct;
            cfg.max_core_failure_rate = 0.01;

            [active, ~] = filter_core_methods(testCase.AllMethods, ft, cfg);

            testCase.verifyTrue(ismember('adc_threshold', active), ...
                'adc_threshold must never be pruned.');
        end

        function testPrunedInfoStructure(testCase)
            % Verify each entry in pruned_info has expected fields.
            cfg = testCase.ConfigStruct;
            cfg.max_core_failure_rate = 0.25;

            [~, pruned] = filter_core_methods(testCase.AllMethods, ...
                testCase.FailureTable, cfg);

            testCase.verifyGreaterThan(numel(pruned), 0, 'Should have pruned methods.');
            for i = 1:numel(pruned)
                testCase.verifyTrue(isfield(pruned(i), 'name'));
                testCase.verifyTrue(isfield(pruned(i), 'reason'));
                testCase.verifyTrue(isfield(pruned(i), 'failure_rate'));
                testCase.verifyTrue(isfield(pruned(i), 'pipeline'));
                testCase.verifyTrue(ischar(pruned(i).name));
                testCase.verifyTrue(ischar(pruned(i).reason));
            end
        end

        function testAllMethodsPrunedExceptAdc(testCase)
            % Very low threshold should prune everything except adc_threshold.
            cfg = testCase.ConfigStruct;
            cfg.max_core_failure_rate = 0.01;

            [active, pruned] = filter_core_methods(testCase.AllMethods, ...
                testCase.FailureTable, cfg);

            testCase.verifyEqual(numel(active), 1, 'Only adc_threshold should remain.');
            testCase.verifyEqual(active{1}, 'adc_threshold');
            testCase.verifyEqual(numel(pruned), 10, 'Should prune 10 methods.');
        end

        function testMinCoreVoxelsPruning(testCase)
            % Method with low median voxels should be pruned.
            ft = testCase.FailureTable;
            % Set spectral to have very low median voxels
            spectral_idx = find(strcmp(ft.method_names, 'spectral'));
            ft.median_core_voxels(spectral_idx, :) = [8 5 3];

            cfg = testCase.ConfigStruct;
            cfg.min_core_voxels = 50;

            [active, pruned] = filter_core_methods(testCase.AllMethods, ft, cfg);

            testCase.verifyFalse(ismember('spectral', active), ...
                'spectral with low voxel count should be pruned.');
            % Find spectral in pruned
            spectral_pruned = false;
            for i = 1:numel(pruned)
                if strcmp(pruned(i).name, 'spectral')
                    testCase.verifyEqual(pruned(i).reason, 'insufficient_voxels');
                    spectral_pruned = true;
                end
            end
            testCase.verifyTrue(spectral_pruned, 'spectral should be in pruned_info.');
        end

        function testCombinedFilters(testCase)
            % Manual exclusion + failure rate + min voxels all active.
            ft = testCase.FailureTable;
            percentile_idx = find(strcmp(ft.method_names, 'percentile'));
            ft.median_core_voxels(percentile_idx, :) = [3 2 1];

            cfg = testCase.ConfigStruct;
            cfg.excluded_core_methods = {'fdm'};
            cfg.max_core_failure_rate = 0.35;
            cfg.min_core_voxels = 10;

            [active, pruned] = filter_core_methods(testCase.AllMethods, ft, cfg);

            % fdm: manual exclusion
            testCase.verifyFalse(ismember('fdm', active));
            % active_contours (50%), spectral (40%): failure rate
            testCase.verifyFalse(ismember('active_contours', active));
            testCase.verifyFalse(ismember('spectral', active));
            % percentile: insufficient voxels
            testCase.verifyFalse(ismember('percentile', active));
            % adc_threshold must remain
            testCase.verifyTrue(ismember('adc_threshold', active));

            % Verify reasons in pruned_info
            reasons = {pruned.reason};
            testCase.verifyTrue(any(strcmp(reasons, 'manual_exclusion')));
            testCase.verifyTrue(any(strcmp(reasons, 'failure_rate')));
            testCase.verifyTrue(any(strcmp(reasons, 'insufficient_voxels')));
        end

    end
end


function ft = buildFailureTable()
% Build synthetic failure_table with known rates.
    ft.method_names = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
    ft.pipeline_names = {'Standard', 'DnCNN', 'IVIMNet'};
    ft.n_patients = 20;
    ft.n_timepoints = 3;

    % any_failure_rate: [11 x 3]
    %  adc_threshold: 0%, d_threshold: 5%, df_intersection: 10%,
    %  otsu: 12%, gmm: 30%, kmeans: 15%, region_growing: 20%,
    %  active_contours: 50%, percentile: 8%, spectral: 40%, fdm: 10%
    ft.any_failure_rate = [
        0.00 0.00 0.00;   % adc_threshold
        0.05 0.05 0.06;   % d_threshold
        0.10 0.10 0.12;   % df_intersection
        0.12 0.12 0.14;   % otsu
        0.28 0.30 0.32;   % gmm
        0.13 0.15 0.17;   % kmeans
        0.18 0.20 0.22;   % region_growing
        0.45 0.50 0.55;   % active_contours
        0.06 0.08 0.10;   % percentile
        0.35 0.40 0.42;   % spectral
        0.08 0.10 0.12;   % fdm
    ];

    ft.fallback_rate = ft.any_failure_rate * 0.6;
    ft.empty_rate = ft.any_failure_rate * 0.2;
    ft.insufficient_rate = ft.any_failure_rate * 0.1;
    ft.all_nan_rate = ft.any_failure_rate * 0.1;

    ft.median_core_voxels = [
        35 32 30;
        30 28 26;
        25 23 21;
        38 35 33;
        40 37 34;
        36 33 30;
        20 18 15;
        22 20 17;
        32 30 28;
        28 25 22;
        30 28 26;
    ];
end
