classdef test_make_grouped_folds < matlab.unittest.TestCase
% TEST_MAKE_GROUPED_FOLDS — Unit tests for make_grouped_folds.m
%
% Validates stratified grouped k-fold cross-validation:
%   - Patient grouping (all rows for same patient in same fold)
%   - Correct number of folds
%   - Stratification of event rate across folds
%   - Edge cases (single patient, more folds than patients)
%   - Competing risk handling
%   - Output dimensions

    properties
        origPath
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.origPath = path;
            baseDir = fullfile(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods (TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.origPath);
        end
    end

    methods (Test)
        function testAllRowsSamePatientSameFold(testCase)
            % 3 patients with multiple rows each; all rows for the same
            % patient must land in the same fold.
            rng(42);
            ids = {'P1'; 'P1'; 'P1'; 'P2'; 'P2'; 'P3'; 'P3'; 'P3'; 'P3'};
            y   = [1; 0; 0; 0; 0; 1; 0; 0; 0];
            evalc('fold_id = make_grouped_folds(ids, y, 3);');

            % All P1 rows should share a fold
            testCase.verifyEqual(numel(unique(fold_id([1 2 3]))), 1, ...
                'All rows for P1 must be in the same fold.');
            % All P2 rows should share a fold
            testCase.verifyEqual(numel(unique(fold_id([4 5]))), 1, ...
                'All rows for P2 must be in the same fold.');
            % All P3 rows should share a fold
            testCase.verifyEqual(numel(unique(fold_id([6 7 8 9]))), 1, ...
                'All rows for P3 must be in the same fold.');
        end

        function testCorrectNumberOfFolds(testCase)
            rng(42);
            n_patients = 10;
            ids = cell(n_patients * 2, 1);
            y   = zeros(n_patients * 2, 1);
            for i = 1:n_patients
                ids{2*i-1} = sprintf('P%d', i);
                ids{2*i}   = sprintf('P%d', i);
                if i <= 5, y(2*i-1) = 1; end
            end
            n_folds = 5;
            evalc('fold_id = make_grouped_folds(ids, y, n_folds);');
            testCase.verifyEqual(sort(unique(fold_id))', 1:n_folds, ...
                'unique(fold_id) should contain exactly 1:n_folds.');
        end

        function testStratificationPreservesEventRate(testCase)
            % With ~50% event rate across 20 patients, each of 5 folds
            % should have at least 1 event and 1 non-event patient.
            rng(42);
            n_patients = 20;
            ids = cell(n_patients, 1);
            y   = zeros(n_patients, 1);
            for i = 1:n_patients
                ids{i} = sprintf('P%d', i);
                if i <= 10, y(i) = 1; end
            end
            n_folds = 5;
            evalc('fold_id = make_grouped_folds(ids, y, n_folds);');

            for f = 1:n_folds
                mask = (fold_id == f);
                n_events = sum(y(mask) == 1);
                n_non    = sum(y(mask) == 0);
                testCase.verifyGreaterThan(n_events, 0, ...
                    sprintf('Fold %d should have at least 1 event.', f));
                testCase.verifyGreaterThan(n_non, 0, ...
                    sprintf('Fold %d should have at least 1 non-event.', f));
            end
        end

        function testSinglePatientReturnsSingleFold(testCase)
            ids = {'P1'; 'P1'; 'P1'};
            y   = [1; 0; 0];
            evalc('fold_id = make_grouped_folds(ids, y, 5);');
            testCase.verifyTrue(all(fold_id == 1), ...
                'Single patient should produce all fold_id = 1.');
        end

        function testMoreFoldsThanPatients(testCase)
            % 3 unique patients but requesting 10 folds; should clamp to 3.
            rng(42);
            ids = {'A'; 'A'; 'B'; 'B'; 'C'; 'C'};
            y   = [1; 0; 0; 0; 1; 0];
            evalc('fold_id = make_grouped_folds(ids, y, 10);');
            n_actual_folds = numel(unique(fold_id));
            testCase.verifyLessThanOrEqual(n_actual_folds, 3, ...
                'Number of folds should be clamped to number of unique patients.');
            testCase.verifyGreaterThanOrEqual(n_actual_folds, 1, ...
                'There should be at least 1 fold.');
        end

        function testCompetingRisksNotCountedAsEvents(testCase)
            % y=2 (competing risk) should be treated like y=0 for
            % stratification; only y=1 counts as an event.
            rng(42);
            n_patients = 20;
            ids = cell(n_patients, 1);
            y   = zeros(n_patients, 1);
            for i = 1:n_patients
                ids{i} = sprintf('P%d', i);
            end
            % 5 events, 5 competing risks, 10 censored
            y(1:5) = 1;
            y(6:10) = 2;

            evalc('fold_id = make_grouped_folds(ids, y, 5);');

            % Verify fold assignment is valid
            testCase.verifyTrue(all(fold_id >= 1 & fold_id <= 5), ...
                'All fold IDs should be in range 1:5.');
            % Competing risk patients (y=2) should be distributed similarly
            % to censored patients (y=0), not lumped with events.
            for f = 1:5
                mask = (fold_id == f);
                % Each fold should have at least 1 event among its patients
                testCase.verifyGreaterThan(sum(y(mask) == 1), 0, ...
                    sprintf('Fold %d should have at least 1 local failure event.', f));
            end
        end

        function testFoldIdDimensions(testCase)
            rng(42);
            ids = {'P1'; 'P1'; 'P2'; 'P2'; 'P3'};
            y   = [1; 0; 0; 0; 1];
            evalc('fold_id = make_grouped_folds(ids, y, 3);');
            testCase.verifyEqual(size(fold_id), [5, 1], ...
                'fold_id should be an n_rows x 1 column vector.');
        end
    end
end
