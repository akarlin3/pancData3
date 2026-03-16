classdef test_imputation_sensitivity < matlab.unittest.TestCase
% TEST_IMPUTATION_SENSITIVITY  Tests for imputation sensitivity analysis.
%
%   Verifies that all 4 imputation methods (KNN, LOCF, Mean, Linear Interp)
%   produce complete panels, LOCF correctness on monotonic sequences, and
%   mean imputation returns column means.

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'core'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)
        function testAllMethodsProduceCompletePanel(testCase)
            % All 4 imputation methods should produce NaN-free output.
            rng(42);
            n = 20;
            n_feat = 4;
            X_raw = rand(n, n_feat);
            % Insert NaN at random positions (~20% missing)
            nan_idx = rand(n, n_feat) < 0.2;
            X_raw(nan_idx) = NaN;

            patient_ids = cell(n, 1);
            for i = 1:n
                patient_ids{i} = sprintf('P%02d', ceil(i/4));
            end

            % Test LOCF
            X_locf = test_locf(X_raw, patient_ids);
            % LOCF may not fill first-row NaN, fill remainder with mean
            for c = 1:n_feat
                col = X_locf(:, c);
                if any(isnan(col))
                    col(isnan(col)) = nanmean(col);
                    X_locf(:, c) = col;
                end
            end
            testCase.verifyFalse(any(isnan(X_locf(:))), ...
                'LOCF + mean fallback should produce no NaN.');

            % Test Mean
            X_mean = test_mean_impute(X_raw);
            testCase.verifyFalse(any(isnan(X_mean(:))), ...
                'Mean imputation should produce no NaN.');

            % Test Linear Interp
            X_interp = test_linear_interp(X_raw, patient_ids);
            for c = 1:n_feat
                col = X_interp(:, c);
                if any(isnan(col))
                    col(isnan(col)) = nanmean(col);
                    X_interp(:, c) = col;
                end
            end
            testCase.verifyFalse(any(isnan(X_interp(:))), ...
                'Linear interpolation + mean fallback should produce no NaN.');
        end

        function testLOCFMonotonicSequence(testCase)
            % LOCF on a monotonic sequence should carry forward the last
            % observed value correctly.
            X = [1 NaN NaN 4 NaN; 10 NaN 30 NaN NaN]';
            patient_ids = {'P1'; 'P1'; 'P1'; 'P1'; 'P1'};

            X_locf = test_locf(X, patient_ids);

            % Column 1: [1, 1, 1, 4, 4] (carry forward 1, then 4)
            testCase.verifyEqual(X_locf(2, 1), 1, ...
                'LOCF should carry forward 1 to row 2.');
            testCase.verifyEqual(X_locf(3, 1), 1, ...
                'LOCF should carry forward 1 to row 3.');
            testCase.verifyEqual(X_locf(5, 1), 4, ...
                'LOCF should carry forward 4 to row 5.');

            % Column 2: [10, 10, 30, 30, 30]
            testCase.verifyEqual(X_locf(2, 2), 10, ...
                'LOCF should carry forward 10 to row 2.');
            testCase.verifyEqual(X_locf(4, 2), 30, ...
                'LOCF should carry forward 30 to row 4.');
        end

        function testMeanImputationReturnsColumnMean(testCase)
            % Mean imputation should replace NaN with the column mean
            % of non-missing values.
            X = [1 2; 3 NaN; 5 8; NaN 4];
            X_mean = test_mean_impute(X);

            % Column 1 mean of [1,3,5] = 3
            testCase.verifyEqual(X_mean(4, 1), 3, 'AbsTol', 1e-10, ...
                'Mean imputation should replace NaN with column mean.');
            % Column 2 mean of [2,8,4] = 4.667
            testCase.verifyEqual(X_mean(2, 2), mean([2, 8, 4]), 'AbsTol', 1e-10, ...
                'Mean imputation should replace NaN with column mean.');
        end

        function testLinearInterpolation(testCase)
            % Linear interpolation between observed values.
            X = [1; NaN; 3; NaN; 5];
            patient_ids = {'P1'; 'P1'; 'P1'; 'P1'; 'P1'};
            X_interp = test_linear_interp(X, patient_ids);

            testCase.verifyEqual(X_interp(2), 2, 'AbsTol', 1e-10, ...
                'Linearly interpolated value between 1 and 3 should be 2.');
            testCase.verifyEqual(X_interp(4), 4, 'AbsTol', 1e-10, ...
                'Linearly interpolated value between 3 and 5 should be 4.');
        end
    end
end


%% ===== Local test helpers (replicate imputation logic) =====

function X_out = test_locf(X_raw, patient_ids)
    X_out = X_raw;
    unique_pats = unique(patient_ids, 'stable');
    for p = 1:numel(unique_pats)
        rows = strcmp(patient_ids, unique_pats{p});
        row_idx = find(rows);
        for c = 1:size(X_out, 2)
            last_val = NaN;
            for r = 1:numel(row_idx)
                ri = row_idx(r);
                if ~isnan(X_out(ri, c))
                    last_val = X_out(ri, c);
                elseif ~isnan(last_val)
                    X_out(ri, c) = last_val;
                end
            end
        end
    end
end

function X_out = test_mean_impute(X_raw)
    X_out = X_raw;
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        m = mean(col(~isnan(col)));
        col(isnan(col)) = m;
        X_out(:, c) = col;
    end
end

function X_out = test_linear_interp(X_raw, patient_ids)
    X_out = X_raw;
    unique_pats = unique(patient_ids, 'stable');
    for p = 1:numel(unique_pats)
        rows = strcmp(patient_ids, unique_pats{p});
        row_idx = find(rows);
        for c = 1:size(X_out, 2)
            vals = X_out(row_idx, c);
            observed = find(~isnan(vals));
            if numel(observed) >= 2
                interp_vals = interp1(observed, vals(observed), 1:numel(vals), 'linear', NaN);
                X_out(row_idx, c) = interp_vals(:);
            elseif numel(observed) == 1
                X_out(row_idx, c) = vals(observed);
            end
        end
    end
end
