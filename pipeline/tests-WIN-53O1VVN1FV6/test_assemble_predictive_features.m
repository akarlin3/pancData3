classdef test_assemble_predictive_features < matlab.unittest.TestCase
% TEST_ASSEMBLE_PREDICTIVE_FEATURES — Unit tests for assemble_predictive_features.m
%
% Validates the 22-column feature matrix builder:
%   - Correct column count for mid-treatment vs post-treatment
%   - NaN column removal
%   - Feature name tracking through filtering
%   - Cell array detection error
%   - Index mapping after removal

    properties
        nPat = 10;
        nTp = 5;
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_full_22_columns_mid_treatment(testCase)
            % Mid-treatment fraction should produce 22 columns (before NaN removal)
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                testCase.makeData();
            valid = true(testCase.nPat, 1);
            target_fx = 3;  % mid-treatment

            [X, names, orig_idx, full_names] = assemble_predictive_features( ...
                valid, target_fx, testCase.nTp, 'Fx3', tempdir, ...
                ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds);

            testCase.verifyEqual(numel(full_names), 22, ...
                'Full feature name list should have 22 entries.');
            testCase.verifyEqual(size(X, 1), testCase.nPat, ...
                'Row count should match number of valid patients.');
            testCase.verifyEqual(size(X, 2), numel(names), ...
                'Column count should match feature name count.');
        end

        function test_post_treatment_excludes_dose(testCase)
            % Post-treatment (last timepoint) should exclude dose columns 13-22
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                testCase.makeData();
            valid = true(testCase.nPat, 1);
            target_fx = testCase.nTp;  % last timepoint

            [~, names, orig_idx, full_names] = assemble_predictive_features( ...
                valid, target_fx, testCase.nTp, 'Post', tempdir, ...
                ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds);

            testCase.verifyEqual(numel(full_names), 12, ...
                'Post-treatment should have 12 feature names (no dose).');
            testCase.verifyLessThanOrEqual(max(orig_idx), 12);
        end

        function test_fx6_also_excludes_dose(testCase)
            % Fraction 6 is explicitly excluded even if nTp > 6
            nTpLarge = 8;
            nP = testCase.nPat;
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                testCase.makeDataN(nP, nTpLarge);
            valid = true(nP, 1);

            [~, ~, orig_idx, full_names] = assemble_predictive_features( ...
                valid, 6, nTpLarge, 'Fx6', tempdir, ...
                ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds);

            testCase.verifyEqual(numel(full_names), 12);
        end

        function test_nan_column_removal(testCase)
            % If one parameter is all-NaN, it should be removed
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                testCase.makeData();
            % Make D_abs all-NaN for all patients at Fx1 and target
            D_abs(:,:) = NaN;
            D_pct(:,:) = NaN;
            valid = true(testCase.nPat, 1);

            [X, names, orig_idx, ~] = assemble_predictive_features( ...
                valid, 3, testCase.nTp, 'Fx3', tempdir, ...
                ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds);

            % D_BL (col 2) and D_Abs (col 6) and D_Pct (col 10) should be removed
            testCase.verifyFalse(any(strcmp(names, 'D_BL')));
            testCase.verifyFalse(any(strcmp(names, 'D_Abs')));
            testCase.verifyFalse(any(strcmp(names, 'D_Pct')));
            testCase.verifyEqual(size(X, 2), numel(names));
            testCase.verifyEqual(numel(orig_idx), numel(names));
        end

        function test_valid_pts_mask(testCase)
            % Only valid patients should be included
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                testCase.makeData();
            valid = false(testCase.nPat, 1);
            valid(1:5) = true;

            [X, ~, ~, ~] = assemble_predictive_features( ...
                valid, 3, testCase.nTp, 'Fx3', tempdir, ...
                ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds);

            testCase.verifyEqual(size(X, 1), 5, ...
                'Only 5 valid patients should appear in output.');
        end

        function test_original_indices_map_correctly(testCase)
            % orig_idx should map each output column back to 1:22
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                testCase.makeData();
            valid = true(testCase.nPat, 1);

            [~, ~, orig_idx, ~] = assemble_predictive_features( ...
                valid, 3, testCase.nTp, 'Fx3', tempdir, ...
                ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                d95g, v50g, d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds);

            testCase.verifyGreaterThanOrEqual(min(orig_idx), 1);
            testCase.verifyLessThanOrEqual(max(orig_idx), 22);
            % Indices should be strictly increasing (order preserved)
            testCase.verifyTrue(issorted(orig_idx));
        end
    end

    methods (Access = private)
        function [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, ...
                  f_delta, Dstar_pct, d95g, v50g, d95a, v50a, d95d, v50d, ...
                  d95f, v50f, d95ds, v50ds] = makeData(testCase)
            [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, ...
             f_delta, Dstar_pct, d95g, v50g, d95a, v50a, d95d, v50d, ...
             d95f, v50f, d95ds, v50ds] = testCase.makeDataN(testCase.nPat, testCase.nTp);
        end

        function [ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, ...
                  f_delta, Dstar_pct, d95g, v50g, d95a, v50a, d95d, v50d, ...
                  d95f, v50f, d95ds, v50ds] = makeDataN(~, nP, nT)
            rng(42);
            ADC_abs   = rand(nP, nT) * 0.002;
            D_abs     = rand(nP, nT) * 0.002;
            f_abs     = rand(nP, nT) * 0.3;
            Dstar_abs = rand(nP, nT) * 0.05;
            ADC_pct   = randn(nP, nT) * 20;
            D_pct     = randn(nP, nT) * 20;
            f_delta   = randn(nP, nT) * 0.05;
            Dstar_pct = randn(nP, nT) * 30;
            d95g      = rand(nP, nT) * 50;
            v50g      = rand(nP, nT);
            d95a      = rand(nP, nT) * 50;
            v50a      = rand(nP, nT);
            d95d      = rand(nP, nT) * 50;
            v50d      = rand(nP, nT);
            d95f      = rand(nP, nT) * 50;
            v50f      = rand(nP, nT);
            d95ds     = rand(nP, nT) * 50;
            v50ds     = rand(nP, nT);
        end
    end
end
