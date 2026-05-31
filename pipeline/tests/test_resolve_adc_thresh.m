classdef test_resolve_adc_thresh < matlab.unittest.TestCase
    % TEST_RESOLVE_ADC_THRESH  Sub-volume threshold-source selector.
    %
    % Verifies the config selector that lets the sub-volume be defined by a
    % pre-specified ADC threshold OR by one of the optimizer tactics, with the
    % critical safety property that Tactic 3 NEVER uses its biased in-sample
    % value -- only the nested-CV (out-of-fold) recommendation, falling back
    % to the preset when unavailable.

    methods(TestMethodSetup)
        function setupPaths(~)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testPresetPassthrough(testCase)
            [t, info] = resolve_adc_thresh('preset', 0.0012, struct());
            testCase.verifyEqual(t, 0.0012);
            testCase.verifyEqual(info.source_used, 'preset');
            testCase.verifyFalse(info.fell_back);
        end

        function testDefaultsToPresetWhenEmpty(testCase)
            t = resolve_adc_thresh('', 0.001, struct());
            testCase.verifyEqual(t, 0.001);
        end

        function testProposedFixedCut(testCase)
            % 'proposed' returns the fixed 0.0016 pre-specified cut,
            % independent of preset and opt_results (no selection bias).
            [t, info] = resolve_adc_thresh('proposed', 0.001, struct());
            testCase.verifyEqual(t, 0.0016);
            testCase.verifyEqual(info.source_used, 'proposed');
            testCase.verifyFalse(info.fell_back);
        end

        function testTactic1And2(testCase)
            opt = struct('optimal_thresh', 0.0013, 'inflection_thresh', 0.0017);
            t1 = resolve_adc_thresh('tactic1', 0.001, opt);
            t2 = resolve_adc_thresh('tactic2', 0.001, opt);
            testCase.verifyEqual(t1, 0.0013);
            testCase.verifyEqual(t2, 0.0017);
        end

        function testTactic3UsesNestedCvNotBiased(testCase)
            % CRITICAL: Tactic 3 must use the nested-CV out-of-fold value, NOT
            % the biased in-sample significance_thresh.
            opt = struct();
            opt.significance_thresh = 0.0009;   % the BIASED value -- must be ignored
            opt.nested_cv = struct('recommended_thresh', 0.0015);
            opt.nested_cv_repeated = struct('recommended_thresh', 0.0016);
            [t, info] = resolve_adc_thresh('tactic3', 0.001, opt);
            % Prefers repeated K-fold recommendation.
            testCase.verifyEqual(t, 0.0016);
            testCase.verifyEqual(info.source_used, 'tactic3');
            testCase.verifyNotEqual(t, 0.0009);  % never the biased value
        end

        function testTactic3FallsBackToLooThenPreset(testCase)
            % Repeated unavailable -> use LOO recommendation.
            opt = struct('significance_thresh', 0.0009, ...
                'nested_cv', struct('recommended_thresh', 0.0015), ...
                'nested_cv_repeated', struct('recommended_thresh', NaN));
            t = resolve_adc_thresh('tactic3', 0.001, opt);
            testCase.verifyEqual(t, 0.0015);

            % Both nested-CV unavailable -> fall back to preset, never biased.
            opt2 = struct('significance_thresh', 0.0009, ...
                'nested_cv', struct('recommended_thresh', NaN), ...
                'nested_cv_repeated', struct('recommended_thresh', NaN));
            [t2, info2] = resolve_adc_thresh('tactic3', 0.001, opt2);
            testCase.verifyEqual(t2, 0.001);
            testCase.verifyTrue(info2.fell_back);
            testCase.verifyNotEqual(t2, 0.0009);
        end

        function testUnavailableTacticFallsBack(testCase)
            [t, info] = resolve_adc_thresh('tactic1', 0.0011, struct());
            testCase.verifyEqual(t, 0.0011);
            testCase.verifyTrue(info.fell_back);
        end

        function testUnknownSourceFallsBack(testCase)
            [t, info] = resolve_adc_thresh('bogus', 0.001, struct());
            testCase.verifyEqual(t, 0.001);
            testCase.verifyTrue(info.fell_back);
        end

    end
end
