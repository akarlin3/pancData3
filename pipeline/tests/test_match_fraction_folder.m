classdef test_match_fraction_folder < matlab.unittest.TestCase
    % TEST_MATCH_FRACTION_FOLDER  Unit tests for match_fraction_folder.m.
    %
    % match_fraction_folder replaces the earlier regex-based fraction-folder
    % matcher used in discover_patient_files and load_dwi_data. These tests
    % lock in the behaviour that was previously fragile under MATLAB regex:
    %   * 'Fx1 - repeatability' is recognised as an Fx1 folder
    %   * 'Fx10' / 'Fx11' are NOT recognised as Fx1
    %   * 'postprocessing' is NOT recognised as 'post'
    %   * case-insensitive matching still works

    methods(TestMethodSetup)
        function addPathsForMatcher(testCase) %#ok<INUSD>
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(Test)
        function testExactMatch(testCase)
            mask = match_fraction_folder({'Fx1'}, 'Fx1');
            testCase.verifyEqual(mask, true);
        end

        function testFx1DashRepeatabilityMatches(testCase)
            % Regression: MSK cohort uses 'Fx1 - repeatability' for the
            % baseline repeat-scan session. This is the exact case that
            % the earlier regex-based matcher silently dropped.
            mask = match_fraction_folder({'Fx1 - repeatability'}, 'Fx1');
            testCase.verifyEqual(mask, true, ...
                '''Fx1 - repeatability'' must match Fx1.');
        end

        function testUnderscoreSuffixMatches(testCase)
            mask = match_fraction_folder({'Fx1_old'}, 'Fx1');
            testCase.verifyEqual(mask, true);
        end

        function testFx10DoesNotMatchFx1(testCase)
            % Fx10 must not be confused with Fx1 — the original intent
            % of the word-boundary guard. This is the regression the
            % previous '\b' regex was trying (and succeeding) to prevent.
            mask = match_fraction_folder({'Fx10'}, 'Fx1');
            testCase.verifyEqual(mask, false);
        end

        function testFx11DoesNotMatchFx1(testCase)
            mask = match_fraction_folder({'Fx11'}, 'Fx1');
            testCase.verifyEqual(mask, false);
        end

        function testCaseInsensitive(testCase)
            mask = match_fraction_folder({'fx1', 'FX1', 'Fx1'}, 'Fx1');
            testCase.verifyEqual(mask, [true, true, true]);
        end

        function testPostprocessingDoesNotMatchPost(testCase)
            % 'postprocessing' starts with 'post' but continues with a
            % letter, so it must NOT be treated as the post-RT scan.
            mask = match_fraction_folder({'postprocessing'}, 'post');
            testCase.verifyEqual(mask, false);
        end

        function testPostMatchesPost(testCase)
            mask = match_fraction_folder({'post'}, 'post');
            testCase.verifyEqual(mask, true);
        end

        function testPostDashVariantMatches(testCase)
            mask = match_fraction_folder({'post - treatment'}, 'post');
            testCase.verifyEqual(mask, true);
        end

        function testMixedBatchReturnsPerEntryMask(testCase)
            names = {'Fx1', 'Fx1 - repeatability', 'Fx10', 'post', 'Fx2'};
            mask = match_fraction_folder(names, 'Fx1');
            testCase.verifyEqual(mask, [true, true, false, false, false]);
        end

        function testEmptyInputReturnsEmpty(testCase)
            mask = match_fraction_folder({}, 'Fx1');
            testCase.verifyEqual(mask, false(1, 0));
        end

        function testNonMatchingPrefixReturnsFalse(testCase)
            mask = match_fraction_folder({'xFx1', 'DWI1', 'rtdose'}, 'Fx1');
            testCase.verifyEqual(mask, [false, false, false]);
        end
    end
end
