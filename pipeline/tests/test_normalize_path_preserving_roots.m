classdef test_normalize_path_preserving_roots < matlab.unittest.TestCase
    % TEST_NORMALIZE_PATH_PRESERVING_ROOTS  Unit tests for UNC / Unix-absolute
    % path preservation used by compute_spatial_repeatability and
    % optimize_adc_threshold.

    methods (TestClassSetup)
        function addPaths(~)
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, '..', 'utils'));
        end
    end

    methods (Test)
        function unc_prefix_preserved(tc)
            % '\\server\share\...' must keep both leading slashes.
            p_in = '\\pensmph6\mpcsresearch1\aliottae\pancreas_dwi\P1\GTV1.mat';
            out = normalize_path_preserving_roots(p_in);
            tc.verifyTrue(strncmp(out, [filesep filesep], 2), ...
                'UNC prefix was stripped');
        end

        function unc_mixed_separators_preserved(tc)
            % Incoming mixed forward/back slashes, UNC must survive.
            p_in = '\\pensmph6/mpcsresearch1\aliottae/pancreas_dwi/P1/GTV1.mat';
            out = normalize_path_preserving_roots(p_in);
            tc.verifyTrue(strncmp(out, [filesep filesep], 2));
            tc.verifyTrue(contains(out, 'pensmph6'));
            tc.verifyTrue(contains(out, 'GTV1.mat'));
        end

        function unix_absolute_preserved(tc)
            % '/abs/path' keeps its leading slash.
            p_in = '/home/user/pancData3/tests/GTV1.mat';
            out = normalize_path_preserving_roots(p_in);
            tc.verifyEqual(out(1), filesep);
            tc.verifyTrue(contains(out, 'GTV1.mat'));
        end

        function relative_path_unchanged_shape(tc)
            % Relative paths remain relative — no leading filesep added.
            p_in = 'pipeline/utils/file.mat';
            out = normalize_path_preserving_roots(p_in);
            tc.verifyFalse(strncmp(out, filesep, 1));
        end

        function empty_input_returns_empty(tc)
            tc.verifyEqual(normalize_path_preserving_roots(''), '');
        end
    end
end
