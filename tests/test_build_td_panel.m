classdef test_build_td_panel < matlab.unittest.TestCase
    % TEST_BUILD_TD_PANEL Unit tests for the counting-process panel builder
    methods(TestMethodSetup)
        function addPath(testCase)
            % Ensure utils is on the path
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testEmptyInputs(testCase)
            % Test with completely empty inputs
            feat_arrays = {};
            feat_names = {};
            lf_vec = [];
            total_time_vec = [];
            nTp = 5;
            scan_days = [0 5 10 15 20];
            decay_half_life = 18;

            % Verify that it warns 'build_td_panel:noData'
            testCase.verifyWarning(@() build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, decay_half_life), ...
                'build_td_panel:noData');

            % Verify outputs are empty
            warning('off', 'build_td_panel:noData');
            [X_td, t_start, t_stop, event_td, pat_id_td, frac_td] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, decay_half_life);
            warning('on', 'build_td_panel:noData');

            testCase.verifyEmpty(X_td);
            testCase.verifyEmpty(t_start);
            testCase.verifyEmpty(t_stop);
            testCase.verifyEmpty(event_td);
            testCase.verifyEmpty(pat_id_td);
            testCase.verifyEmpty(frac_td);
        end

        function testAllNaNInputs(testCase)
            % Test with NaN inputs which should result in no valid intervals
            n_pts = 3;
            nTp = 3;
            feat_arrays = {nan(n_pts, nTp)};
            feat_names = {'Feature1'};
            lf_vec = [0; 0; 1];
            % If total_time_vec is NaN, patient is invalid
            total_time_vec = [NaN; NaN; NaN];
            scan_days = [0 5 10];
            decay_half_life = 18;

            testCase.verifyWarning(@() build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, decay_half_life), ...
                'build_td_panel:noData');

             warning('off', 'build_td_panel:noData');
            [X_td, ~, ~, ~, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, decay_half_life);
            warning('on', 'build_td_panel:noData');

            testCase.verifyEmpty(X_td);
        end

        function testPartialNaNInputs(testCase)
            % Test where some patients are valid and some are NaN
            n_pts = 2;
            nTp = 2;
            % Patient 1: valid data
            % Patient 2: NaN time
            feat_arrays = {[1 2; 3 4]}; % 2 pts, 2 tp
            feat_names = {'F1'};
            lf_vec = [0; 1];
            total_time_vec = [100; NaN]; % Patient 2 has NaN time -> invalid
            scan_days = [0 10];

            % Should NOT warn because at least one patient (Pt 1) is valid
            [X_td, ~, ~, ~, pat_id_td, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            testCase.verifyNotEmpty(X_td);
            testCase.verifyEqual(size(X_td, 1), 2); % 2 intervals for patient 1
            testCase.verifyEqual(unique(pat_id_td), 1); % Only patient 1
        end
    end
end
