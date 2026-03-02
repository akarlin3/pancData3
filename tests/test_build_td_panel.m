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

        function testIntervalBoundariesMatchScanSchedule(testCase)
            % Verify that interval start/stop times exactly match the
            % scan schedule and survival time.
            nTp = 3;
            scan_days = [0, 10, 20];
            feat_arrays = {[5, 6, 7]};  % 1 patient, 3 timepoints
            feat_names = {'F1'};
            lf_vec = 1;            % event at survival time
            total_time_vec = 50;   % survives to day 50

            [~, t_start, t_stop, event_td, ~, frac_td] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            % Expected: 3 intervals
            %   [0, 10), [10, 20), [20, 50]
            testCase.verifyEqual(numel(t_start), 3, 'Should have 3 intervals.');
            testCase.verifyEqual(t_start, [0; 10; 20], 'AbsTol', 1e-10);
            testCase.verifyEqual(t_stop, [10; 20; 50], 'AbsTol', 1e-10);

            % Event should be assigned to the last interval only
            testCase.verifyEqual(event_td, [0; 0; 1]);

            % Fraction indices
            testCase.verifyEqual(frac_td, [1; 2; 3]);
        end

        function testEventAssignedToFinalIntervalOnly(testCase)
            % For a patient with an event, only the last interval should
            % carry the event indicator.
            nTp = 4;
            scan_days = [0, 5, 10, 15];
            feat_arrays = {[1, 2, 3, 4]};
            feat_names = {'F1'};
            lf_vec = 1;
            total_time_vec = 12;  % dies at day 12, between scan 3 and 4

            [~, t_start, t_stop, event_td, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            % Intervals: [0, 5), [5, 10), [10, 12]
            % Scan 4 at day 15 >= 12 so not included
            testCase.verifyEqual(numel(event_td), 3);
            testCase.verifyEqual(event_td(end), 1, ...
                'Event should be at the final interval.');
            testCase.verifyEqual(sum(event_td(1:end-1)), 0, ...
                'Non-final intervals should have event=0.');
        end

        function testCompetingRiskCodePropagated(testCase)
            % Event code 2 (competing risk) should be passed through to
            % the terminal interval.
            nTp = 2;
            scan_days = [0, 10];
            feat_arrays = {[1, 2]};
            feat_names = {'F1'};
            lf_vec = 2;            % competing risk event
            total_time_vec = 30;

            [~, ~, ~, event_td, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            % The terminal interval should carry event code 2
            testCase.verifyEqual(event_td(end), 2, ...
                'Competing risk code should be propagated to terminal interval.');
        end

        function testCensoredPatientHasNoEvents(testCase)
            % A censored patient (lf=0) should have event=0 on all intervals.
            nTp = 3;
            scan_days = [0, 10, 20];
            feat_arrays = {[1, 2, 3]};
            feat_names = {'F1'};
            lf_vec = 0;
            total_time_vec = 50;

            [~, ~, ~, event_td, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            testCase.verifyTrue(all(event_td == 0), ...
                'Censored patient should have event=0 on all intervals.');
        end

        function testMultiplePatientsCorrectPatientIDs(testCase)
            % Each patient's intervals should have the correct patient ID.
            nTp = 2;
            scan_days = [0, 10];
            feat_arrays = {[1, 2; 3, 4; 5, 6]};  % 3 patients
            feat_names = {'F1'};
            lf_vec = [1; 0; 1];
            total_time_vec = [25; 30; 15];

            [~, ~, ~, event_td, pat_id_td, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            % Each patient should appear in pat_id_td
            testCase.verifyTrue(ismember(1, pat_id_td));
            testCase.verifyTrue(ismember(2, pat_id_td));
            testCase.verifyTrue(ismember(3, pat_id_td));

            % Events only for patients 1 and 3
            for pid = [1, 3]
                pid_events = event_td(pat_id_td == pid);
                testCase.verifyEqual(pid_events(end), lf_vec(pid), ...
                    sprintf('Patient %d terminal event should match lf_vec.', pid));
            end

            % Patient 2 should have no events
            p2_events = event_td(pat_id_td == 2);
            testCase.verifyTrue(all(p2_events == 0), ...
                'Censored patient 2 should have no events.');
        end

        function testScanDayAfterSurvivalTimeExcluded(testCase)
            % If a scan day >= survival time, that interval should not be generated.
            nTp = 4;
            scan_days = [0, 5, 10, 100];
            feat_arrays = {[1, 2, 3, 4]};
            feat_names = {'F1'};
            lf_vec = 1;
            total_time_vec = 8;  % Dies at day 8, before scan 3

            [~, t_start, t_stop, ~, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            % Only 2 intervals: [0, 5), [5, 8]
            testCase.verifyEqual(numel(t_start), 2);
            testCase.verifyEqual(t_stop(end), 8);
        end

        function testMultipleFeaturesCorrectColumns(testCase)
            % When multiple features are passed, X_td columns should match
            nTp = 2;
            scan_days = [0, 10];
            feat_arrays = {[1, 2], [10, 20], [100, 200]};  % 3 features
            feat_names = {'F1', 'F2', 'F3'};
            lf_vec = 0;
            total_time_vec = 50;

            [X_td, ~, ~, ~, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            testCase.verifyEqual(size(X_td, 2), 3, ...
                'X_td should have 3 feature columns.');
            % First interval uses tp=1 values
            testCase.verifyEqual(X_td(1, :), [1, 10, 100], 'AbsTol', 1e-10);
            % Second interval uses tp=2 values
            testCase.verifyEqual(X_td(2, :), [2, 20, 200], 'AbsTol', 1e-10);
        end
    end
end
