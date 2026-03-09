classdef test_build_td_panel < matlab.unittest.TestCase
    % TEST_BUILD_TD_PANEL Unit tests for the counting-process panel builder.
    %
    % build_td_panel converts per-patient, per-timepoint feature arrays and
    % survival information into a counting-process (start, stop] panel
    % suitable for time-dependent Cox regression. Each patient contributes
    % one row per scan interval, with the event indicator placed only on the
    % terminal interval.
    %
    % These tests cover:
    %   - Empty and all-NaN inputs (noData warning)
    %   - Partial NaN (some patients invalid, others valid)
    %   - Interval boundary alignment with the scan schedule
    %   - Event assignment to the final interval only
    %   - Competing risk code propagation (lf=2)
    %   - Censored patients (lf=0, no events anywhere)
    %   - Multi-patient panel with correct patient IDs
    %   - Scan days beyond survival time (interval exclusion)
    %   - Multiple feature columns in X_td

    methods(TestMethodSetup)
        function addPath(testCase)
            % Ensure utils is on the path
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testEmptyInputs(testCase)
            % Completely empty inputs (no patients, no features) should
            % trigger a 'build_td_panel:noData' warning and return empty
            % arrays for all six outputs.
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
            % When all patients have NaN survival times, none are valid for
            % panel construction. Should warn and return empty X_td.
            n_pts = 3;
            nTp = 3;
            feat_arrays = {nan(n_pts, nTp)};  % feature values are NaN too
            feat_names = {'Feature1'};
            lf_vec = [0; 0; 1];
            % NaN total_time_vec makes all patients invalid
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
            % Mixed validity: Patient 1 has valid data, Patient 2 has NaN
            % survival time. Only Patient 1 should appear in the panel.
            n_pts = 2;
            nTp = 2;
            feat_arrays = {[1 2; 3 4]}; % 2 patients, 2 timepoints
            feat_names = {'F1'};
            lf_vec = [0; 1];
            total_time_vec = [100; NaN]; % Patient 2 has NaN time -> invalid
            scan_days = [0 10];

            % Should NOT warn because at least one patient (Pt 1) is valid
            [X_td, ~, ~, ~, pat_id_td, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            testCase.verifyNotEmpty(X_td);
            % Patient 1 with 2 timepoints and survival time 100 generates 2 intervals
            testCase.verifyEqual(size(X_td, 1), 2);
            % Only patient 1 should be present in the panel
            testCase.verifyEqual(unique(pat_id_td), 1);
        end

        function testIntervalBoundariesMatchScanSchedule(testCase)
            % Verify that interval start/stop times exactly match the
            % scan schedule and survival time for a single patient with
            % 3 scans at days 0, 10, 20 and survival to day 50.
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

            % Fraction indices should correspond to timepoint order
            testCase.verifyEqual(frac_td, [1; 2; 3]);
        end

        function testEventAssignedToFinalIntervalOnly(testCase)
            % Patient dies at day 12, between scan 3 (day 10) and scan 4
            % (day 15). Scan 4 is after death so its interval is excluded.
            % Only the last valid interval [10, 12] carries the event.
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
            % Event code 2 (competing risk, e.g., non-cancer death) should
            % be passed through unchanged to the terminal interval, enabling
            % downstream competing risks regression.
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
            % A censored patient (lf=0, no event observed) should have
            % event=0 on every interval in the counting-process panel.
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
            % Three patients with different event statuses: Pt1 (event),
            % Pt2 (censored), Pt3 (event). Verify that pat_id_td correctly
            % labels each interval's owner and that terminal events match.
            nTp = 2;
            scan_days = [0, 10];
            feat_arrays = {[1, 2; 3, 4; 5, 6]};  % 3 patients, 2 timepoints
            feat_names = {'F1'};
            lf_vec = [1; 0; 1];
            total_time_vec = [25; 30; 15];

            [~, ~, ~, event_td, pat_id_td, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            % Each patient should appear in pat_id_td
            testCase.verifyTrue(ismember(1, pat_id_td));
            testCase.verifyTrue(ismember(2, pat_id_td));
            testCase.verifyTrue(ismember(3, pat_id_td));

            % Terminal events for patients 1 and 3 should match their lf_vec
            for pid = [1, 3]
                pid_events = event_td(pat_id_td == pid);
                testCase.verifyEqual(pid_events(end), lf_vec(pid), ...
                    sprintf('Patient %d terminal event should match lf_vec.', pid));
            end

            % Patient 2 (censored) should have no events on any interval
            p2_events = event_td(pat_id_td == 2);
            testCase.verifyTrue(all(p2_events == 0), ...
                'Censored patient 2 should have no events.');
        end

        function testScanDayAfterSurvivalTimeExcluded(testCase)
            % Patient dies at day 8. Scan 3 (day 10) and scan 4 (day 100)
            % occur after death, so their intervals should not be generated.
            % Only intervals [0, 5) and [5, 8] should exist.
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
            % Three features passed as separate arrays should produce
            % X_td with 3 columns, each carrying the correct timepoint value.
            nTp = 2;
            scan_days = [0, 10];
            feat_arrays = {[1, 2], [10, 20], [100, 200]};  % 3 features, 1 patient
            feat_names = {'F1', 'F2', 'F3'};
            lf_vec = 0;
            total_time_vec = 50;

            [X_td, ~, ~, ~, ~, ~] = build_td_panel( ...
                feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, 18);

            testCase.verifyEqual(size(X_td, 2), 3, ...
                'X_td should have 3 feature columns.');
            % First interval uses timepoint 1 values: [1, 10, 100]
            testCase.verifyEqual(X_td(1, :), [1, 10, 100], 'AbsTol', 1e-10);
            % Second interval uses timepoint 2 values: [2, 20, 200]
            testCase.verifyEqual(X_td(2, :), [2, 20, 200], 'AbsTol', 1e-10);
        end
    end
end
