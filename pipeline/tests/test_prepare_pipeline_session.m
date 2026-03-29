classdef test_prepare_pipeline_session < matlab.unittest.TestCase
    % TEST_PREPARE_PIPELINE_SESSION Unit tests for prepare_pipeline_session.
    %
    % Validates session setup: config parsing, DWI type resolution, output
    % folder creation, file path construction, config-gated step injection
    % (compare_cores, cross_pipeline_dice, core_failure_rates,
    % core_method_outcomes), and error handling for invalid config.

    properties
        TmpDir
        PipelineDir
        ConfigFile
    end

    methods(TestMethodSetup)
        function setupTempDir(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);

            % Suppress pipeline progress GUIs during testing
            setenv('SUPPRESS_PIPELINE_GUI', '1');

            % Create a minimal pipeline directory structure
            testCase.PipelineDir = fullfile(testCase.TmpDir, 'pipeline');
            mkdir(testCase.PipelineDir);
            mkdir(fullfile(testCase.PipelineDir, 'core'));
            mkdir(fullfile(testCase.PipelineDir, 'utils'));

            % Add utils to path so prepare_pipeline_session can call helpers
            current_dir = fileparts(mfilename('fullpath'));
            addpath(fullfile(current_dir, '..', 'utils'));
            addpath(fullfile(current_dir, '..', 'core'));

            % Create a minimal valid config
            testCase.ConfigFile = fullfile(testCase.TmpDir, 'config.json');
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = false;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);
        end
    end

    methods(TestMethodTeardown)
        function cleanupTmpDir(testCase)
            diary off;
            % Restore figure visibility in case a test left it off
            set(0, 'DefaultFigureVisible', 'on');
            if isfolder(testCase.TmpDir)
                rmdir(testCase.TmpDir, 's');
            end
        end
    end

    methods(Test)
        function test_returns_session_struct(testCase)
            % Session should contain all required fields.
            steps = {'load', 'sanity'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            % Restore figure visibility
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyFalse(session.abort);
            testCase.verifyTrue(isstruct(session.config_struct));
            testCase.verifyEqual(session.current_name, 'Standard');
            testCase.verifyEqual(session.current_dtype, 1);
            testCase.verifyTrue(session.log_fid > 0);
            % Clean up log handle
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_file_paths_contain_dwi_type(testCase)
            % All artifact paths should contain the DWI type name.
            steps = {'load'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifySubstring(session.dwi_vectors_file, 'Standard');
            testCase.verifySubstring(session.summary_metrics_file, 'Standard');
            testCase.verifySubstring(session.results_file, 'Standard');
            testCase.verifySubstring(session.baseline_results_file, 'Standard');
            testCase.verifySubstring(session.dosimetry_results_file, 'Standard');
            testCase.verifySubstring(session.predictive_results_file, 'Standard');
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_compare_cores_injected_when_config_true(testCase)
            % When run_compare_cores is true, compare_cores should be
            % injected after metrics_baseline in steps_to_run.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = true;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            steps = {'load', 'metrics_baseline', 'metrics_longitudinal'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyTrue(ismember('compare_cores', session.steps_to_run));
            % compare_cores should appear after metrics_baseline
            idx_base = find(strcmp(session.steps_to_run, 'metrics_baseline'));
            idx_comp = find(strcmp(session.steps_to_run, 'compare_cores'));
            testCase.verifyGreaterThan(idx_comp, idx_base);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_compare_cores_not_duplicated(testCase)
            % If compare_cores is already in steps, it should not be added again.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = true;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            steps = {'load', 'metrics_baseline', 'compare_cores'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            n_compare = sum(strcmp(session.steps_to_run, 'compare_cores'));
            testCase.verifyEqual(n_compare, 1);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_cross_pipeline_dice_injected(testCase)
            % When run_cross_pipeline_dice is true, cross_pipeline_dice
            % should be injected after metrics_baseline.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = false;
            cfg.run_cross_pipeline_dice = true;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            steps = {'load', 'metrics_baseline', 'metrics_longitudinal'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyTrue(ismember('cross_pipeline_dice', session.steps_to_run));
            idx_base = find(strcmp(session.steps_to_run, 'metrics_baseline'));
            idx_cpd = find(strcmp(session.steps_to_run, 'cross_pipeline_dice'));
            testCase.verifyGreaterThan(idx_cpd, idx_base);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_core_failure_rates_injected(testCase)
            % When run_core_failure_rates is true, core_failure_rates
            % should be injected after metrics_baseline.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = false;
            cfg.run_core_failure_rates = true;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            steps = {'load', 'metrics_baseline', 'metrics_longitudinal'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyTrue(ismember('core_failure_rates', session.steps_to_run));
            idx_base = find(strcmp(session.steps_to_run, 'metrics_baseline'));
            idx_cfr = find(strcmp(session.steps_to_run, 'core_failure_rates'));
            testCase.verifyGreaterThan(idx_cfr, idx_base);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_core_method_outcomes_injected(testCase)
            % When run_core_method_outcomes is true, core_method_outcomes
            % should be injected after metrics_dosimetry.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = false;
            cfg.run_core_method_outcomes = true;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            steps = {'load', 'metrics_baseline', 'metrics_dosimetry', 'metrics_survival'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyTrue(ismember('core_method_outcomes', session.steps_to_run));
            idx_dos = find(strcmp(session.steps_to_run, 'metrics_dosimetry'));
            idx_cmo = find(strcmp(session.steps_to_run, 'core_method_outcomes'));
            testCase.verifyGreaterThan(idx_cmo, idx_dos);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_all_config_gated_steps_injected(testCase)
            % When all config-gated steps are enabled, they should all
            % appear in the correct order.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = true;
            cfg.run_cross_pipeline_dice = true;
            cfg.run_core_failure_rates = true;
            cfg.run_core_method_outcomes = true;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            steps = {'load', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_survival'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            s = session.steps_to_run;
            testCase.verifyTrue(ismember('compare_cores', s));
            testCase.verifyTrue(ismember('cross_pipeline_dice', s));
            testCase.verifyTrue(ismember('core_failure_rates', s));
            testCase.verifyTrue(ismember('core_method_outcomes', s));

            % Verify ordering
            idx_base = find(strcmp(s, 'metrics_baseline'));
            idx_cc   = find(strcmp(s, 'compare_cores'));
            idx_cpd  = find(strcmp(s, 'cross_pipeline_dice'));
            idx_cfr  = find(strcmp(s, 'core_failure_rates'));
            idx_dos  = find(strcmp(s, 'metrics_dosimetry'));
            idx_cmo  = find(strcmp(s, 'core_method_outcomes'));

            testCase.verifyGreaterThan(idx_cc, idx_base);
            testCase.verifyGreaterThan(idx_cpd, idx_cc);
            testCase.verifyGreaterThan(idx_cfr, idx_cpd);
            testCase.verifyGreaterThan(idx_cmo, idx_dos);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_invalid_config_sets_abort(testCase)
            % Bad config file should set abort=true without throwing.
            % The function should also restore DefaultFigureVisible on
            % abort so the caller does not need to clean up.
            bad_config = fullfile(testCase.TmpDir, 'bad_config.json');
            fid = fopen(bad_config, 'w');
            fprintf(fid, '{{invalid json');
            fclose(fid);

            set(0, 'DefaultFigureVisible', 'on');
            session = prepare_pipeline_session(testCase.PipelineDir, bad_config, '', {'load'});

            testCase.verifyTrue(session.abort);
            % Figure visibility should have been restored by the catch block
            testCase.verifyEqual(char(get(0, 'DefaultFigureVisible')), 'on');
        end

        function test_type_output_folder_created(testCase)
            % The DWI-type subfolder should be created.
            steps = {'load'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyTrue(isfolder(session.type_output_folder));
            [~, type_name] = fileparts(session.type_output_folder);
            testCase.verifyEqual(type_name, 'Standard');
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_explicit_output_folder_reused(testCase)
            % When master_output_folder is provided, it should be reused.
            explicit = fullfile(testCase.TmpDir, 'my_output');
            mkdir(explicit);
            steps = {'load'};
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, explicit, steps);
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyEqual(session.master_output_folder, explicit);
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_dncnn_type_resolution(testCase)
            % dwi_type 'dnCNN' should resolve to dtype=2, name='dnCNN'.
            cfg = struct();
            cfg.dataloc = testCase.TmpDir;
            cfg.dwi_type = 'dnCNN';
            cfg.run_compare_cores = false;
            cfg.run_all_core_methods = false;
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', {'load'});
            set(0, 'DefaultFigureVisible', session.prev_fig_vis);
            diary off;

            testCase.verifyEqual(session.current_dtype, 2);
            testCase.verifyEqual(session.current_name, 'dnCNN');
            if session.log_fid > 0, fclose(session.log_fid); end
        end

        function test_post_config_error_restores_state(testCase)
            % Verify that prepare_pipeline_session correctly handles state
            % cleanup. When the function succeeds, the caller owns cleanup
            % (via session.prev_fig_vis). When it throws, the catch block
            % inside the function restores DefaultFigureVisible.
            cfg = struct();
            cfg.dataloc = fullfile(testCase.TmpDir, 'nonexistent_data_dir');
            cfg.dwi_type = 'Standard';
            cfg.run_compare_cores = false;
            cfg.run_all_core_methods = false;
            cfg.clear_cache = true;
            bad_cfg_file = fullfile(testCase.TmpDir, 'post_config_fail.json');
            fid = fopen(bad_cfg_file, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            set(0, 'DefaultFigureVisible', 'on');
            threw = false;
            session = struct('prev_fig_vis', 'on', 'log_fid', -1);
            try
                session = prepare_pipeline_session(testCase.PipelineDir, bad_cfg_file, '', {'load'});
            catch
                threw = true;
            end
            diary off;

            if threw
                % The catch-and-rethrow path restores visibility
                testCase.verifyEqual(char(get(0, 'DefaultFigureVisible')), 'on');
            else
                % Function succeeded — caller owns cleanup; verify prev_fig_vis
                testCase.verifyEqual(char(session.prev_fig_vis), 'on', ...
                    'prev_fig_vis should preserve the original visibility state.');
                set(0, 'DefaultFigureVisible', session.prev_fig_vis);
                if session.log_fid > 0, fclose(session.log_fid); end
            end
        end

        function test_figure_visibility_set_to_off(testCase)
            % prepare_pipeline_session should set figure visibility to 'off'.
            set(0, 'DefaultFigureVisible', 'on');
            session = prepare_pipeline_session(testCase.PipelineDir, testCase.ConfigFile, '', {'load'});
            diary off;

            testCase.verifyEqual(char(get(0, 'DefaultFigureVisible')), 'off');
            testCase.verifyEqual(char(session.prev_fig_vis), 'on');
            set(0, 'DefaultFigureVisible', 'on');
            if session.log_fid > 0, fclose(session.log_fid); end
        end
    end
end
