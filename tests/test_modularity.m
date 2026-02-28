classdef test_modularity < matlab.unittest.TestCase
    % TEST_MODULARITY Unit tests for skipped pipeline steps and modularity
    properties
        TempDir
        ConfigPath
        ConfigStruct
        RepoRoot
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.dataloc = [testCase.TempDir filesep];
            testCase.ConfigStruct.ivim_bthr = 100;
            testCase.ConfigStruct.skip_to_reload = true;
            testCase.ConfigStruct.use_checkpoints = true;
            % These are needed by load_dwi_data
            testCase.ConfigStruct.adc_thresh = 0.00115;
            testCase.ConfigStruct.high_adc_thresh = 0.001;
            testCase.ConfigStruct.min_vox_hist = 1;
            testCase.ConfigStruct.adc_max = 3e-3;
            testCase.ConfigStruct.dcm2nii_call = 'dcm2niix';
            testCase.ConfigStruct.clinical_data_sheet = 'clin.xlsx';

            testCase.ConfigPath = fullfile(testCase.TempDir, 'config.json');

            % Write config file
            fid = fopen(testCase.ConfigPath, 'w');
            fprintf(fid, '%s', jsonencode(testCase.ConfigStruct));
            fclose(fid);

            % Add repo root to path so run_dwi_pipeline is found
            % `pwd` might be inside `tests/`, so use fileparts
            [testCase.RepoRoot, ~, ~] = fileparts(pwd);
            addpath(testCase.RepoRoot);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            rmdir(testCase.TempDir, 's');
        end
    end

    methods(Test)
        function testSkipLoad(testCase)
            % Create dummy data files to simulate existing load
            data_vectors_gtvp = struct('adc_vector', [], 'd_vector', []);
            data_vectors_gtvn = struct();
            summary_metrics = struct('id_list', {{'Test'}}, 'mrn_list', {{'123'}}, ...
                'adc_mean', 1, 'd_mean', 1, 'f_mean', 1, 'dstar_mean', 1, ...
                'd95_gtvp', 1, 'dmean_gtvp', 1, 'lf', 0, 'lf_date', [], 'censor_date', [], 'total_time', [], 'total_follow_up_time', []);

            type_dir = fullfile(testCase.TempDir, 'Standard');
            mkdir(type_dir);
            save(fullfile(testCase.TempDir, 'dwi_vectors.mat'), 'data_vectors_gtvp', 'data_vectors_gtvn');
            save(fullfile(testCase.TempDir, 'summary_metrics.mat'), 'summary_metrics');
            save(fullfile(type_dir, 'summary_metrics_Standard.mat'), 'summary_metrics');

            % Run pipeline skipping 'load'
            % We run 'sanity' which should pass with minimal data or fail gracefully.
            % We catch any error, but check the output for the "Skipping Load" message.

            cmd = sprintf("run_dwi_pipeline('%s', {'sanity'})", testCase.ConfigPath);
            try
                T = evalc(cmd);
            catch
                % If sanity checks fail (e.g. alignment checks), that's fine,
                % we just want to verify loading logic.
                % But evalc might not capture output if error is thrown immediately.
                % So we should wrap run_dwi_pipeline to not error?
                % Or inspect the exception?
                % Actually, if error occurs, T is not assigned.
            end

            % Re-run capturing output by wrapping in try-catch block inside evalc string
            cmd_safe = sprintf("addpath('%s'); try, clear global MASTER_OUTPUT_FOLDER; run_dwi_pipeline('%s', {'sanity'}, '%s'); diary off; catch, diary off; end", testCase.RepoRoot, testCase.ConfigPath, testCase.TempDir);
            T = evalc(cmd_safe);

            testCase.verifyTrue(contains(T, 'Skipping Load Step'), 'Should skip load step');
            testCase.verifyTrue(contains(T, 'Loaded data from disk'), 'Should load data from disk');
        end

        function testSkipMetrics(testCase)
             type_dir = fullfile(testCase.TempDir, 'Standard');
             if ~exist(type_dir, 'dir'), mkdir(type_dir); end
             calculated_results = struct();
             save(fullfile(type_dir, 'calculated_results_Standard.mat'), 'calculated_results');

             % Also need dwi_vectors and summary_metrics for visualize arg list
             data_vectors_gtvp = struct('adc_vector', []);
             data_vectors_gtvn = struct();
             summary_metrics = struct('id_list', {{'Test'}}, 'mrn_list', {{'123'}}, ...
                'adc_mean', 1, 'd_mean', 1, 'f_mean', 1, 'dstar_mean', 1, ...
                'd95_gtvp', 1, 'dmean_gtvp', 1, 'lf', 0);
             save(fullfile(testCase.TempDir, 'dwi_vectors.mat'), 'data_vectors_gtvp', 'data_vectors_gtvn');
             save(fullfile(type_dir, 'summary_metrics_Standard.mat'), 'summary_metrics');

             cmd_safe = sprintf("addpath('%s'); try, clear global MASTER_OUTPUT_FOLDER; run_dwi_pipeline('%s', {'visualize'}, '%s'); diary off; catch, diary off; end", testCase.RepoRoot, testCase.ConfigPath, testCase.TempDir);
             T = evalc(cmd_safe);

             testCase.verifyTrue(contains(T, 'Skipping metrics_longitudinal'), 'Should skip metrics');
             testCase.verifyTrue(contains(T, 'Loaded calculated_results from disk'), 'Should load results from disk');
             testCase.verifyTrue(contains(T, 'Visualizing results'), 'Should run visualization');
        end

        function testLoadSavesSummary(testCase)
             % Test that running 'load' (with skip_to_reload=true) creates summary_metrics.mat

             % 1. Create dwi_vectors.mat with sufficient dummy data for load_dwi_data Section 5
             id_list = {'Test'};
             mrn_list = {'123'};
             fx_dates = {'20230101'};
             dwi_locations = cell(1,6,1);
             rtdose_locations = cell(1,6);
             gtv_locations = cell(1,6,1);
             gtvn_locations = cell(1,6,1);
             dmean_gtvp = ones(1,6); dmean_gtvn = ones(1,6);
             d95_gtvp = ones(1,6); d95_gtvn = ones(1,6);
             v50gy_gtvp = ones(1,6); v50gy_gtvn = ones(1,6);
             bad_dwi_locations = {}; bad_dwi_count = 0;
             lf = 0; immuno = 0;
             data_vectors_gtvn = struct();

             % data_vectors_gtvp(1,6,1)
             for i=1:6
                 data_vectors_gtvp(1,i,1).adc_vector = [1; 1];
                 data_vectors_gtvp(1,i,1).d_vector = [1; 1];
                 data_vectors_gtvp(1,i,1).f_vector = [0.1; 0.1];
                 data_vectors_gtvp(1,i,1).dstar_vector = [0.01; 0.01];
                 data_vectors_gtvp(1,i,1).vox_vol = 1;
                 data_vectors_gtvp(1,i,1).dose_vector = [1; 1];

                 % dnCNN and ivimnet fields
                 data_vectors_gtvp(1,i,1).adc_vector_dncnn = [1; 1];
                 data_vectors_gtvp(1,i,1).d_vector_dncnn = [1; 1];
                 data_vectors_gtvp(1,i,1).f_vector_dncnn = [0.1; 0.1];
                 data_vectors_gtvp(1,i,1).dstar_vector_dncnn = [0.01; 0.01];

                 data_vectors_gtvp(1,i,1).d_vector_ivimnet = [1; 1];
                 data_vectors_gtvp(1,i,1).f_vector_ivimnet = [0.1; 0.1];
                 data_vectors_gtvp(1,i,1).dstar_vector_ivimnet = [0.01; 0.01];
             end

             save(fullfile(testCase.TempDir, 'dwi_vectors.mat'), 'data_vectors_gtvn','data_vectors_gtvp','lf','immuno','mrn_list','id_list','fx_dates','dwi_locations','rtdose_locations','gtv_locations','gtvn_locations','dmean_gtvp','dmean_gtvn','d95_gtvp','d95_gtvn','v50gy_gtvp','v50gy_gtvn','bad_dwi_locations','bad_dwi_count');

             cmd_safe = sprintf("addpath('%s'); try, clear global MASTER_OUTPUT_FOLDER; run_dwi_pipeline('%s', {'load'}, '%s'); diary off; catch ME, diary off; disp(ME.message); end", testCase.RepoRoot, testCase.ConfigPath, testCase.TempDir);
             T = evalc(cmd_safe);

             fprintf('DEBUG T:\n%s\n', T);
             % Check if summary_metrics_Standard.mat exists in the isolated folder
             type_dir = fullfile(testCase.TempDir, 'Standard');
             testCase.verifyTrue(exist(fullfile(type_dir, 'summary_metrics_Standard.mat'), 'file') == 2, 'summary_metrics_Standard.mat should be created');
        end
    end
end
