classdef test_dwi_pipeline < matlab.unittest.TestCase
    % TESTDWIPIPELINE Integration test suite for the DWI Analysis Pipeline.
    %
    % This class uses the matlab.unittest framework to ensure that the
    % data loading and metric calculations produce valid, physical results,
    % avoiding workspace collisions and providing robust assertions.
    %
    % For algorithmic tests (BH, Holm, LOOCV, etc.), see test_statistical_methods.m
    % For static code analysis tests, see test_source_code_standards.m
    %
    % Run tests with:
    %   results = runtests('tests/test_dwi_pipeline.m');

    properties
        % Define properties that can be shared across tests if needed
        MockDataDir
        ConfigStruct
    end

    methods(TestMethodSetup)
        % Setup for each test
        function createMockConfig(testCase)
            % Suppress figure pop-ups during pipeline execution
            set(0, 'DefaultFigureVisible', 'off');

            % Create a dummy configuration structure for testing
            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.skip_to_reload = false;
            testCase.ConfigStruct.ivim_bthr = 100;
            testCase.ConfigStruct.dataloc = fullfile(pwd, 'mock_data');
            testCase.ConfigStruct.dcm2nii_call = 'dummy_dcm2niix';
            testCase.ConfigStruct.adc_thresh = 0.00115;
            testCase.ConfigStruct.high_adc_thresh = 0.001;
            testCase.ConfigStruct.d_thresh = 0.001;
            testCase.ConfigStruct.f_thresh = 0.1;
            testCase.ConfigStruct.dstar_thresh = 0.01;
            testCase.ConfigStruct.min_vox_hist = 10; % Lowered for small mock
            testCase.ConfigStruct.adc_max = 3.0e-3;
            testCase.ConfigStruct.clinical_data_sheet = 'mock_sheet.xlsx';
            testCase.ConfigStruct.patient_ids = {'P99-MOCK'};

            % Create a temporary directory for mock data
            testCase.MockDataDir = testCase.ConfigStruct.dataloc;
            if ~exist(testCase.MockDataDir, 'dir')
                mkdir(testCase.MockDataDir);
            end

            % Add necessary paths for tests
            testCase.ConfigStruct.orig_path = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(baseDir);  % repo root for run_dwi_pipeline
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods(TestMethodTeardown)
        % Cleanup after each test
        function removeMockData(testCase)
            % Restore path BEFORE deleting mock_data to avoid
            % "Removed ... from the MATLAB path" warnings.
            if isfield(testCase.ConfigStruct, 'orig_path')
                % Filter out entries that no longer exist on disk (e.g.
                % mock_data was deleted by a prior test in this run).
                entries = strsplit(testCase.ConfigStruct.orig_path, pathsep);
                entries = entries(cellfun(@(p) isempty(p) || isfolder(p), entries));
                path(strjoin(entries, pathsep));
            end

            if exist(testCase.MockDataDir, 'dir')
                rmdir(testCase.MockDataDir, 's');
            end

            set(0, 'DefaultFigureVisible', 'on');
        end
    end

    methods(Access = private)
        function [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = generateValidMockData(testCase)
            % Helper function to generate mock data structures that pass sanity_checks
            % and have sufficient dimensionality for downstream metrics

            nPat = 3;
            nTp = 4; % Need at least 2 timepoints for longitudinal, say Fx1, Fx2, Fx3, Post
            nDwiType = 3; % Standard, DnCNN, IVIMnet

            % 1. Create Data Vectors
            % Generate dummy voxel data
            nvox = 50;
            base_adc = 0.001;
            base_d = 0.0008;
            base_f = 0.15;
            base_dstar = 0.05;
            base_dose = 2.0;

            % Pre-allocate struct arrays with matching field signature to
            % avoid 'heterogeneousStrucAssignment' error on indexed assignment.
            s_template = struct( ...
                'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
                'adc_vector_dncnn', [], 'd_vector_dncnn', [], 'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
                'd_vector_ivimnet', [], 'f_vector_ivimnet', [], 'dstar_vector_ivimnet', [], ...
                'dose_vector', [], 'ID', '', 'MRN', '', 'LF', 0, 'Immuno', 0, ...
                'Fraction', 0, 'Repeatability_index', 0, 'vox_vol', 0, 'dvh', [], 'd95', 0, 'v50gy', 0);
            data_vectors_gtvp = repmat(s_template, nPat, nTp);
            data_vectors_gtvn = repmat(s_template, nPat, nTp);

            for j = 1:nPat
                for k = 1:nTp
                    for d_type = 1:nDwiType
                        for rpi = 1:1 % Just 1 repeat
                            % Add some noise, keep values physical
                            v_adc = abs(base_adc + 1e-4 * randn(nvox, 1));
                            v_d = abs(base_d + 1e-4 * randn(nvox, 1));
                            v_f = abs(base_f + 0.02 * randn(nvox, 1));
                            v_dstar = abs(base_dstar + 0.01 * randn(nvox, 1));
                            v_dose = abs(base_dose * k + 0.5 * randn(nvox, 1));

                            % Make sure no exact zeros
                            v_adc(v_adc == 0) = 1e-5;
                            v_d(v_d == 0) = 1e-5;

                            s_p = struct();
                            s_p.adc_vector = v_adc;
                            s_p.d_vector = v_d;
                            s_p.f_vector = v_f;
                            s_p.dstar_vector = v_dstar;
                            s_p.adc_vector_dncnn = v_adc;
                            s_p.d_vector_dncnn = v_d;
                            s_p.f_vector_dncnn = v_f;
                            s_p.dstar_vector_dncnn = v_dstar;

                            s_p.d_vector_ivimnet = v_d;
                            s_p.f_vector_ivimnet = v_f;
                            s_p.dstar_vector_ivimnet = v_dstar;
                            s_p.dose_vector = v_dose;

                            s_p.ID = sprintf('P%02d', j);
                            s_p.MRN = sprintf('M%02d', j);
                            s_p.LF = mod(j, 2); % Alternating LF
                            s_p.Immuno = mod(j+1, 2);
                            s_p.Fraction = k;
                            s_p.Repeatability_index = rpi;
                            s_p.vox_vol = 0.02;
                            s_p.dvh = rand(100,2);
                            s_p.d95 = base_dose * k * 0.9;
                            s_p.v50gy = 10;

                            data_vectors_gtvp(j, k, rpi) = s_p;

                            % Same for GTVn
                            s_n = s_p;
                            s_n.adc_vector = v_adc * 1.1;
                            data_vectors_gtvn(j, k, rpi) = s_n;
                        end
                    end
                end
            end

            % 2. Create Summary Metrics
            summary_metrics = struct();
            summary_metrics.id_list = arrayfun(@(x) sprintf('P%02d', x), 1:nPat, 'UniformOutput', false);
            summary_metrics.mrn_list = arrayfun(@(x) sprintf('M%02d', x), 1:nPat, 'UniformOutput', false);

            summary_metrics.adc_mean = base_adc + 1e-4 * randn(nPat, nTp, 3);
            summary_metrics.d_mean = base_d + 1e-4 * randn(nPat, nTp, nDwiType);
            summary_metrics.f_mean = base_f + 0.02 * randn(nPat, nTp, nDwiType);
            summary_metrics.dstar_mean = base_dstar + 0.01 * randn(nPat, nTp, nDwiType);

            summary_metrics.d95_gtvp = base_dose * repmat(1:nTp, nPat, 1) * 0.9;
            summary_metrics.d95_gtvn = base_dose * repmat(1:nTp, nPat, 1) * 0.8;
            summary_metrics.v50gy_gtvp = 10 + rand(nPat, nTp);
            summary_metrics.v50gy_gtvn = 5 + rand(nPat, nTp);

            summary_metrics.lf = mod((1:nPat)', 2);
            summary_metrics.dmean_gtvp = base_dose * repmat(1:nTp, nPat, 1);
            summary_metrics.gtv_vol = 10 + rand(nPat, nTp);
            summary_metrics.gtv_locations = cell(nPat, nTp, 1);
            summary_metrics.dwi_locations = cell(nPat, nTp, 1);
            summary_metrics.adc_sd = 1e-4 * rand(nPat, nTp, nDwiType);

            % Repeatability fields (nPat x nRpt x nDwiType)
            nRpt = 2;
            summary_metrics.adc_mean_rpt = abs(base_adc + 1e-4 * randn(nPat, nRpt, nDwiType)) + 1e-5;
            summary_metrics.adc_sub_rpt = abs(base_adc + 1e-4 * randn(nPat, nRpt, nDwiType)) + 1e-5;
            summary_metrics.d_mean_rpt = abs(base_d + 1e-4 * randn(nPat, nRpt, nDwiType)) + 1e-5;
            summary_metrics.f_mean_rpt = abs(base_f + 0.02 * randn(nPat, nRpt, nDwiType)) + 1e-5;
            summary_metrics.dstar_mean_rpt = abs(base_dstar + 0.01 * randn(nPat, nRpt, nDwiType)) + 1e-5;
            summary_metrics.n_rpt = nRpt * ones(nPat, 1);

            % Make sure all metrics are positive to pass sanity checks
            summary_metrics.adc_mean = abs(summary_metrics.adc_mean) + 1e-5;
            summary_metrics.d_mean = abs(summary_metrics.d_mean) + 1e-5;
            summary_metrics.f_mean = abs(summary_metrics.f_mean) + 1e-5;
            summary_metrics.dstar_mean = abs(summary_metrics.dstar_mean) + 1e-5;

            % Repeatability fields required by metrics_baseline
            nRpt = 2;
            summary_metrics.n_rpt = nRpt * ones(nPat, 1);
            summary_metrics.adc_mean_rpt = abs(base_adc + 1e-4 * randn(nPat, nRpt, nDwiType));
            summary_metrics.adc_sub_rpt = abs(base_adc + 1e-4 * randn(nPat, nRpt, nDwiType));
            summary_metrics.d_mean_rpt = abs(base_d + 1e-4 * randn(nPat, nRpt, nDwiType));
            summary_metrics.f_mean_rpt = abs(base_f + 0.02 * randn(nPat, nRpt, nDwiType));
            summary_metrics.dstar_mean_rpt = abs(base_dstar + 0.01 * randn(nPat, nRpt, nDwiType));

            % Volume and dose fields
            summary_metrics.gtv_vol = 10 + rand(nPat, nTp);
            summary_metrics.dmean_gtvp = 50 + rand(nPat, nTp) * 5;

            % LF field for visualize_results
            summary_metrics.lf = mod((1:nPat)', 2);
        end
    end

    methods(Test)

        function testRunDwiPipelineEndToEnd(testCase)
            % Create mock data
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = testCase.generateValidMockData();

            % Create mock config.json
            config_file = fullfile(testCase.MockDataDir, 'config.json');

            % Setup configuration struct to be written as json
            cfg = struct();
            cfg.dataloc = testCase.MockDataDir;
            cfg.dcm2nii_call = 'dummy';
            cfg.clinical_data_sheet = 'mock_clinical.xlsx';
            cfg.skip_to_reload = true;
            cfg.dwi_type = 'IVIMnet';

            % Write json
            fid = fopen(config_file, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            % Create mock clinical Excel file matching patient IDs
            nPat = 3;
            id_list_xls = arrayfun(@(x) sprintf('P%02d', x), (1:nPat)', 'UniformOutput', false);
            lf_vals = mod((1:nPat)', 2);
            dt_event  = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop  = repmat(datetime('2022-03-01'), nPat, 1);
            dt_reg_censor = repmat(datetime('2023-06-01'), nPat, 1);
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, ...
                dt_reg_censor, dt_rtstart, dt_rtstop, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate'});
            writetable(T_clin, fullfile(testCase.MockDataDir, 'mock_clinical.xlsx'));

            % Create mock dwi_vectors_IVIMnet.mat and summary_metrics_IVIMnet.mat
            % in the locations expected by run_dwi_pipeline

            % The pipeline uses output_folder for summary_metrics
            output_folder = fullfile(testCase.MockDataDir, 'IVIMnet');
            if ~exist(output_folder, 'dir')
                mkdir(output_folder);
            end

            % Save dummy data vectors to dataloc
            dwi_vectors_file = fullfile(testCase.MockDataDir, 'dwi_vectors_IVIMnet.mat');
            save(dwi_vectors_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            % Save summary metrics to output_folder
            summary_metrics_file = fullfile(output_folder, 'summary_metrics_IVIMnet.mat');
            save(summary_metrics_file, 'summary_metrics');

            % Also save the legacy format just in case
            legacy_dwi_file = fullfile(testCase.MockDataDir, 'dwi_vectors.mat');
            save(legacy_dwi_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            % Define steps to run, skipping 'load' since we mocked the files directly
            steps_to_run = {'sanity', 'visualize', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};

            % Keep original directory
            orig_dir = pwd;

            % Try to run the pipeline
            try
                % The function might use global variable MASTER_OUTPUT_FOLDER
                clear global MASTER_OUTPUT_FOLDER;
                % Skip preflight test suite to avoid recursive test execution
                setenv('SKIP_PIPELINE_PREFLIGHT', '1');
                % Change into MockDataDir since some functions might rely on pwd
                cd(testCase.MockDataDir);
                % Call the master orchestrator function
                run_dwi_pipeline('config.json', steps_to_run, testCase.MockDataDir);
                passed = true;
                cd(orig_dir);
                setenv('SKIP_PIPELINE_PREFLIGHT', '');
            catch ME
                cd(orig_dir);
                setenv('SKIP_PIPELINE_PREFLIGHT', '');
                if exist('OCTAVE_VERSION', 'builtin')
                    fprintf('Error: %s\n', ME.message);
                else
                    disp(ME.getReport());
                end
                passed = false;
            end

            % Assert pipeline ran successfully
            testCase.verifyTrue(passed, 'run_dwi_pipeline failed to execute end-to-end with mock data.');

            % Verify expected outputs were generated
            sanity_results = fullfile(output_folder, 'sanity_checks_results_IVIMnet.txt');
            testCase.verifyTrue(exist(sanity_results, 'file') == 2, 'Sanity check results file not found.');

            visualize_results = fullfile(output_folder, 'visualize_results_state_IVIMnet.txt');
            testCase.verifyTrue(exist(visualize_results, 'file') == 2, 'Visualize results file not found.');

            baseline_results = fullfile(output_folder, 'metrics_baseline_results_IVIMnet.mat');
            testCase.verifyTrue(exist(baseline_results, 'file') == 2, 'Baseline results file not found.');

            longitudinal_results = fullfile(output_folder, 'metrics_longitudinal_results_IVIMnet.txt');
            testCase.verifyTrue(exist(longitudinal_results, 'file') == 2, 'Longitudinal results file not found.');

            dosimetry_results = fullfile(output_folder, 'metrics_dosimetry_results_IVIMnet.mat');
            testCase.verifyTrue(exist(dosimetry_results, 'file') == 2, 'Dosimetry results file not found.');

            comparisons_results = fullfile(output_folder, 'metrics_stats_comparisons_results_IVIMnet.txt');
            testCase.verifyTrue(exist(comparisons_results, 'file') == 2, 'Comparisons results file not found.');

            predictive_results = fullfile(output_folder, 'metrics_stats_predictive_results_IVIMnet.mat');
            testCase.verifyTrue(exist(predictive_results, 'file') == 2, 'Predictive results file not found.');

            survival_results = fullfile(output_folder, 'metrics_survival_results_IVIMnet.txt');
            testCase.verifyTrue(exist(survival_results, 'file') == 2, 'Survival results file not found.');
        end

        function testMockNiftiPipelineOutputs(testCase)
            % -----------------------------------------------------------------
            % BOILERPLATE MOCK 3D NIFTI ARRAY TEST
            % -----------------------------------------------------------------
            % 1. Create a tiny mock 4D DWI tensor (X, Y, Z, b-values)
            % We simulate a 5x5x5 volume with 4 b-values (0, 30, 150, 550)
            sz = [5, 5, 5, 4];

            % Generate physically plausible exponential decay signal
            bvals = [0, 30, 150, 550];
            true_adc = 0.001; % Typical tissue ADC

            mockSignal = zeros(sz);
            for b = 1:4
                mockSignal(:,:,:,b) = 1500 * exp(-bvals(b) * true_adc) + 10*randn(5,5,5);
            end
            mockSignal(mockSignal < 1) = 1; % Prevent negative noise from causing log(0)

            % Generate a corresponding 3D GTV mask
            mockMask = zeros(5, 5, 5);
            mockMask(2:4, 2:4, 2:4) = 1; % Active central core

            % -----------------------------------------------------------------
            % 2. Execute a simulated pipeline segment
            % Normally, you would pass these into your load_dwi_data or fit_adc_mono
            % algorithms. Here we simulate the monoexponential ADC fit logic inline
            % to demonstrate the assertion framework.

            S_2d = reshape(mockSignal, [prod(sz(1:3)), length(bvals)]);
            adc_vec = zeros(numel(mockMask), 1);

            % Simple OLS log-linear fit for demonstration
            S_a = S_2d;
            % (-b(2:end) \ log(S(b>0)/S(b=0)))
            adc_vec = (-bvals(2:end)' \ log(S_a(:,2:end) ./ S_a(:,1))')';

            % Extract GTV only
            final_adc_features = adc_vec(mockMask == 1);

            % -----------------------------------------------------------------
            % 3. Automated Assertions using verifyTrue / verifyEqual

            % Is the output physically valid? (No negative ADCs)
            testCase.verifyTrue(all(final_adc_features > 0), ...
                'Negative ADC values calculated. Model fit is non-physical.');

            % Are there any NaNs generated during the process?
            testCase.verifyFalse(any(isnan(final_adc_features)), ...
                'NaNs detected in the output metric arrays.');

            % Does it accurately approximate our strictly physical ground truth?
            mean_calc_adc = mean(final_adc_features);
            testCase.verifyEqual(mean_calc_adc, true_adc, 'RelTol', 0.15, ...
                'Calculated ADC deviates from physical truth by >15%.');

        end

    end
end
