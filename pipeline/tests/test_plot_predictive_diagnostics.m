classdef test_plot_predictive_diagnostics < matlab.unittest.TestCase
% TEST_PLOT_PREDICTIVE_DIAGNOSTICS — Unit tests for plot_predictive_diagnostics.m
%
% Validates diagnostic figure generation including:
%   - ROC curve creation and file output
%   - Sanity check panels (volume, ADC SD, signal vs noise)
%   - 2D scatter plots with decision boundary
%   - Graceful handling of insufficient data
%   - Output file naming (asterisk → 'star' sanitization)

    properties
        OutputFolder
        OldFigVis
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));
        end

        function setupOutputFolder(testCase)
            testCase.OutputFolder = fullfile(tempdir, ['test_pred_diag_' char(java.util.UUID.randomUUID())]);
            mkdir(testCase.OutputFolder);
            % Suppress figure display during tests
            testCase.OldFigVis = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods (TestMethodTeardown)
        function cleanupOutputFolder(testCase)
            set(0, 'DefaultFigureVisible', testCase.OldFigVis);
            if exist(testCase.OutputFolder, 'dir')
                rmdir(testCase.OutputFolder, 's');
            end
        end
    end

    methods (Test)
        function test_roc_curve_saved(testCase)
            % ROC curve PNG should be created when valid data is provided
            [args] = testCase.makeValidArgs(2, testCase.OutputFolder);

            plot_predictive_diagnostics(args{:});

            roc_file = fullfile(testCase.OutputFolder, 'ROC_OOF_Risk_Score_Fx3_Standard.png');
            testCase.verifyTrue(exist(roc_file, 'file') > 0, ...
                'ROC curve PNG should be saved.');
        end

        function test_sanity_check_panels_saved(testCase)
            % One sanity check PNG per selected feature
            [args] = testCase.makeValidArgs(2, testCase.OutputFolder);

            plot_predictive_diagnostics(args{:});

            % Check that sanity check files exist for each feature
            files = dir(fullfile(testCase.OutputFolder, 'Sanity_Checks_*.png'));
            testCase.verifyGreaterThanOrEqual(numel(files), 1, ...
                'At least one sanity check PNG should be saved.');
        end

        function test_2d_scatter_with_two_features(testCase)
            % With 2+ features, 2D scatter plots should be generated
            [args] = testCase.makeValidArgs(2, testCase.OutputFolder);

            plot_predictive_diagnostics(args{:});

            files = dir(fullfile(testCase.OutputFolder, '2D_Space_*.png'));
            testCase.verifyGreaterThanOrEqual(numel(files), 1, ...
                '2D scatter PNG should be generated with 2+ features.');
        end

        function test_no_2d_scatter_with_one_feature(testCase)
            % With only 1 feature, 2D scatter should be skipped
            [args] = testCase.makeValidArgs(1, testCase.OutputFolder);

            plot_predictive_diagnostics(args{:});

            files = dir(fullfile(testCase.OutputFolder, '2D_Space_*.png'));
            testCase.verifyEqual(numel(files), 0, ...
                '2D scatter should not be generated with only 1 feature.');
        end

        function test_asterisk_sanitized_in_filename(testCase)
            % D* feature name should produce 'Dstar' in filename, not '*'
            [args] = testCase.makeValidArgs(1, testCase.OutputFolder);
            % Override feature name to include asterisk
            args{6} = {'D*'};  % sig_names

            plot_predictive_diagnostics(args{:});

            files = dir(fullfile(testCase.OutputFolder, 'Sanity_Checks_*Dstar*.png'));
            testCase.verifyGreaterThanOrEqual(numel(files), 1, ...
                'Asterisk in feature name should be replaced with "star".');
        end

        function test_no_crash_with_empty_roc_data(testCase)
            % When all risk scores are NaN, should print message but not crash
            [args] = testCase.makeValidArgs(1, testCase.OutputFolder);
            % Set all risk scores to NaN
            args{12} = nan(size(args{12}));  % risk_scores_all_target

            plot_predictive_diagnostics(args{:});

            % If we got here without error, the test passes
            testCase.verifyTrue(true, 'Should handle NaN risk scores gracefully.');
        end

        function test_no_open_figures_after_completion(testCase)
            % All figures should be closed after function completes
            figs_before = findall(0, 'Type', 'figure');
            n_before = numel(figs_before);

            [args] = testCase.makeValidArgs(2, testCase.OutputFolder);
            plot_predictive_diagnostics(args{:});

            figs_after = findall(0, 'Type', 'figure');
            testCase.verifyEqual(numel(figs_after), n_before, ...
                'No new figures should remain open after completion.');
        end
    end

    methods (Static, Access = private)
        function [args] = makeValidArgs(n_sig, output_folder)
            % Generate a minimal valid argument set for plot_predictive_diagnostics
            n_pts = 30;
            rng(42);

            % Build feature data
            sig_data_selected = cell(1, n_sig);
            sig_names = cell(1, n_sig);
            sig_is_abs = false(1, n_sig);
            sig_is_pct_imaging = true(1, n_sig);
            sig_disp_names = cell(1, n_sig);
            sig_units = cell(1, n_sig);
            sig_col_idx = ones(1, n_sig);
            sig_abs_data = cell(1, n_sig);
            sig_pct_data = cell(1, n_sig);
            selected_indices = zeros(1, n_sig);

            for i = 1:n_sig
                sig_data_selected{i} = randn(n_pts, 3);
                sig_names{i} = sprintf('ADC_feat%d', i);
                sig_disp_names{i} = sprintf('\\Delta ADC %d', i);
                sig_units{i} = '%';
                sig_abs_data{i} = randn(n_pts, 3);
                sig_pct_data{i} = randn(n_pts, 3);
                selected_indices(i) = i;  % base_idx will be 1 (ADC)
            end

            lf_group = zeros(n_pts, 1);
            lf_group(1:10) = 1;         % local failure
            lf_group(11:15) = 0;        % local control
            lf_group(16:25) = 0;
            lf_group(26:30) = 2;        % competing risk (excluded from ROC)

            risk_scores = randn(n_pts, 1);
            valid_pts = true(n_pts, 1);
            m_gtv_vol = abs(randn(n_pts, 3)) * 10 + 1;
            adc_sd = abs(randn(n_pts, 3, 3)) * 0.001;
            ADC_abs = abs(randn(n_pts, 3)) * 0.001 + 0.001;

            target_fx = 2;
            fx_label = 'Fx3';
            dtype_label = 'Standard';
            dtype = 1;
            use_firth = false;

            args = {selected_indices, n_sig, sig_data_selected, sig_names, ...
                sig_is_abs, sig_is_pct_imaging, sig_disp_names, sig_units, ...
                sig_col_idx, sig_abs_data, sig_pct_data, ...
                risk_scores, lf_group, valid_pts, ...
                m_gtv_vol, adc_sd, ADC_abs, ...
                target_fx, fx_label, dtype_label, dtype, output_folder, use_firth};
        end
    end
end
