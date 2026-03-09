classdef test_cross_dwi_subvolume < matlab.unittest.TestCase
    % TEST_CROSS_DWI_SUBVOLUME Tests for plot_cross_dwi_subvolume_comparison.
    %
    % plot_cross_dwi_subvolume_comparison generates scatter plots comparing
    % ADC subvolume percentage across DWI processing types (Standard, dnCNN,
    % IVIMnet). It loads saved summary_metrics checkpoint files from disk
    % and overlays data from all available types on a single figure.
    %
    % Tests verify:
    %   1. Two types available -> figure is created
    %   2. Only one type -> graceful skip (no figure)
    %   3. Repeat scan data -> figure includes a repeat-scan panel
    %   4. Standard + IVIMnet (no dnCNN) -> figure is created
    %   5. All-NaN subvolume data -> graceful skip (no figure)
    %
    % Run tests with:
    %   results = runtests('tests/test_cross_dwi_subvolume.m');

    properties
        TempDir       % Temporary directory for mock checkpoint files
        ConfigStruct  % Mock pipeline configuration
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Configure paths to point at the temp directory where mock
            % checkpoint .mat files will be saved.
            testCase.ConfigStruct.dataloc = testCase.TempDir;
            testCase.ConfigStruct.output_folder = testCase.TempDir;
            testCase.ConfigStruct.master_output_folder = testCase.TempDir;

            pancDataPath = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(pancDataPath, 'core'));
            addpath(fullfile(pancDataPath, 'utils'));

            % Suppress figure windows during automated testing
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function cleanUp(testCase)
            % Turn off diary (module may have opened one) before rmdir
            % to avoid file-lock issues on Windows.
            diary off;
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testTwoTypesProducesFigure(testCase)
            % With both Standard and dnCNN checkpoint files present,
            % the function should generate a scatter plot comparing ADC
            % subvolume percentages across the two types.
            nPat = 10;
            sm_std = make_mock_sm(nPat);
            % Populate Fx1 (timepoint 1) subvolume for DWI type 1 (Standard)
            sm_std.adc_sub_vol_pc(:,1,1) = linspace(0.2, 0.5, nPat)';
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_Standard.mat'), sm_std);

            sm_dncnn = make_mock_sm(nPat);
            % Populate Fx1 subvolume for DWI type 2 (dnCNN) with slightly different values
            sm_dncnn.adc_sub_vol_pc(:,1,2) = linspace(0.15, 0.45, nPat)';
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_dnCNN.mat'), sm_dncnn);

            plot_cross_dwi_subvolume_comparison(sm_std, testCase.ConfigStruct);

            out_file = fullfile(testCase.TempDir, 'Cross_DWI_ADC_Subvolume_Fx1.png');
            testCase.verifyTrue(exist(out_file, 'file') > 0, 'Expected output PNG to be created.');
        end

        function testSingleTypeSkips(testCase)
            % Cross-DWI comparison requires at least 2 DWI types. When
            % only Standard is available (no dnCNN or IVIMnet checkpoint),
            % the function should skip gracefully without creating a figure.
            nPat = 5;
            sm_std = make_mock_sm(nPat);
            sm_std.adc_sub_vol_pc(:,1,1) = 0.3 * ones(nPat, 1);
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_Standard.mat'), sm_std);

            plot_cross_dwi_subvolume_comparison(sm_std, testCase.ConfigStruct);

            out_file = fullfile(testCase.TempDir, 'Cross_DWI_ADC_Subvolume_Fx1.png');
            testCase.verifyFalse(exist(out_file, 'file') > 0, 'Should not create figure with only 1 type.');
        end

        function testRepeatPanelShown(testCase)
            % Checkpoints with repeat scan data populated
            nPat = 8;
            nRpt = 3;

            sm_std = make_mock_sm(nPat, nRpt);
            sm_std.adc_sub_vol_pc(:,1,1) = linspace(0.2, 0.5, nPat)';
            sm_std.adc_sub_vol_pc_rpt = nan(nPat, nRpt, 3);
            sm_std.adc_sub_vol_pc_rpt(:,1,1) = linspace(0.20, 0.50, nPat)';
            sm_std.adc_sub_vol_pc_rpt(:,2,1) = linspace(0.18, 0.48, nPat)';
            sm_std.adc_sub_vol_pc_rpt(:,3,1) = linspace(0.22, 0.52, nPat)';
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_Standard.mat'), sm_std);

            sm_dncnn = make_mock_sm(nPat, nRpt);
            sm_dncnn.adc_sub_vol_pc(:,1,2) = linspace(0.15, 0.45, nPat)';
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_dnCNN.mat'), sm_dncnn);

            plot_cross_dwi_subvolume_comparison(sm_std, testCase.ConfigStruct);

            out_file = fullfile(testCase.TempDir, 'Cross_DWI_ADC_Subvolume_Fx1.png');
            testCase.verifyTrue(exist(out_file, 'file') > 0, 'Expected output PNG with repeat panel.');
        end

        function testIdenticalStandardIVIMnet(testCase)
            % Standard and IVIMnet with identical ADC subvolume
            nPat = 6;
            shared_vals = linspace(0.2, 0.4, nPat)';

            sm_std = make_mock_sm(nPat);
            sm_std.adc_sub_vol_pc(:,1,1) = shared_vals;
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_Standard.mat'), sm_std);

            sm_ivimnet = make_mock_sm(nPat);
            sm_ivimnet.adc_sub_vol_pc(:,1,3) = shared_vals;
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_IVIMnet.mat'), sm_ivimnet);

            % Should run without error and create the figure
            plot_cross_dwi_subvolume_comparison(sm_std, testCase.ConfigStruct);

            out_file = fullfile(testCase.TempDir, 'Cross_DWI_ADC_Subvolume_Fx1.png');
            testCase.verifyTrue(exist(out_file, 'file') > 0, 'Expected output PNG for Standard+IVIMnet.');
        end

        function testAllNaN(testCase)
            % All-NaN subvolume data — should skip gracefully
            nPat = 5;
            sm_std = make_mock_sm(nPat);
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_Standard.mat'), sm_std);

            sm_dncnn = make_mock_sm(nPat);
            save_mock(fullfile(testCase.TempDir, 'summary_metrics_dnCNN.mat'), sm_dncnn);

            plot_cross_dwi_subvolume_comparison(sm_std, testCase.ConfigStruct);

            out_file = fullfile(testCase.TempDir, 'Cross_DWI_ADC_Subvolume_Fx1.png');
            testCase.verifyFalse(exist(out_file, 'file') > 0, 'Should not create figure with all-NaN data.');
        end
    end
end

function sm = make_mock_sm(nPat, nRpt)
% Create a minimal mock summary_metrics struct (not wrapped).
    if nargin < 2, nRpt = 1; end
    nTp = 6;
    sm.id_list = arrayfun(@(x) sprintf('P%02d', x), 1:nPat, 'UniformOutput', false);
    sm.adc_sub_vol_pc = nan(nPat, nTp, 3);
    sm.adc_sub_vol = nan(nPat, nTp, 3);
    sm.adc_mean_rpt = nan(nPat, nRpt, 3);
    sm.adc_sub_vol_rpt = nan(nPat, nRpt, 3);
end

function save_mock(filepath, sm)
% Save sm as 'summary_metrics' variable in a .mat file.
    summary_metrics = sm; %#ok<NASGU>
    save(filepath, 'summary_metrics');
end
