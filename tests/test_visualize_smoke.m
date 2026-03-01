classdef test_visualize_smoke < matlab.unittest.TestCase
    % TEST_VISUALIZE_SMOKE Smoke tests for visualize_results
    properties
        TempDir
        ConfigStruct
        SummaryMetrics
        DataVectors
        CalculatedResults
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Create temp dir
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Setup patient structure
            patID = 'P01';
            niiDir = fullfile(testCase.TempDir, patID, 'nii');
            mkdir(niiDir);

            % Create valid 4D NIfTI for DWI
            % dimensions: 10x10x10x4, data type: double
            dwi_img = zeros(10, 10, 10, 4, 'double');
            dwi_img(:,:,:,1) = 1000; % b=0
            dwi_img(:,:,:,2) = 800;  % b=30
            dwi_img(:,:,:,3) = 500;  % b=150
            dwi_img(:,:,:,4) = 100;  % b=550
            dwi_file = fullfile(niiDir, 'fx1_dwi1.nii');
            niftiwrite(dwi_img, dwi_file);
            % compress to .nii.gz
            gzip(dwi_file);
            delete(dwi_file);

            % Create valid 3D NIfTI for GTV
            gtv_img = zeros(10, 10, 10, 'double');
            gtv_img(4:7, 4:7, 4:7) = 1; % simple central box
            gtv_file = fullfile(niiDir, 'fx1_gtv1.nii');
            niftiwrite(gtv_img, gtv_file);
            gzip(gtv_file);
            delete(gtv_file);

            % Create valid .bval
            fid = fopen(fullfile(niiDir, 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 30 150 550');
            fclose(fid);

            % Setup inputs
            testCase.ConfigStruct.dataloc = testCase.TempDir;
            testCase.ConfigStruct.output_folder = fullfile(testCase.TempDir, 'saved_figures');
            testCase.ConfigStruct.dwi_types_to_run = 1; % Run 'Standard'

            testCase.SummaryMetrics.id_list = {patID};
            testCase.SummaryMetrics.mrn_list = {'MRN01'};
            testCase.SummaryMetrics.lf = [0];

            % Initialize summary arrays (1 valid patient, 1 timepoint, 1 dwi_type)
            testCase.SummaryMetrics.adc_mean = 1.0e-3 * ones(1,1,1);
            testCase.SummaryMetrics.d_mean = 1.0e-3 * ones(1,1,1);
            testCase.SummaryMetrics.f_mean = 0.1 * ones(1,1,1);
            testCase.SummaryMetrics.dstar_mean = 0.05 * ones(1,1,1);
            testCase.SummaryMetrics.d95_gtvp = 40 * ones(1,1);
            testCase.SummaryMetrics.dmean_gtvp = 50 * ones(1,1);

            testCase.CalculatedResults = struct();

            % Setup DataVectors
            testCase.DataVectors = struct('adc_vector', {ones(10,1)});

            % Add paths
            addpath(fullfile(pwd, 'core'));
            addpath(fullfile(pwd, 'utils'));
            addpath(fullfile(pwd, 'dependencies'));

            % Ensure no figures are visible during testing
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            % Close any left-over figures
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testSmokeValidData(testCase)
            % Should run without error and generate all files
            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(exist(fullfile(outputDir, 'Parameter_Maps_1.png'), 'file') > 0, ...
                'Expected Parameter_Maps_1.png to be created');
            testCase.verifyTrue(exist(fullfile(outputDir, 'Feature_Histograms_Standard.png'), 'file') > 0, ...
                'Expected Feature_Histograms_Standard.png to be created');
            testCase.verifyTrue(exist(fullfile(outputDir, 'Feature_BoxPlots_Standard.png'), 'file') > 0, ...
                'Expected Feature_BoxPlots_Standard.png to be created');
            testCase.verifyTrue(exist(fullfile(outputDir, 'Dose_vs_Diffusion_Standard.png'), 'file') > 0, ...
                'Expected Dose_vs_Diffusion_Standard.png to be created');
        end

        function testSmokeMissingBval(testCase)
            % Remove bval file
            delete(fullfile(testCase.TempDir, 'P01', 'nii', 'fx1_dwi1.bval'));

            % Should not crash, but skip Parameter Maps (since file missing)
            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(exist(fullfile(outputDir, 'Parameter_Maps_1.png'), 'file') == 0, ...
                'Parameter_Maps_1.png should NOT be created if bval is missing');

            % Section 2 & 3 still run because they rely on SummaryMetrics, not file loading
            testCase.verifyTrue(exist(fullfile(outputDir, 'Feature_Histograms_Standard.png'), 'file') > 0, ...
                'Expected Feature_Histograms_Standard.png to be created');
        end

        function testSmokeProtocolDeviation(testCase)
            % Write incorrect bval
            fid = fopen(fullfile(testCase.TempDir, 'P01', 'nii', 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 50 150 550');
            fclose(fid);

            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(exist(fullfile(outputDir, 'Parameter_Maps_1.png'), 'file') == 0, ...
                'Parameter_Maps_1.png should NOT be created if protocol deviates');
        end
    end
end
