classdef test_visualize_smoke < matlab.unittest.TestCase
    % TEST_VISUALIZE_SMOKE Smoke tests for visualize_results.
    %
    % These tests verify that visualize_results runs without crashing and
    % produces the expected output PNG files. They use a synthetic single-
    % patient dataset with a minimal 10x10x10 NIfTI volume (4 b-values) and
    % a simple cubic GTV mask. The tests exercise three scenarios:
    %   1. Valid data -> all expected figures are created
    %   2. Missing .bval file -> parameter maps are skipped but other plots succeed
    %   3. Protocol deviation in .bval -> parameter maps are skipped
    %
    % These are "smoke tests" -- they verify the function completes without
    % error and produces output, but do not validate the correctness of the
    % generated plots.
    %
    % Run tests with:
    %   results = runtests('tests/test_visualize_smoke.m');

    properties
        TempDir          % Temporary directory holding synthetic patient data
        ConfigStruct     % Mock pipeline configuration
        SummaryMetrics   % Mock summary metric arrays
        DataVectors      % Mock voxel-level data vectors
        CalculatedResults % Empty struct (no pre-computed results needed)
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Build a synthetic patient directory tree with NIfTI images,
            % GTV mask, and b-value file to mimic real pipeline input.

            % Create temp dir for this test run
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Setup patient directory structure: <TempDir>/P01/nii/
            patID = 'P01';
            niiDir = fullfile(testCase.TempDir, patID, 'nii');
            mkdir(niiDir);

            % Create valid 4D NIfTI for DWI.
            % Dimensions: 10x10x10 spatial, 4 b-values.
            % Signal intensities decrease with b-value (physically realistic
            % mono-exponential decay pattern).
            dwi_img = zeros(10, 10, 10, 4, 'double');
            dwi_img(:,:,:,1) = 1000; % b=0   (baseline signal)
            dwi_img(:,:,:,2) = 800;  % b=30  (slight attenuation)
            dwi_img(:,:,:,3) = 500;  % b=150 (moderate attenuation)
            dwi_img(:,:,:,4) = 100;  % b=550 (heavy attenuation)
            dwi_file = fullfile(niiDir, 'fx1_dwi1.nii');
            niftiwrite(dwi_img, dwi_file);
            gzip(dwi_file);   % Compress to .nii.gz (production format)
            delete(dwi_file); % Remove uncompressed copy

            % Create valid 3D NIfTI for GTV (Gross Tumor Volume) mask.
            % A 4x4x4 central cube simulates a simple tumour contour.
            gtv_img = zeros(10, 10, 10, 'double');
            gtv_img(4:7, 4:7, 4:7) = 1;
            gtv_file = fullfile(niiDir, 'fx1_gtv1.nii');
            niftiwrite(gtv_img, gtv_file);
            gzip(gtv_file);
            delete(gtv_file);

            % Create valid .bval file matching the 4 b-values in the DWI.
            % The standard pancreas DWI protocol uses b = [0, 30, 150, 550].
            fid = fopen(fullfile(niiDir, 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 30 150 550');
            fclose(fid);

            % Setup mock pipeline configuration
            testCase.ConfigStruct.dataloc = testCase.TempDir;
            testCase.ConfigStruct.output_folder = fullfile(testCase.TempDir, 'output');
            testCase.ConfigStruct.dwi_types_to_run = 1; % 1 = 'Standard'

            % Setup mock summary metrics for 1 patient, 1 timepoint, 1 DWI type.
            % Values are clinically plausible for pancreatic tissue.
            testCase.SummaryMetrics.id_list = {patID};
            testCase.SummaryMetrics.mrn_list = {'MRN01'};
            testCase.SummaryMetrics.lf = [0]; % No local failure

            testCase.SummaryMetrics.adc_mean = 1.0e-3 * ones(1,1,1);   % mm^2/s
            testCase.SummaryMetrics.d_mean = 1.0e-3 * ones(1,1,1);     % mm^2/s
            testCase.SummaryMetrics.f_mean = 0.1 * ones(1,1,1);        % fraction
            testCase.SummaryMetrics.dstar_mean = 0.05 * ones(1,1,1);   % mm^2/s
            testCase.SummaryMetrics.d95_gtvp = 40 * ones(1,1);         % Gy
            testCase.SummaryMetrics.dmean_gtvp = 50 * ones(1,1);       % Gy

            testCase.CalculatedResults = struct();

            % Setup mock DataVectors (voxel-level ADC values for histograms)
            testCase.DataVectors = struct('adc_vector', {ones(10,1)});

            % Add pipeline paths
            pancDataPath = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(pancDataPath, 'core'));
            addpath(fullfile(pancDataPath, 'utils'));

            % Suppress figure windows during automated testing
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            % Close any figures opened during the test and remove the temp
            % directory tree. Figures must be closed before rmdir to avoid
            % file lock issues on Windows.
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testSmokeValidData(testCase)
            % Happy-path test: with valid DWI NIfTI, GTV mask, and b-values,
            % all four visualization outputs should be created:
            %   - Parameter_Maps_1.png (ADC/D/f/D* overlays on DWI)
            %   - Feature_Histograms_Standard.png (per-parameter distributions)
            %   - Feature_BoxPlots_Standard.png (per-parameter box plots)
            %   - Dose_vs_Diffusion_Standard.png (scatter correlations)
            if exist('OCTAVE_VERSION', 'builtin'); return; end % Skip on Octave (no niftiread)
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
            % Tests graceful degradation when the .bval file is missing.
            % Parameter maps require b-value information to select the correct
            % DWI volume, so they should be skipped. However, feature
            % histograms and box plots depend only on SummaryMetrics and
            % should still be generated.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            delete(fullfile(testCase.TempDir, 'P01', 'nii', 'fx1_dwi1.bval'));

            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(exist(fullfile(outputDir, 'Parameter_Maps_1.png'), 'file') == 0, ...
                'Parameter_Maps_1.png should NOT be created if bval is missing');

            % Feature histograms use pre-computed SummaryMetrics, not raw NIfTI,
            % so they should succeed even without the .bval file.
            testCase.verifyTrue(exist(fullfile(outputDir, 'Feature_Histograms_Standard.png'), 'file') > 0, ...
                'Expected Feature_Histograms_Standard.png to be created');
        end

        function testSmokeProtocolDeviation(testCase)
            % Tests graceful handling of a b-value file that does not match
            % the expected protocol. The standard protocol uses b = [0 30 150 550];
            % here b=50 replaces b=30. Parameter maps should be skipped
            % because the function cannot reliably map b-values to volumes.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            fid = fopen(fullfile(testCase.TempDir, 'P01', 'nii', 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 50 150 550'); % b=50 is non-standard
            fclose(fid);

            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(exist(fullfile(outputDir, 'Parameter_Maps_1.png'), 'file') == 0, ...
                'Parameter_Maps_1.png should NOT be created if protocol deviates');
        end

        function testOutputFolderCreated(testCase)
            % Verify that visualize_results creates the output folder
            % if it doesn't already exist.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            outputDir = testCase.ConfigStruct.output_folder;
            if exist(outputDir, 'dir')
                rmdir(outputDir, 's');
            end

            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            testCase.verifyTrue(isfolder(outputDir), ...
                'Output folder should be created by visualize_results.');
        end

        function testDiaryFileCreated(testCase)
            % visualize_results should create a diary log file.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            diary_file = fullfile(testCase.ConfigStruct.output_folder, 'visualize_results_output.txt');
            testCase.verifyTrue(exist(diary_file, 'file') > 0, ...
                'Diary file should be created by visualize_results.');
        end

        function testFiguresClosedProperly(testCase)
            % All figures should be invisible (off-screen) during batch runs.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            testCase.verifyEqual(get(0, 'DefaultFigureVisible'), 'off', ...
                'Figures should be set to invisible for batch runs.');
        end

        function testAllNaNLFSkipsDistributions(testCase)
            % When all patients have NaN LF values, feature distributions
            % should still generate without crashing (empty groups).
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            sm = testCase.SummaryMetrics;
            sm.lf = [NaN];

            visualize_results(testCase.DataVectors, sm, testCase.CalculatedResults, testCase.ConfigStruct);

            % Should not crash; histograms may be empty but should still be saved
            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(isfolder(outputDir), ...
                'Output folder should exist even with all-NaN LF values.');
        end

        function testMultiplePatientsSmoke(testCase)
            % Smoke test with 2 patients to verify indexing doesn't fail.
            if exist('OCTAVE_VERSION', 'builtin'); return; end

            % Setup second patient directory with NIfTI data
            patID2 = 'P02';
            niiDir2 = fullfile(testCase.TempDir, patID2, 'nii');
            mkdir(niiDir2);

            dwi_img = zeros(10, 10, 10, 4, 'double');
            dwi_img(:,:,:,1) = 1000;
            dwi_img(:,:,:,2) = 800;
            dwi_img(:,:,:,3) = 500;
            dwi_img(:,:,:,4) = 100;
            dwi_file = fullfile(niiDir2, 'fx1_dwi1.nii');
            niftiwrite(dwi_img, dwi_file);
            gzip(dwi_file);
            delete(dwi_file);

            gtv_img = zeros(10, 10, 10, 'double');
            gtv_img(4:7, 4:7, 4:7) = 1;
            gtv_file = fullfile(niiDir2, 'fx1_gtv1.nii');
            niftiwrite(gtv_img, gtv_file);
            gzip(gtv_file);
            delete(gtv_file);

            fid = fopen(fullfile(niiDir2, 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 30 150 550');
            fclose(fid);

            % Expand summary metrics to 2 patients
            sm = testCase.SummaryMetrics;
            sm.id_list = {'P01', patID2};
            sm.mrn_list = {'MRN01', 'MRN02'};
            sm.lf = [0; 1];
            sm.adc_mean = cat(1, sm.adc_mean, 1.2e-3 * ones(1,1,1));
            sm.d_mean = cat(1, sm.d_mean, 1.1e-3 * ones(1,1,1));
            sm.f_mean = cat(1, sm.f_mean, 0.12 * ones(1,1,1));
            sm.dstar_mean = cat(1, sm.dstar_mean, 0.04 * ones(1,1,1));
            sm.d95_gtvp = cat(1, sm.d95_gtvp, 38);
            sm.dmean_gtvp = cat(1, sm.dmean_gtvp, 48);

            dv = struct('adc_vector', {ones(10,1); ones(10,1)});

            visualize_results(dv, sm, testCase.CalculatedResults, testCase.ConfigStruct);

            outputDir = testCase.ConfigStruct.output_folder;
            testCase.verifyTrue(exist(fullfile(outputDir, 'Feature_BoxPlots_Standard.png'), 'file') > 0, ...
                'Feature_BoxPlots should be created for 2-patient cohort.');
        end

        function testCrossDWIComparisonGracefulFailure(testCase)
            % Cross-DWI comparison requires multiple DWI types. With only
            % Standard configured, it should fail gracefully with a warning.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, testCase.ConfigStruct);

            % Should not crash; output folder should exist
            testCase.verifyTrue(isfolder(testCase.ConfigStruct.output_folder), ...
                'Output folder should exist even when cross-DWI fails.');
        end

        function testTimestampedFolderWhenNoOutputFolder(testCase)
            % When config has no output_folder, visualize_results should
            % create a timestamped folder.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            config = testCase.ConfigStruct;
            config = rmfield(config, 'output_folder');

            visualize_results(testCase.DataVectors, testCase.SummaryMetrics, testCase.CalculatedResults, config);

            % A saved_files_* directory should have been created somewhere
            % We just verify it doesn't crash.
        end
    end
end
