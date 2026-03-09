classdef test_visualize_refactor < matlab.unittest.TestCase
    % TEST_VISUALIZE_REFACTOR Unit tests for visualization path construction.
    %
    % This test file validates that visualize_results.m constructs file paths
    % correctly using fullfile() rather than string concatenation. A common bug
    % pattern is [dataloc id '/'] which fails when dataloc lacks a trailing
    % separator, producing malformed paths like ".../TempDirP01/" instead of
    % ".../TempDir/P01/". The test creates a minimal patient directory with
    % dummy NIfTI files and verifies that the visualization module can locate
    % them, confirming correct path construction.
    properties
        TempDir
        ConfigStruct
        SummaryMetrics
        DataVectors
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Create a temporary directory tree mimicking the patient data layout
            % expected by visualize_results: dataloc/PatientID/nii/*.nii.gz
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Setup patient structure
            patID = 'P01';
            niiDir = fullfile(testCase.TempDir, patID, 'nii');
            mkdir(niiDir);

            % Create dummy .nii.gz files (empty)
            fclose(fopen(fullfile(niiDir, 'fx1_dwi1.nii.gz'), 'w'));
            fclose(fopen(fullfile(niiDir, 'fx1_gtv1.nii.gz'), 'w'));

            % Create a dummy .bval file with 4 b-values (required by visualize_results
            % to determine the number of diffusion weightings)
            fid = fopen(fullfile(niiDir, 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 30 150 550');
            fclose(fid);

            % Setup config inputs.
            % IMPORTANT: We deliberately do NOT add a trailing separator to test fullfile robustness.
            % This is the crux of the test: if the code uses string concatenation instead of
            % fullfile(), the path will be malformed.
            testCase.ConfigStruct.dataloc = testCase.TempDir;
            testCase.ConfigStruct.output_folder = fullfile(testCase.TempDir, 'output');
            testCase.ConfigStruct.dwi_types_to_run = 1;

            testCase.SummaryMetrics.id_list = {patID};
            testCase.SummaryMetrics.mrn_list = {'MRN01'};
            testCase.SummaryMetrics.lf = [0];
            % Initialize summary arrays with ones: dimensions are [nPatients, nTimepoints, nDWITypes]
            % These placeholder values allow visualize_results to proceed past data checks.
            testCase.SummaryMetrics.adc_mean = ones(1,6,3);
            testCase.SummaryMetrics.d_mean = ones(1,6,3);
            testCase.SummaryMetrics.f_mean = ones(1,6,3);
            testCase.SummaryMetrics.dstar_mean = ones(1,6,3);
            testCase.SummaryMetrics.d95_gtvp = ones(1,6);
            testCase.SummaryMetrics.dmean_gtvp = ones(1,6);

            % Setup DataVectors with a non-empty struct so visualize_results does not
            % short-circuit on the isempty() guard at the top of the function.
            testCase.DataVectors = struct('adc_vector', {ones(10,1)});

            % Add core/utils/dependencies to the MATLAB path so visualize_results
            % and its helpers are accessible.
            baseDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            % Close any diary opened by visualize_results before removing the temp
            % directory; otherwise the diary file lock (especially on Windows)
            % prevents rmdir from succeeding.
            diary off;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testVisualizeRuns(testCase)
            % This test verifies that visualize_results constructs paths correctly.
            %
            % If paths are constructed using [dataloc id '/'], and dataloc lacks a trailing slash,
            % the path becomes ".../TempDirP01/..." which does not exist.
            % The code will then check exist(), find nothing, and 'continue' (skip).
            % The output diary will say "Plotted 0 patients".
            %
            % If paths are constructed using fullfile(dataloc, id), the path becomes
            % ".../TempDir/P01/..." which DOES exist.
            % The code will check exist(), find the file, and proceed to niftiinfo().
            % Since we created empty dummy files, niftiinfo/niftiread will fail.
            %
            % Therefore:
            % - Success = Caught NIfTI error (means path was found).
            % - Failure = "Plotted 0 patients" (means path was not found).

            try
                visualize_results(testCase.DataVectors, testCase.SummaryMetrics, struct(), testCase.ConfigStruct);
            catch ME
                % We expect an error because the .nii.gz files are empty (0 bytes)
                % This confirms the code successfully found the files and tried to read them.
                if contains(ME.message, 'File') || contains(ME.identifier, 'images') || contains(ME.identifier, 'nifti')
                     % Success!
                     return;
                end
                % If it's some other error, rethrow it
                rethrow(ME);
            end

            % If we reach here, it means no error was thrown.
            % This implies either:
            % 1. It successfully read empty files (unlikely for niftiread)
            % 2. It skipped the patient (path check failed)

            % Check the diary output in the configured output folder
            outputDir = testCase.ConfigStruct.output_folder;
            diaryFile = fullfile(outputDir, 'visualize_results_output.txt');

            if exist(diaryFile, 'file')
                content = fileread(diaryFile);
                if contains(content, 'Plotted 0 patients') || contains(content, 'No patients with complete NIfTI data found')
                     error('TestFailed:PathConstruction', ...
                         'visualize_results failed to find the patient files. Path construction likely failed (missing separator?).');
                end
            end
        end
    end
end
