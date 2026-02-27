classdef test_visualize_refactor < matlab.unittest.TestCase
    properties
        TempDir
        ConfigStruct
        SummaryMetrics
        DataVectors
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

            % Create dummy .nii.gz files (empty)
            fclose(fopen(fullfile(niiDir, 'fx1_dwi1.nii.gz'), 'w'));
            fclose(fopen(fullfile(niiDir, 'fx1_gtv1.nii.gz'), 'w'));

            % Create dummy .bval
            fid = fopen(fullfile(niiDir, 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 30 150 550');
            fclose(fid);

            % Setup inputs
            % IMPORTANT: We deliberately do NOT add a trailing separator to test fullfile robustness
            testCase.ConfigStruct.dataloc = testCase.TempDir;

            testCase.SummaryMetrics.id_list = {patID};
            testCase.SummaryMetrics.mrn_list = {'MRN01'};
            testCase.SummaryMetrics.lf = [0];
            % Initialize summary arrays with NaNs/Ones
            testCase.SummaryMetrics.adc_mean = ones(1,6,3);
            testCase.SummaryMetrics.d_mean = ones(1,6,3);
            testCase.SummaryMetrics.f_mean = ones(1,6,3);
            testCase.SummaryMetrics.dstar_mean = ones(1,6,3);
            testCase.SummaryMetrics.d95_gtvp = ones(1,6);
            testCase.SummaryMetrics.dmean_gtvp = ones(1,6);

            % Setup DataVectors with a non-empty struct to pass the "isempty" check
            testCase.DataVectors = struct('adc_vector', {ones(10,1)});

            % Add core/utils to path
            addpath(fullfile(pwd, 'core'));
            addpath(fullfile(pwd, 'utils'));
            addpath(fullfile(pwd, 'dependencies'));
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
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

            % Check the diary output
            outputDir = fullfile(pwd, 'saved_figures');
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
