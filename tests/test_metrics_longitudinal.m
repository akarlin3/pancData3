classdef test_metrics_longitudinal < matlab.unittest.TestCase
    % TEST_METRICS_LONGITUDINAL Smoke tests for the metrics_longitudinal function.
    %
    % Verifies that:
    %   - The longitudinal figure is created and saved to the output folder
    %   - The function handles single-patient datasets without error
    %   - All-NaN data does not crash the function
    %   - nTp values from 2 to 6 are all handled correctly

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testOutputFigureCreated(testCase)
            % Standard call with 10 patients and 4 timepoints.
            % Verifies the expected PNG is written to output_folder.
            rng(1);
            n   = 10;
            nTp = 4;
            dtype_label = 'Standard';

            ADC_abs   = rand(n, nTp) * 2e-3;
            D_abs     = rand(n, nTp) * 1e-3;
            f_abs     = rand(n, nTp) * 0.3;
            Dstar_abs = rand(n, nTp) * 0.05;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_pct     = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_pct, Dstar_pct, nTp, dtype_label, testCase.TempDir);

            expected_file = fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '.png']);
            testCase.verifyTrue(exist(expected_file, 'file') > 0, ...
                'Longitudinal figure PNG should be saved to output_folder.');
        end

        function testSinglePatientNoError(testCase)
            % Edge case: only one patient. Mean = the patient's value; SEM = 0.
            nTp = 4;
            ADC_abs   = rand(1, nTp) * 2e-3;
            D_abs     = rand(1, nTp) * 1e-3;
            f_abs     = rand(1, nTp) * 0.3;
            Dstar_abs = rand(1, nTp) * 0.05;
            ADC_pct   = rand(1, nTp) * 10;
            D_pct     = rand(1, nTp) * 10;
            f_pct     = rand(1, nTp) * 0.05;
            Dstar_pct = rand(1, nTp) * 10;

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_pct, Dstar_pct, nTp, 'SinglePt', testCase.TempDir);
            testCase.verifyTrue(true, 'Single-patient case should not throw.');
        end

        function testAllNaNDataNoError(testCase)
            % All-NaN input should not crash — NaN-aware mean/std handle it.
            nTp = 4;
            n   = 8;
            ADC_abs   = nan(n, nTp);
            D_abs     = nan(n, nTp);
            f_abs     = nan(n, nTp);
            Dstar_abs = nan(n, nTp);
            ADC_pct   = nan(n, nTp);
            D_pct     = nan(n, nTp);
            f_pct     = nan(n, nTp);
            Dstar_pct = nan(n, nTp);

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_pct, Dstar_pct, nTp, 'AllNaN', testCase.TempDir);
            testCase.verifyTrue(true, 'All-NaN data should not throw.');
        end

        function testSixTimepointsNoError(testCase)
            % nTp = 6 (maximum in this study) should generate a valid figure.
            rng(5);
            n   = 15;
            nTp = 6;

            ADC_abs   = rand(n, nTp) * 2e-3;
            D_abs     = rand(n, nTp) * 1e-3;
            f_abs     = rand(n, nTp) * 0.3;
            Dstar_abs = rand(n, nTp) * 0.05;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_pct     = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_pct, Dstar_pct, nTp, 'SixTp', testCase.TempDir);

            expected_file = fullfile(testCase.TempDir, 'Longitudinal_Mean_Metrics_SixTp.png');
            testCase.verifyTrue(exist(expected_file, 'file') > 0, ...
                'Six-timepoint figure should be saved.');
        end

        function testDtypeLabelAppearsInFilename(testCase)
            % The dtype_label is embedded in the saved filename.
            nTp = 3;
            n   = 5;
            lbl = 'dnCNN';

            metrics_longitudinal( ...
                rand(n, nTp), rand(n, nTp), rand(n, nTp), rand(n, nTp), ...
                rand(n, nTp), rand(n, nTp), rand(n, nTp), rand(n, nTp), ...
                nTp, lbl, testCase.TempDir);

            expected_file = fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' lbl '.png']);
            testCase.verifyTrue(exist(expected_file, 'file') > 0, ...
                'Filename should include the dtype_label.');
        end

    end
end
