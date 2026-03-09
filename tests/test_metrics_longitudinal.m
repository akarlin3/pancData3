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
            % Verifies the standard case: 10 patients, 4 timepoints.
            % The function should save a PNG file named
            % Longitudinal_Mean_Metrics_{dtype_label}.png to output_folder.
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
            f_delta     = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, dtype_label, testCase.TempDir);

            expected_file = fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '.png']);
            testCase.verifyTrue(exist(expected_file, 'file') > 0, ...
                'Longitudinal figure PNG should be saved to output_folder.');
        end

        function testSinglePatientNoError(testCase)
            % Edge case: single patient. nanmean equals the patient's own
            % values and SEM is NaN (undefined for n=1). The function must
            % handle this without dividing by zero or erroring on NaN SEM.
            nTp = 4;
            ADC_abs   = rand(1, nTp) * 2e-3;
            D_abs     = rand(1, nTp) * 1e-3;
            f_abs     = rand(1, nTp) * 0.3;
            Dstar_abs = rand(1, nTp) * 0.05;
            ADC_pct   = rand(1, nTp) * 10;
            D_pct     = rand(1, nTp) * 10;
            f_delta     = rand(1, nTp) * 0.05;
            Dstar_pct = rand(1, nTp) * 10;

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, 'SinglePt', testCase.TempDir);
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
            f_delta     = nan(n, nTp);
            Dstar_pct = nan(n, nTp);

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, 'AllNaN', testCase.TempDir);
            testCase.verifyTrue(true, 'All-NaN data should not throw.');
        end

        function testSixTimepointsNoError(testCase)
            % Verifies that nTp=6 (the maximum timepoint count in the
            % pancreatic DWI study protocol) produces a valid figure.
            % Tests that x-axis label generation handles all 6 fractions.
            rng(5);
            n   = 15;
            nTp = 6;

            ADC_abs   = rand(n, nTp) * 2e-3;
            D_abs     = rand(n, nTp) * 1e-3;
            f_abs     = rand(n, nTp) * 0.3;
            Dstar_abs = rand(n, nTp) * 0.05;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_delta     = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, 'SixTp', testCase.TempDir);

            expected_file = fullfile(testCase.TempDir, 'Longitudinal_Mean_Metrics_SixTp.png');
            testCase.verifyTrue(exist(expected_file, 'file') > 0, ...
                'Six-timepoint figure should be saved.');
        end

        function testDtypeLabelAppearsInFilename(testCase)
            % Verifies that the dtype_label string ('dnCNN' here) is
            % correctly embedded in the output PNG filename, enabling
            % automated identification of DWI type in batch output folders.
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

        function testOutcomeStratifiedFigures(testCase)
            % Verifies that when the optional m_lf argument is provided,
            % separate per-outcome figures are created for each group
            % present in the data: LC (lf=0), LF (lf=1), and Competing
            % Risk (lf=2). Each group gets its own PNG file.
            rng(2);
            n   = 12;
            nTp = 4;
            dtype_label = 'Standard';

            ADC_abs   = rand(n, nTp) * 2e-3;
            D_abs     = rand(n, nTp) * 1e-3;
            f_abs     = rand(n, nTp) * 0.3;
            Dstar_abs = rand(n, nTp) * 0.05;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_delta   = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            % 6 LC, 4 LF, 2 Competing Risk
            m_lf = [0; 0; 0; 0; 0; 0; 1; 1; 1; 1; 2; 2];

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, dtype_label, ...
                testCase.TempDir, m_lf);

            % Verify per-outcome PNGs
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_LC.png']), 'file') > 0, ...
                'LC figure should be saved.');
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_LF.png']), 'file') > 0, ...
                'LF figure should be saved.');
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_CompetingRisk.png']), 'file') > 0, ...
                'Competing risk figure should be saved.');
        end

        function testCombinedStratifiedFigure(testCase)
            % Verifies that a combined overlay figure (all outcome groups
            % on the same axes for comparison) is saved as _ByOutcome.png
            % in addition to the individual per-outcome figures.
            rng(3);
            n   = 10;
            nTp = 4;
            dtype_label = 'Standard';

            ADC_abs   = rand(n, nTp) * 2e-3;
            D_abs     = rand(n, nTp) * 1e-3;
            f_abs     = rand(n, nTp) * 0.3;
            Dstar_abs = rand(n, nTp) * 0.05;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_delta   = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            m_lf = [0; 0; 0; 0; 0; 1; 1; 1; 1; 1];

            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, dtype_label, ...
                testCase.TempDir, m_lf);

            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_ByOutcome.png']), 'file') > 0, ...
                'Combined by-outcome figure should be saved.');
        end

        function testBackwardCompatibilityNoOutcome(testCase)
            % Verifies backward compatibility: calling with 11 arguments
            % (omitting m_lf) should produce only the all-patients figure
            % and NO outcome-stratified figures. This ensures older code
            % that doesn't pass m_lf continues to work.
            rng(4);
            n   = 8;
            nTp = 4;
            dtype_label = 'Standard';

            ADC_abs   = rand(n, nTp) * 2e-3;
            D_abs     = rand(n, nTp) * 1e-3;
            f_abs     = rand(n, nTp) * 0.3;
            Dstar_abs = rand(n, nTp) * 0.05;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_delta   = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            % Call with 11 args only (no m_lf)
            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ...
                ADC_pct, D_pct, f_delta, Dstar_pct, nTp, dtype_label, ...
                testCase.TempDir);

            % All-patients figure should exist
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '.png']), 'file') > 0, ...
                'All-patients figure should be saved without m_lf.');
            % No outcome-stratified figures should exist
            testCase.verifyFalse(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_ByOutcome.png']), 'file') > 0, ...
                'No outcome figures should be generated without m_lf.');
        end

        function testMissingOutcomeGroupSkipped(testCase)
            % Verifies that outcome groups not present in the data are
            % gracefully skipped. With m_lf containing only 0 (LC) and 1
            % (LF) but no 2 (competing risk), only LC and LF figures
            % should be created — no CompetingRisk PNG should exist.
            rng(5);
            n   = 6;
            nTp = 3;
            dtype_label = 'Standard';

            m_lf = [0; 0; 0; 1; 1; 1];  % No competing risk

            metrics_longitudinal( ...
                rand(n, nTp), rand(n, nTp), rand(n, nTp), rand(n, nTp), ...
                rand(n, nTp), rand(n, nTp), rand(n, nTp), rand(n, nTp), ...
                nTp, dtype_label, testCase.TempDir, m_lf);

            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_LC.png']), 'file') > 0, ...
                'LC figure should exist.');
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_LF.png']), 'file') > 0, ...
                'LF figure should exist.');
            testCase.verifyFalse(exist(fullfile(testCase.TempDir, ...
                ['Longitudinal_Mean_Metrics_' dtype_label '_CompetingRisk.png']), 'file') > 0, ...
                'Competing risk figure should NOT exist when no CR patients.');
        end

    end
end
