classdef test_masking_bug < matlab.unittest.TestCase
    methods(Test)
        function testWCV_Masking_BugReproduction(testCase)
            % Reproduces the logic bug found in metrics.m
            % adc_wCV(n_rpt<2) = nan;

            n_patients = 2;
            n_pipelines = 2;
            % 2 patients, 2 pipelines
            adc_wCV = ones(n_patients, n_pipelines);

            % Patient 2 has < 2 repeats
            n_rpt = [2; 1];

            % Apply the buggy logic found in the codebase
            % This uses linear indexing on the matrix with a column vector mask
            adc_wCV_buggy = adc_wCV;
            adc_wCV_buggy(n_rpt < 2) = nan;

            % Expectation: Only the first element corresponding to linear index of n_rpt<2 is masked
            % Linear index of (2,1) is 2.
            % Linear index of (2,2) is 4.
            % The mask (0; 1) only hits index 2.

            testCase.verifyTrue(isnan(adc_wCV_buggy(2,1)), 'Patient 2 Col 1 should be NaN');
            testCase.verifyFalse(isnan(adc_wCV_buggy(2,2)), 'Patient 2 Col 2 should NOT be NaN (Bug detected)');
        end

        function testWCV_Masking_Fix(testCase)
            % Verifies the proposed fix
            % adc_wCV(n_rpt<2, :) = nan;

            n_patients = 2;
            n_pipelines = 2;
            adc_wCV = ones(n_patients, n_pipelines);
            n_rpt = [2; 1];

            % Apply fixed logic
            adc_wCV_fixed = adc_wCV;
            adc_wCV_fixed(n_rpt < 2, :) = nan;

            testCase.verifyTrue(isnan(adc_wCV_fixed(2,1)), 'Patient 2 Col 1 should be NaN');
            testCase.verifyTrue(isnan(adc_wCV_fixed(2,2)), 'Patient 2 Col 2 should be NaN (Fix verified)');
        end

        function testIVIM_Masking_BugReproduction(testCase)
            % Reproduces the logic bug/redundancy for IVIM metrics
            % d_wCV(repmat(n_rpt,[1,3])<2) = nan;

            n_patients = 2;
            n_pipelines = 3; % IVIM has 3 pipelines
            d_wCV = ones(n_patients, n_pipelines);
            n_rpt = [2; 1];

            % Current code uses repmat, which works but is verbose
            d_wCV_current = d_wCV;
            d_wCV_current(repmat(n_rpt,[1,3])<2) = nan;

            testCase.verifyTrue(all(isnan(d_wCV_current(2,:))), 'Patient 2 should be all NaN with repmat');

            % Proposed fix: clean column indexing
            d_wCV_fixed = d_wCV;
            d_wCV_fixed(n_rpt < 2, :) = nan;

            testCase.verifyTrue(all(isnan(d_wCV_fixed(2,:))), 'Patient 2 should be all NaN with fixed logic');
            testCase.verifyEqual(d_wCV_fixed, d_wCV_current, 'Fixed logic should produce identical result to repmat');
        end
    end
end
