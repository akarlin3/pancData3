classdef test_accumarray < matlab.unittest.TestCase
    % Spot-check of MATLAB's accumarray for patient-level aggregation

    methods (Test)
        function testAccumarrayPatientAggregation(testCase)
            ic = [1; 1; 2; 2; 3; 3];
            y = [0; 0; 0; 1; 1; 1];
            n_unique = 3;
            pt_y = double(accumarray(ic, double(y > 0), [n_unique, 1], @any));
            testCase.verifyEqual(pt_y, [0; 1; 1]);
        end
    end
end
