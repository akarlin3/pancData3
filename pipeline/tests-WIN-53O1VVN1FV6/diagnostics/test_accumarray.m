classdef test_accumarray < matlab.unittest.TestCase
% TEST_ACCUMARRAY  Spot-check that MATLAB's accumarray correctly aggregates
%   binary outcomes to the patient level using @any.
%
%   This validates the core operation used in make_grouped_folds to reduce
%   multi-row-per-patient outcomes to a single binary flag per patient
%   (1 if any row is positive, 0 otherwise).

    methods (Test)
        function testAccumarrayPatientAggregation(testCase)
        %TESTACCUMARRAYPATIENTAGGREGATION Verify accumarray with @any
        %   correctly produces per-group binary outcomes.
        %
        %   Groups: Patient 1 (rows 1-2, all 0) -> 0
        %           Patient 2 (rows 3-4, one 1) -> 1
        %           Patient 3 (rows 5-6, all 1) -> 1
            ic = [1; 1; 2; 2; 3; 3];       % Group index (patient membership)
            y = [0; 0; 0; 1; 1; 1];        % Row-level binary outcomes
            n_unique = 3;
            pt_y = double(accumarray(ic, double(y > 0), [n_unique, 1], @any));
            testCase.verifyEqual(pt_y, [0; 1; 1]);
        end
    end
end
