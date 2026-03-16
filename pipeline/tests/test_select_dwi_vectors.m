classdef test_select_dwi_vectors < matlab.unittest.TestCase
% TEST_SELECT_DWI_VECTORS — Unit tests for select_dwi_vectors.m
%
% Validates DWI vector extraction for all three processing types:
%   - Standard (type 1): raw ADC, D, f, D*
%   - dnCNN (type 2): denoised variants (_dncnn suffix)
%   - IVIMnet (type 3): standard ADC + IVIMnet D/f/D* (_ivimnet suffix)

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
        end
    end

    methods (Test)
        function test_standard_type(testCase)
            s = testCase.buildStruct();
            [adc, d, f, dstar] = select_dwi_vectors(s, 1, 1, 1, 1);
            testCase.verifyEqual(adc, [1; 2; 3]);
            testCase.verifyEqual(d, [4; 5; 6]);
            testCase.verifyEqual(f, [0.1; 0.2; 0.3]);
            testCase.verifyEqual(dstar, [10; 20; 30]);
        end

        function test_dncnn_type(testCase)
            s = testCase.buildStruct();
            [adc, d, f, dstar] = select_dwi_vectors(s, 1, 1, 1, 2);
            testCase.verifyEqual(adc, [11; 12; 13]);
            testCase.verifyEqual(d, [14; 15; 16]);
            testCase.verifyEqual(f, [0.11; 0.12; 0.13]);
            testCase.verifyEqual(dstar, [110; 120; 130]);
        end

        function test_ivimnet_type(testCase)
            s = testCase.buildStruct();
            [adc, d, f, dstar] = select_dwi_vectors(s, 1, 1, 1, 3);
            % IVIMnet reuses standard ADC
            testCase.verifyEqual(adc, [1; 2; 3], ...
                'IVIMnet should reuse standard ADC.');
            testCase.verifyEqual(d, [24; 25; 26]);
            testCase.verifyEqual(f, [0.21; 0.22; 0.23]);
            testCase.verifyEqual(dstar, [210; 220; 230]);
        end

        function test_multi_index(testCase)
            % Test accessing different patient/timepoint/repeat indices
            s(1,1,1) = testCase.buildEntry(1);
            s(2,1,1) = testCase.buildEntry(2);
            s(1,2,1) = testCase.buildEntry(3);

            [adc1, ~, ~, ~] = select_dwi_vectors(s, 1, 1, 1, 1);
            [adc2, ~, ~, ~] = select_dwi_vectors(s, 2, 1, 1, 1);
            [adc3, ~, ~, ~] = select_dwi_vectors(s, 1, 2, 1, 1);

            testCase.verifyNotEqual(adc1, adc2);
            testCase.verifyNotEqual(adc1, adc3);
        end

        function test_output_types(testCase)
            s = testCase.buildStruct();
            for dtype = 1:3
                [adc, d, f, dstar] = select_dwi_vectors(s, 1, 1, 1, dtype);
                testCase.verifyTrue(isnumeric(adc));
                testCase.verifyTrue(isnumeric(d));
                testCase.verifyTrue(isnumeric(f));
                testCase.verifyTrue(isnumeric(dstar));
            end
        end
    end

    methods (Static, Access = private)
        function s = buildStruct()
            s(1,1,1) = test_select_dwi_vectors.buildEntry(1);
        end

        function e = buildEntry(seed)
            b = seed * 10;
            e.adc_vector = [1; 2; 3] * seed;
            e.d_vector = [4; 5; 6] * seed;
            e.f_vector = [0.1; 0.2; 0.3] * seed;
            e.dstar_vector = [10; 20; 30] * seed;
            e.adc_vector_dncnn = [11; 12; 13] * seed;
            e.d_vector_dncnn = [14; 15; 16] * seed;
            e.f_vector_dncnn = [0.11; 0.12; 0.13] * seed;
            e.dstar_vector_dncnn = [110; 120; 130] * seed;
            e.d_vector_ivimnet = [24; 25; 26] * seed;
            e.f_vector_ivimnet = [0.21; 0.22; 0.23] * seed;
            e.dstar_vector_ivimnet = [210; 220; 230] * seed;
        end
    end
end
