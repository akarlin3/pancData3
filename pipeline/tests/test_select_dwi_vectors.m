classdef test_select_dwi_vectors < matlab.unittest.TestCase
% TEST_SELECT_DWI_VECTORS — Unit tests for select_dwi_vectors.m
%
% Validates DWI vector extraction including:
%   - Standard (type 1): uses base field names
%   - dnCNN (type 2): uses _dncnn suffixed fields
%   - IVIMnet (type 3): uses standard ADC but _ivimnet for D/f/D*
%   - Correct indexing into 3D struct array

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_standard_type1(testCase)
            dv = testCase.makeStruct();

            [adc, d, f, ds] = select_dwi_vectors(dv, 1, 1, 1, 1);

            testCase.verifyEqual(adc, [1; 2; 3], 'Standard ADC should use adc_vector.');
            testCase.verifyEqual(d, [4; 5; 6], 'Standard D should use d_vector.');
            testCase.verifyEqual(f, [7; 8; 9], 'Standard f should use f_vector.');
            testCase.verifyEqual(ds, [10; 11; 12], 'Standard D* should use dstar_vector.');
        end

        function test_dncnn_type2(testCase)
            dv = testCase.makeStruct();

            [adc, d, f, ds] = select_dwi_vectors(dv, 1, 1, 1, 2);

            testCase.verifyEqual(adc, [21; 22; 23], 'dnCNN ADC should use adc_vector_dncnn.');
            testCase.verifyEqual(d, [24; 25; 26], 'dnCNN D should use d_vector_dncnn.');
            testCase.verifyEqual(f, [27; 28; 29], 'dnCNN f should use f_vector_dncnn.');
            testCase.verifyEqual(ds, [30; 31; 32], 'dnCNN D* should use dstar_vector_dncnn.');
        end

        function test_ivimnet_type3_uses_standard_adc(testCase)
            dv = testCase.makeStruct();

            [adc, d, f, ds] = select_dwi_vectors(dv, 1, 1, 1, 3);

            testCase.verifyEqual(adc, [1; 2; 3], ...
                'IVIMnet ADC should reuse standard adc_vector (mono-exponential fit).');
            testCase.verifyEqual(d, [34; 35; 36], 'IVIMnet D should use d_vector_ivimnet.');
            testCase.verifyEqual(f, [37; 38; 39], 'IVIMnet f should use f_vector_ivimnet.');
            testCase.verifyEqual(ds, [40; 41; 42], 'IVIMnet D* should use dstar_vector_ivimnet.');
        end

        function test_indexing_into_3d_struct(testCase)
            % Build a 2x2x2 struct array and verify correct indexing
            dv = testCase.makeStruct();
            dv(2,2,2) = dv(1,1,1);  % preallocate
            dv(2,1,1).adc_vector = [100; 200];
            dv(2,1,1).d_vector = [101; 201];
            dv(2,1,1).f_vector = [102; 202];
            dv(2,1,1).dstar_vector = [103; 203];

            [adc, d, f, ds] = select_dwi_vectors(dv, 2, 1, 1, 1);

            testCase.verifyEqual(adc, [100; 200]);
            testCase.verifyEqual(d, [101; 201]);
        end

        function test_empty_vectors(testCase)
            % Struct with empty vectors should return empty
            dv(1,1,1).adc_vector = [];
            dv(1,1,1).d_vector = [];
            dv(1,1,1).f_vector = [];
            dv(1,1,1).dstar_vector = [];
            dv(1,1,1).adc_vector_dncnn = [];
            dv(1,1,1).d_vector_dncnn = [];
            dv(1,1,1).f_vector_dncnn = [];
            dv(1,1,1).dstar_vector_dncnn = [];
            dv(1,1,1).d_vector_ivimnet = [];
            dv(1,1,1).f_vector_ivimnet = [];
            dv(1,1,1).dstar_vector_ivimnet = [];

            [adc, d, f, ds] = select_dwi_vectors(dv, 1, 1, 1, 1);

            testCase.verifyTrue(isempty(adc));
            testCase.verifyTrue(isempty(d));
            testCase.verifyTrue(isempty(f));
            testCase.verifyTrue(isempty(ds));
        end
    end

    methods (Static, Access = private)
        function dv = makeStruct()
            % Create a 1x1x1 struct with all DWI vector fields populated
            dv(1,1,1).adc_vector = [1; 2; 3];
            dv(1,1,1).d_vector = [4; 5; 6];
            dv(1,1,1).f_vector = [7; 8; 9];
            dv(1,1,1).dstar_vector = [10; 11; 12];
            dv(1,1,1).adc_vector_dncnn = [21; 22; 23];
            dv(1,1,1).d_vector_dncnn = [24; 25; 26];
            dv(1,1,1).f_vector_dncnn = [27; 28; 29];
            dv(1,1,1).dstar_vector_dncnn = [30; 31; 32];
            dv(1,1,1).d_vector_ivimnet = [34; 35; 36];
            dv(1,1,1).f_vector_ivimnet = [37; 38; 39];
            dv(1,1,1).dstar_vector_ivimnet = [40; 41; 42];
        end
    end
end
