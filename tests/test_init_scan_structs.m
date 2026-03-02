classdef test_init_scan_structs < matlab.unittest.TestCase
    % TEST_INIT_SCAN_STRUCTS Unit tests for the scan struct initializer.
    %
    % Validates:
    %   - Output has expected fields
    %   - Dimensions match requested n_fx x n_rp
    %   - All fields initialize to empty []

    methods(TestMethodSetup)
        function addPaths(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testExpectedFieldNames(testCase)
            [gtvp, gtvn] = init_scan_structs(1, 1);

            expected_fields = {'adc_vector', 'd_vector', 'f_vector', 'dstar_vector', ...
                'dose_vector', 'dvh', 'd95', 'v50gy', ...
                'd_vector_dncnn', 'f_vector_dncnn', 'dstar_vector_dncnn', ...
                'adc_vector_dncnn', ...
                'd_vector_ivimnet', 'f_vector_ivimnet', 'dstar_vector_ivimnet', ...
                'ID', 'MRN', 'LF', 'Immuno', ...
                'Fraction', 'Repeatability_index', 'vox_vol'};

            actual_fields = fieldnames(gtvp);
            testCase.verifyEqual(sort(actual_fields), sort(expected_fields'), ...
                'GTVp struct should have exactly the expected fields.');

            actual_fields_n = fieldnames(gtvn);
            testCase.verifyEqual(sort(actual_fields_n), sort(expected_fields'), ...
                'GTVn struct should have exactly the expected fields.');
        end

        function testDimensionsMatchRequest(testCase)
            n_fx = 5;
            n_rp = 3;
            [gtvp, gtvn] = init_scan_structs(n_fx, n_rp);

            testCase.verifyEqual(size(gtvp), [n_fx, n_rp], ...
                'GTVp struct array should be n_fx x n_rp.');
            testCase.verifyEqual(size(gtvn), [n_fx, n_rp], ...
                'GTVn struct array should be n_fx x n_rp.');
        end

        function testAllFieldsInitToEmpty(testCase)
            [gtvp, gtvn] = init_scan_structs(3, 2);

            fields = fieldnames(gtvp);
            for fi = 1:numel(fields)
                for i = 1:3
                    for j = 1:2
                        testCase.verifyEmpty(gtvp(i, j).(fields{fi}), ...
                            sprintf('GTVp(%d,%d).%s should be empty.', i, j, fields{fi}));
                        testCase.verifyEmpty(gtvn(i, j).(fields{fi}), ...
                            sprintf('GTVn(%d,%d).%s should be empty.', i, j, fields{fi}));
                    end
                end
            end
        end

        function testSingleElement(testCase)
            [gtvp, gtvn] = init_scan_structs(1, 1);

            testCase.verifyEqual(size(gtvp), [1, 1]);
            testCase.verifyEqual(size(gtvn), [1, 1]);
            testCase.verifyTrue(isstruct(gtvp));
            testCase.verifyTrue(isstruct(gtvn));
        end

        function testFieldCount(testCase)
            [gtvp, ~] = init_scan_structs(1, 1);

            n_fields = numel(fieldnames(gtvp));
            testCase.verifyEqual(n_fields, 22, ...
                'Struct should have exactly 22 fields.');
        end

    end
end
