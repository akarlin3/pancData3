classdef test_init_scan_structs < matlab.unittest.TestCase
    % TEST_INIT_SCAN_STRUCTS Unit tests for the scan struct initializer.
    %
    % init_scan_structs creates two struct arrays (GTVp and GTVn) with
    % dimensions [n_fx x n_rp] (fractions x repeat scans). Each element
    % contains fields for diffusion parameter vectors (Standard, dnCNN,
    % IVIMnet), dose data, patient metadata, and voxel geometry. All fields
    % are initialized to empty ([]).
    %
    % Validates:
    %   - Both GTVp and GTVn have exactly the expected 23 field names
    %   - Struct array dimensions match the requested n_fx x n_rp
    %   - Every field in every element is initialized to empty []
    %   - Single-element (1x1) case works correctly
    %   - Field count is exactly 23

    methods(TestMethodSetup)
        function addPaths(testCase)
            % Add the utils directory containing init_scan_structs.
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testExpectedFieldNames(testCase)
            % Verify that both GTVp and GTVn structs contain exactly the
            % expected set of 23 fields covering: Standard DWI vectors
            % (adc, d, f, dstar), dose metrics (dose_vector, dvh, d95, v50gy),
            % dnCNN-denoised vectors, IVIMnet vectors, patient metadata
            % (ID, MRN, LF, Immuno), and scan geometry (Fraction,
            % Repeatability_index, vox_vol, vox_dims).
            [gtvp, gtvn] = init_scan_structs(1, 1);

            expected_fields = {'adc_vector', 'd_vector', 'f_vector', 'dstar_vector', ...
                'dose_vector', 'dvh', 'd95', 'v50gy', ...
                'd_vector_dncnn', 'f_vector_dncnn', 'dstar_vector_dncnn', ...
                'adc_vector_dncnn', ...
                'd_vector_ivimnet', 'f_vector_ivimnet', 'dstar_vector_ivimnet', ...
                'ID', 'MRN', 'LF', 'Immuno', ...
                'Fraction', 'Repeatability_index', 'vox_vol', 'vox_dims'};

            actual_fields = fieldnames(gtvp);
            testCase.verifyEqual(sort(actual_fields), sort(expected_fields'), ...
                'GTVp struct should have exactly the expected fields.');

            actual_fields_n = fieldnames(gtvn);
            testCase.verifyEqual(sort(actual_fields_n), sort(expected_fields'), ...
                'GTVn struct should have exactly the expected fields.');
        end

        function testDimensionsMatchRequest(testCase)
            % Verify that a 5x3 request produces struct arrays of size [5, 3],
            % representing 5 fractions and 3 repeat scans per fraction.
            n_fx = 5;
            n_rp = 3;
            [gtvp, gtvn] = init_scan_structs(n_fx, n_rp);

            testCase.verifyEqual(size(gtvp), [n_fx, n_rp], ...
                'GTVp struct array should be n_fx x n_rp.');
            testCase.verifyEqual(size(gtvn), [n_fx, n_rp], ...
                'GTVn struct array should be n_fx x n_rp.');
        end

        function testAllFieldsInitToEmpty(testCase)
            % Every field in every element of a 3x2 struct array should be
            % empty ([]), ensuring no residual data from previous pipeline runs.
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
            % The minimal case: 1 fraction, 1 repeat scan. Verify that the
            % output is a 1x1 struct (not a cell or other type).
            [gtvp, gtvn] = init_scan_structs(1, 1);

            testCase.verifyEqual(size(gtvp), [1, 1]);
            testCase.verifyEqual(size(gtvn), [1, 1]);
            testCase.verifyTrue(isstruct(gtvp));
            testCase.verifyTrue(isstruct(gtvn));
        end

        function testFieldCount(testCase)
            % Guard against accidental field additions or removals.
            % The struct should have exactly 23 fields.
            [gtvp, ~] = init_scan_structs(1, 1);

            n_fields = numel(fieldnames(gtvp));
            testCase.verifyEqual(n_fields, 23, ...
                'Struct should have exactly 23 fields.');
        end

    end
end
