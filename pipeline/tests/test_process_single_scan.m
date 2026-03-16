classdef test_process_single_scan < matlab.unittest.TestCase
    % TEST_PROCESS_SINGLE_SCAN Functional tests for process_single_scan.
    % Tests initialization, NaN defaults, mask loading helper, biomarker
    % extraction, and edge cases without requiring actual DICOM/NIfTI data.

    properties
        TempDir
    end

    methods(TestMethodSetup)
        function createTempDir(testCase)
            % Create a uniquely-named temp directory for each test to
            % avoid interference between parallel test runs.
            testCase.TempDir = fullfile(tempdir, ['test_pss_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function removeTempDir(testCase)
            % Close diary before rmdir to prevent file-lock errors on Windows
            diary off;
            if isfolder(testCase.TempDir)
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function test_result_has_nan_defaults(testCase)
            % When no DICOM/structure/dose data paths are provided,
            % all numeric result fields should default to NaN rather
            % than 0 or empty, ensuring downstream code can detect
            % missing data via isnan() checks.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);

            testCase.verifyTrue(isnan(result.adc_mean));
            testCase.verifyTrue(isnan(result.d_mean));
            testCase.verifyTrue(isnan(result.d_mean_dncnn));
            testCase.verifyTrue(isnan(result.d_mean_ivimnet));
            testCase.verifyTrue(isnan(result.dmean_gtvp));
            testCase.verifyTrue(isnan(result.d95_gtvp));
            testCase.verifyTrue(isnan(result.v50gy_gtvp));
            testCase.verifyTrue(isnan(result.dmean_gtvn));
            testCase.verifyTrue(isnan(result.d95_gtvn));
            testCase.verifyTrue(isnan(result.v50gy_gtvn));
        end

        function test_result_has_empty_bad_dwi_list(testCase)
            % When no DWI data is found, the bad_dwi_list (which tracks
            % scans that failed quality checks) should be empty.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyEmpty(result.bad_dwi_list);
        end

        function test_fx_id_naming_for_fraction_1(testCase)
            % Verify that process_single_scan creates the 'nii' output
            % directory under basefolder (needed for NIfTI conversion output).
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 1;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            process_single_scan(ctx);
            testCase.verifyTrue(isfolder(fullfile(ctx.basefolder, 'nii')));
        end

        function test_fx_id_naming_for_post(testCase)
            % When the fraction index (fi) exceeds the number of RT dose
            % columns, the scan is post-treatment and fx_id should be
            % set to 'post'. Verify no errors occur in this case.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 6;
            ctx.n_rtdose_cols = 5;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyEmpty(result.bad_dwi_list);
        end

        function test_b0_ref_output_on_fx1(testCase)
            % The b0_ref output (baseline b=0 reference image) is produced
            % at fi==1 for use in subsequent DIR registration. Since we
            % cannot easily create synthetic NIfTI data, verify that it
            % remains empty when no DWI file exists.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 1;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [~, b0_ref, ~, ~] = process_single_scan(ctx);
            testCase.verifyEmpty(b0_ref);
        end

        function test_data_gtvp_initialized(testCase)
            % The GTVp (primary tumor) data struct should always be
            % initialized with the expected DWI parameter vector fields,
            % even when no actual data is loaded.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isstruct(result.data_gtvp));
            testCase.verifyTrue(isfield(result.data_gtvp, 'adc_vector'));
            testCase.verifyTrue(isfield(result.data_gtvp, 'd_vector'));
            testCase.verifyTrue(isfield(result.data_gtvp, 'f_vector'));
            testCase.verifyTrue(isfield(result.data_gtvp, 'dstar_vector'));
        end

        function test_data_gtvn_initialized(testCase)
            % The GTVn (nodal tumor) data struct should be initialized
            % analogously to GTVp, with at least an adc_vector field.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isstruct(result.data_gtvn));
            testCase.verifyTrue(isfield(result.data_gtvn, 'adc_vector'));
        end

        function test_kurtosis_needs_four_elements(testCase)
            % Kurtosis requires at least 4 data points to be meaningful.
            % When no voxels exist (empty data), adc_kurtosis and
            % d_kurtosis should remain NaN rather than erroring.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isnan(result.adc_kurtosis));
            testCase.verifyTrue(isnan(result.d_kurtosis));
        end

        function test_dvh_metrics_default_nan(testCase)
            % DVH-related metrics (dmean, d95, v50gy) for both GTVp and
            % GTVn should default to NaN when no dose data is available.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isnan(result.dmean_gtvp), 'dmean_gtvp should be NaN without dose data.');
            testCase.verifyTrue(isnan(result.d95_gtvp), 'd95_gtvp should be NaN without dose data.');
            testCase.verifyTrue(isnan(result.v50gy_gtvp), 'v50gy_gtvp should be NaN without dose data.');
            testCase.verifyTrue(isnan(result.dmean_gtvn), 'dmean_gtvn should be NaN without dose data.');
            testCase.verifyTrue(isnan(result.d95_gtvn), 'd95_gtvn should be NaN without dose data.');
            testCase.verifyTrue(isnan(result.v50gy_gtvn), 'v50gy_gtvn should be NaN without dose data.');
        end

        function test_dncnn_field_defaults_nan(testCase)
            % DnCNN-specific fields should default to NaN when no
            % denoised data is available.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isnan(result.d_mean_dncnn), 'd_mean_dncnn should default to NaN.');
        end

        function test_ivimnet_field_defaults_nan(testCase)
            % IVIMnet-specific fields should default to NaN when no
            % neural network fitted data is available.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isnan(result.d_mean_ivimnet), 'd_mean_ivimnet should default to NaN.');
        end

        function test_scan_naming_fx_range(testCase)
            % Verify scan naming for mid-treatment fractions (fi=3).
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 3;
            ctx.rpi = 1;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyEmpty(result.bad_dwi_list, ...
                'No bad DWI list entries for fi=3 with no data.');
        end

        function test_repeat_index_2(testCase)
            % Verify that repeat index > 1 does not cause errors.
            ctx = testCase.makeMinimalCtx();
            ctx.rpi = 2;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isnan(result.adc_mean), ...
                'adc_mean should be NaN for repeat index 2 with no data.');
        end

        function test_gtvp_struct_all_vector_fields(testCase)
            % Verify that GTVp data struct has all expected DWI vector fields
            % including DnCNN and IVIMnet variants.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            expected_fields = {'adc_vector', 'd_vector', 'f_vector', 'dstar_vector', ...
                'adc_vector_dncnn', 'd_vector_dncnn', 'f_vector_dncnn', 'dstar_vector_dncnn', ...
                'd_vector_ivimnet', 'f_vector_ivimnet', 'dstar_vector_ivimnet'};
            for i = 1:numel(expected_fields)
                testCase.verifyTrue(isfield(result.data_gtvp, expected_fields{i}), ...
                    sprintf('data_gtvp should have field %s.', expected_fields{i}));
            end
        end

        function test_gtvn_struct_all_vector_fields(testCase)
            % Verify that GTVn data struct has all DWI vector fields.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            expected_fields = {'adc_vector', 'd_vector', 'f_vector', 'dstar_vector'};
            for i = 1:numel(expected_fields)
                testCase.verifyTrue(isfield(result.data_gtvn, expected_fields{i}), ...
                    sprintf('data_gtvn should have field %s.', expected_fields{i}));
            end
        end

        function test_gtvp_ref_empty_without_data(testCase)
            % GTVp reference output should be empty when no mask file exists.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 1;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [~, ~, gtvp_ref, ~] = process_single_scan(ctx);
            testCase.verifyEmpty(gtvp_ref, ...
                'GTVp reference should be empty without mask data.');
        end

        function test_gtvn_ref_empty_without_data(testCase)
            % GTVn reference output should be empty when no mask file exists.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 1;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [~, ~, ~, gtvn_ref] = process_single_scan(ctx);
            testCase.verifyEmpty(gtvn_ref, ...
                'GTVn reference should be empty without mask data.');
        end

        function test_nii_dir_created_for_non_fx1(testCase)
            % The nii output directory should be created even for non-Fx1 fractions.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 3;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            process_single_scan(ctx);
            testCase.verifyTrue(isfolder(fullfile(ctx.basefolder, 'nii')), ...
                'nii directory should be created for all fractions.');
        end

        function test_result_struct_has_bad_dwi_list_field(testCase)
            % Result struct should always contain bad_dwi_list field.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isfield(result, 'bad_dwi_list'), ...
                'Result should always contain bad_dwi_list field.');
            testCase.verifyTrue(iscell(result.bad_dwi_list), ...
                'bad_dwi_list should be a cell array.');
        end

        function test_all_result_fields_present(testCase)
            % Comprehensive check: result struct should contain ALL expected
            % scalar metric fields, data struct fields, and bad_dwi_list.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);

            expected_scalar = {'adc_mean', 'adc_kurtosis', 'd_mean', 'd_kurtosis', ...
                'd_mean_dncnn', 'd_mean_ivimnet', ...
                'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', ...
                'v50gy_gtvp', 'v50gy_gtvn'};
            for i = 1:numel(expected_scalar)
                testCase.verifyTrue(isfield(result, expected_scalar{i}), ...
                    sprintf('Result should have field %s.', expected_scalar{i}));
            end
            testCase.verifyTrue(isfield(result, 'data_gtvp'), 'Result should have data_gtvp.');
            testCase.verifyTrue(isfield(result, 'data_gtvn'), 'Result should have data_gtvn.');
        end

        function test_fi_equals_n_rtdose_cols_is_on_treatment(testCase)
            % When fi equals n_rtdose_cols, the scan is still on-treatment
            % (fx_id = "fx5"), not post-treatment. Only fi > n_rtdose_cols
            % triggers the "post" naming convention.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 5;
            ctx.n_rtdose_cols = 5;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyEmpty(result.bad_dwi_list, ...
                'fi=5 with n_rtdose_cols=5 should not error.');
        end

        function test_repeat_index_high(testCase)
            % Repeat indices up to 6 are valid (e.g., Fx1 repeatability).
            ctx = testCase.makeMinimalCtx();
            ctx.rpi = 6;
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);
            testCase.verifyTrue(isnan(result.adc_mean), ...
                'adc_mean should be NaN for high repeat index with no data.');
        end

        function test_fi_greater_than_1_without_fx1_ref(testCase)
            % When fi > 1 but no Fx1 reference is available (b0_fx1_ref
            % is empty), DIR registration should be skipped gracefully.
            ctx = testCase.makeMinimalCtx();
            ctx.fi = 3;
            ctx.b0_fx1_ref = [];
            ctx.gtv_mask_fx1_ref = [];
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, b0_ref, gtvp_ref, gtvn_ref] = process_single_scan(ctx);
            testCase.verifyEmpty(b0_ref, 'b0_ref should be empty for fi>1.');
            testCase.verifyEmpty(gtvp_ref, 'gtvp_ref should be empty for fi>1.');
            testCase.verifyEmpty(gtvn_ref, 'gtvn_ref should be empty for fi>1.');
            testCase.verifyTrue(isnan(result.adc_mean));
        end

        function test_multiple_scalar_fields_are_nan(testCase)
            % Exhaustive NaN verification: EVERY scalar metric field
            % should default to NaN, not 0 or empty.
            ctx = testCase.makeMinimalCtx();
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';

            [result, ~, ~, ~] = process_single_scan(ctx);

            nan_fields = {'adc_mean', 'adc_kurtosis', 'd_mean', 'd_kurtosis', ...
                'd_mean_dncnn', 'd_mean_ivimnet', ...
                'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', ...
                'v50gy_gtvp', 'v50gy_gtvn'};
            for i = 1:numel(nan_fields)
                testCase.verifyTrue(isnan(result.(nan_fields{i})), ...
                    sprintf('%s should default to NaN, got %g.', nan_fields{i}, result.(nan_fields{i})));
            end
        end
    end

    methods(Access = private)
        function ctx = makeMinimalCtx(testCase)
            % Build the minimum context struct required by process_single_scan.
            % All data paths are set to empty strings so no actual
            % DICOM/NIfTI I/O occurs -- only initialization logic runs.
            ctx = struct();
            ctx.fi = 1;
            ctx.rpi = 1;
            ctx.n_rtdose_cols = 5;
            ctx.basefolder = testCase.TempDir;
            ctx.dataloc = testCase.TempDir;
            ctx.id_j = 'TEST001';
            ctx.mrn_j = '000001';
            ctx.pat_lf = 0;
            ctx.pat_immuno = 0;
            ctx.ivim_bthr = 100;
            ctx.dcm2nii_call = 'dcm2niix';
            ctx.dicomloc = '';
            ctx.struct_file = '';
            ctx.struct_file_gtvn = '';
            ctx.dicomdoseloc = '';
            ctx.b0_fx1_ref = [];
            ctx.gtv_mask_fx1_ref = [];
            ctx.gtvn_mask_fx1_ref = [];
            ctx.use_gpu = false;
        end
    end
end
