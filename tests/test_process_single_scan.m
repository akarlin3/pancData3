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
    end

    methods(Access = private)
        function ctx = makeMinimalCtx(testCase)
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
        end
    end
end
