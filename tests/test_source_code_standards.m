classdef test_source_code_standards < matlab.unittest.TestCase
    % TEST_SOURCE_CODE_STANDARDS Static code analysis tests extracted from
    % test_dwi_pipeline.m. These tests verify source-code patterns, naming
    % conventions, and architectural invariants without executing pipeline
    % logic.
    %
    % Run tests with:
    %   results = runtests('tests/test_source_code_standards.m');

    methods(TestMethodSetup)
        function addProjectPaths(testCase)
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods(Test)

        % ---- Elastic Net ------------------------------------------------

        function testElasticNet_AlphaIsHalf(testCase)
            % The refactored code must use Alpha = 0.5 for Elastic Net.
            code = metrics_code;
            matches = regexp(code, '''Alpha''\s*,\s*0\.5', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should contain Alpha = 0.5 for Elastic Net');
        end

        function testElasticNet_NotPureLasso(testCase)
            % Confirm no lassoglm call uses pure LASSO (Alpha = 1).
            code = metrics_code;
            matches = regexp(code, '''Alpha''\s*,\s*1(\s|,|\))', 'match');
            testCase.verifyTrue(isempty(matches), ...
            'No lassoglm call should use pure LASSO (Alpha=1)');
        end

        % ---- Loop bounds ------------------------------------------------

        function testLoop_NotHardcodedTo2And3(testCase)
            % The loop "for target_fx = [2, 3]" should have been replaced.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'target_fx = [2, 3]'), ...
            'Loop should NOT be hardcoded to [2, 3]');
        end

        function testLoop_UsesNTpUpperBound(testCase)
            % The loop should use nTp as the upper bound.
            code = metrics_code;
            matches = regexp(code, 'target_fx\s*=\s*2\s*:\s*nTp', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'Loop should iterate from 2 to nTp');
        end

        % ---- Removed sections -------------------------------------------

        function testRemoved_TargetedM8(testCase)
            % The "Targeted Statistical Analysis (m = 8 tests)" section
            % should no longer exist in metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'Targeted Statistical Analysis (m = 8'), ...
            'Targeted m=8 section should be removed');
        end

        function testRemoved_TargetedM4(testCase)
            % The "Targeted FDR (Relative Change Only, m = 4)" section
            % should no longer exist in metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'Targeted FDR (Relative Change Only, m = 4)'), ...
            'Targeted m=4 FDR section should be removed');
        end

        % ---- ADC patterns (metrics) -------------------------------------

        function testADC_NoInlineOLS(testCase)
            % The old inline OLS pattern "beta = X \ log(s_signal)" should
            % be gone from metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'X \ log(s_signal)'), ...
            'Inline OLS fitting should be removed from metrics.m');
        end

        function testADC_NoNestedPixelLoops(testCase)
            % The old nested "for y = 1:ny / for x = 1:nx" pixel loops
            % should be gone from metrics.m (replaced by fit_adc_mono call).
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'slice_adc(y,x) = beta(2)'), ...
            'Nested pixel-wise OLS loops should be removed');
        end

        % ---- FDR --------------------------------------------------------

        function testRetained_PerTimepointFDR(testCase)
            % The FDR section should now apply BH per-timepoint.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'Per-Timepoint') || contains(code, 'PER-TIMEPOINT'), ...
            'FDR section should indicate per-timepoint operation');
        end

        function testFDR_NoGlobalPooling(testCase)
            % The old global FDR sweep header should no longer exist.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'Full Metric Sweep'), ...
            'Global FDR sweep should be replaced with per-timepoint');
        end

        function testPerTimepointFDR_NoGlobalSweepInSection9(testCase)
            % Section 9 should no longer pool across timepoints.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'FDR Correction') || ...
            contains(code, 'Per-Timepoint'), ...
            'Section 9 FDR should be per-timepoint, not global');
        end

        % ---- B-value ----------------------------------------------------

        function testBval_NoTruncation(testCase)
            % The old truncation pattern should be gone from metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'min_dim = min(n_imgs, n_bvals)'), ...
            'Arbitrary b-value truncation (min_dim) should be removed');
        end

        % ---- Imputation (static) ----------------------------------------

        function testImpute_NoListwiseDeletion(testCase)
            % The old listwise deletion pattern should be gone from metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'valid_lasso_mask = all(~isnan(X_lasso), 2)'), ...
            'Listwise deletion pattern should be removed from metrics.m');
        end

        function testImpute_UsesKNN(testCase)
            % metrics.m should use knn_impute_train_test for imputation.
            code = metrics_code;
            matches = regexp(code, 'knn_impute_train_test\s*\(', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should use knn_impute_train_test for imputation');
        end

        function testImpute_ExcludesAllNaNRows(testCase)
            % metrics.m should exclude patients with ALL imaging data missing.
            code = metrics_code;
            matches = regexp(code, 'any\s*\(\s*~isnan\s*\(\s*X_lasso_all', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should check for patients missing all imaging data');
        end

        function testImpute_BeforeStandardization(testCase)
            % Verify imputation occurs before the Elastic Net standardization.
            % knn_impute_train_test must appear before lassoglm(..., 'Standardize', true)
            code = metrics_code;
            pos_fill = regexp(code, 'knn_impute_train_test\s*\(', 'start', 'once');
            pos_std  = regexp(code, '''Standardize''\s*,\s*true', 'start', 'once');
            testCase.verifyTrue(~isempty(pos_fill) && ~isempty(pos_std), ...
            'Both knn_impute_train_test and Standardize should exist in metrics.m');
            testCase.verifyTrue(pos_fill < pos_std, ...
            'Imputation (knn_impute_train_test) must occur before standardization');
        end

        % ---- ROC --------------------------------------------------------

        function testOOFROC_NoInSampleROC(testCase)
            % The old in-sample ROC pattern should be gone.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'EXPLORATORY ROC ANALYSIS'), ...
            'In-sample exploratory ROC should be replaced with OOF ROC');
        end

        function testOOFROC_LabeledAsPrimary(testCase)
            % ROC section should say PRIMARY ROC ANALYSIS in the fprintf.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'PRIMARY ROC ANALYSIS'), ...
            'ROC section should be labeled Primary');
        end

        function testOOFROC_NoStaleGLMRefit(testCase)
            % The stale fitglm call inside the ROC block should be gone.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'mdl_roc = fitglm'), ...
            'Stale fitglm call in ROC block should be removed');
        end

        function testOOFROC_YoudenCutoffFromLOOCV(testCase)
            % The Youden optimal cutoff (roc_opt_thresh) must come from perfcurve,
            % not from any logistic regression fitted probabilities.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'roc_opt_thresh'), ...
            'roc_opt_thresh must be defined (Youden cutoff from LOOCV ROC)');
            testCase.verifyTrue(~contains(code, 'mdl_roc.Fitted'), ...
            'Optimal cutoff must not use mdl_roc.Fitted (no in-sample refitting)');
        end

        % ---- Collinearity / feature -------------------------------------

        function testCollinearity_NoGlobalLeakage(testCase)
            % filter_collinear_features should NOT be called on global matrices like X_clean_all or X_impute.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'filter_collinear_features(X_clean_all'), ...
            'filter_collinear_features should not be called on the full cohort X_clean_all');
            testCase.verifyTrue(~contains(code, 'filter_collinear_features(X_impute'), ...
            'filter_collinear_features should not be called on the full imputed matrix X_impute');
            % It SHOULD be called on training folds
            testCase.verifyTrue(contains(code, 'filter_collinear_features(X_tr_imp'), ...
            'filter_collinear_features must be applied to training folds (X_tr_imp)');
        end

        function testSubVol_DIRReenabled(testCase)
            % The old keep_policy = 1:10 exclusion should be gone.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'keep_policy = 1:10'), ...
            'Sub-volume features should no longer be excluded now that DIR is implemented');
        end

        function testSubVol_PostRTStillExcluded(testCase)
            % Post-RT dose exclusion (target_fx == 6) should still be present.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'target_fx == 6'), ...
            'Post-RT dose exclusion should be retained');
        end

        % ---- Sanity -----------------------------------------------------

        function testSanityCheck_NoKurtosisPlot(testCase)
            % The sanity check subplot must NOT use adc_kurt as a metric.
            % It should use adc_sd instead.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'kurt_fx1 = adc_kurt'), ...
            'Sanity check should use ADC SD, not trace-average ADC Kurtosis');
            testCase.verifyTrue(contains(code, 'sd_fx1  = adc_sd'), ...
            'Sanity check heterogeneity subplot should use adc_sd');
        end

        % ---- Visualization ----------------------------------------------

        function testVis_NoTruncation(testCase)
            % The old truncation pattern should be gone from visualize_results.m.
            code = vis_code;
            testCase.verifyTrue(~contains(code, 'min_dim = min('), ...
            'Arbitrary b-value truncation (min_dim) should be removed from visualize_results.m');
        end

        function testVis_ExpectedProtocolDefined(testCase)
            % plot_parameter_maps.m should define expected_bvals for protocol validation.
            code = param_maps_code;
            matches = regexp(code, 'expected_bvals\s*=\s*\[0;\s*30;\s*150;\s*550\]', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'plot_parameter_maps.m should define expected_bvals = [0; 30; 150; 550]');
        end

        function testVis_ProtocolDeviationFlagged(testCase)
            % plot_parameter_maps.m should flag protocol deviations with a warning message.
            code = param_maps_code;
            testCase.verifyTrue(contains(lower(code), 'protocol deviation'), ...
            'Protocol deviation flagging should be present in plot_parameter_maps.m');
        end

        function testVis_DeviationExcludesPatient(testCase)
            % The validation block should use 'continue' to exclude
            % deviating patients from the comparative mapping.
            code = param_maps_code;
            idx_dev = regexpi(code, 'protocol deviation', 'start');
            idx_continue = strfind(code, 'continue');
            testCase.verifyTrue(~isempty(idx_dev) && ~isempty(idx_continue), ...
            'Both protocol deviation flag and continue must exist in plot_parameter_maps.m');
            if isempty(idx_dev) || isempty(idx_continue)
                return;
            end
            % At least one continue must follow the protocol deviation flag
            % (within 200 chars, i.e. within the same if-block)
            found = false;
            for ci = 1:length(idx_continue)
            if idx_continue(ci) > idx_dev(1) && (idx_continue(ci) - idx_dev(1)) < 200
            found = true;
            break;
            end
            end
            testCase.verifyTrue(found, ...
            'A continue statement should follow the protocol deviation flag to exclude the patient');
        end

        % ---- DIR source patterns ----------------------------------------

        function testDIR_FunctionExistsInCodebase(testCase)
            % load_dwi_data_forAvery.m should call apply_dir_mask_propagation.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'apply_dir_mask_propagation'), ...
            'load_dwi_data_forAvery.m should call apply_dir_mask_propagation');
        end

        function testDIR_ReturnsDfieldAndRef(testCase)
            % apply_dir_mask_propagation should return 3 outputs (mask, D_forward, ref3d).
            code = loaddwi_code();
            % Check signature via file reading
            apPath = fullfile(fileparts(mfilename('fullpath')), '..', 'utils', 'apply_dir_mask_propagation.m');
            fid = fopen(apPath, 'r');
            if fid == -1, return; end
            ap_code = fread(fid,'*char')'; fclose(fid);
            testCase.verifyTrue(contains(ap_code, 'D_forward, ref3d] = apply_dir_mask_propagation'), ...
            'apply_dir_mask_propagation should return D_forward and ref3d');
        end

        function testDIR_DoseMapRemainsRigid(testCase)
            % load_dwi_data_forAvery.m should NOT apply imwarp to the dose map.
            code = loaddwi_code;
            testCase.verifyTrue(~contains(code, 'imwarp(dose_map, D_forward_cur'), ...
            'Dose map must remain rigidly aligned and should NOT be warped via imwarp');
        end

        function testDIR_DoseWarpFallback(testCase)
            % At baseline (fi=1) or when DIR fails, code must fall back to rigid dose_map.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'dose_map_dvh     = dose_map'), ...
            'Rigid dose fallback ''dose_map_dvh = dose_map'' should be initialised before DIR block');
        end

        function testDIR_DfieldCached(testCase)
            % D_forward and ref3d should be written to the .mat cache file.
            code = loaddwi_code;
            testCase.verifyTrue( ...
                contains(code, 'parsave_dir_cache(dir_cache_file, gtv_mask_warped, D_forward, ref3d)') || ...
                contains(code, 'parsave_dir_cache(dir_cache_file, gtv_mask_warped, D_forward_cur, ref3d_cur)'), ...
            'D_forward and ref3d must be included in the cache save call');
        end

        function testDIR_GTVnDoseRemainsRigid(testCase)
            % GTVn dose block must also use the rigid dose (dose_map_dvh_n).
            code = loaddwi_code;
            testCase.verifyTrue(~contains(code, 'dose_map_dvh_n = imwarp(dose_map'), ...
            'GTVn dose block must not warp dose_map into dose_map_dvh_n');
        end

        function testDIR_LinearMaskWarping(testCase)
            % apply_dir_mask_propagation.m should use 'linear' interpolation for mask warping.
            apPath = fullfile(fileparts(mfilename('fullpath')), '..', 'utils', 'apply_dir_mask_propagation.m');
            fid = fopen(apPath, 'r');
            if fid == -1, return; end
            ap_code = fread(fid,'*char')'; fclose(fid);
            testCase.verifyTrue(contains(ap_code, "'Interp', 'linear'"), ...
            "apply_dir_mask_propagation should use 'linear' interpolation for mask warping");
            testCase.verifyTrue(~contains(ap_code, "'Interp', 'nearest'"), ...
            "apply_dir_mask_propagation should NOT use 'nearest' interpolation for mask warping");
        end

        function testDIR_SymmetricWorkflow(testCase)
            % apply_dir_mask_propagation.m should implement the symmetric halfway-space logic.
            apPath = fullfile(fileparts(mfilename('fullpath')), '..', 'utils', 'apply_dir_mask_propagation.m');
            fid = fopen(apPath, 'r');
            if fid == -1, return; end
            ap_code = fread(fid,'*char')'; fclose(fid);
            testCase.verifyTrue(contains(ap_code, 'mid_img_refined'), ...
            "apply_dir_mask_propagation should use a refined midpoint image");
            testCase.verifyTrue(contains(ap_code, 'D_forward = D_forward_mid - D_backward_mid'), ...
            "apply_dir_mask_propagation should compute the symmetric proxy field");
        end

        % ---- ADC source patterns ----------------------------------------

        function testADC_NegativeClampedToNaN(testCase)
            % load_dwi_data_forAvery.m should assign NaN to negative ADC.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'adc_vals(adc_vals < 0) = nan') || ...
            contains(code, 'adc_vec_d(adc_vec_d < 0) = nan') || ...
            contains(code, 'adc_vals(adc_vals<0)=nan'), ...
            'Negative ADC should be set to NaN, not zero');
        end

        function testADC_NotClampedToZero(testCase)
            % The old clamping pattern should not exist.
            code = loaddwi_code;
            testCase.verifyTrue(~contains(code, 'adc_vec(adc_vec < 0) = 0'), ...
            'Negative ADC should NOT be clamped to zero');
        end

        % ---- IVIM source patterns ---------------------------------------

        function testIVIM_Bthr100(testCase)
            % Pipeline uses config struct for bthr.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'ivim_bthr = config_struct.ivim_bthr') || contains(code, 'ivim_bthr = 100'), ...
            'ivim_bthr should be read from config_struct in load_dwi_data.m');
        end

        function testIVIM_fUpperBound(testCase)
            % IVIMmodelfit.m should cap perfusion fraction f at 0.4.
            code = ivim_code;
            testCase.verifyTrue(contains(code, '0.4'), ...
            'IVIMmodelfit.m should set f upper bound to 0.4');
        end

        function testIVIM_DstarUpperBound(testCase)
            % IVIMmodelfit.m should cap D* at 0.1 mm²/s.
            code = ivim_code;
            testCase.verifyTrue(contains(code, '0.1'), ...
            'IVIMmodelfit.m should set D* upper bound to 0.1 mm²/s');
        end

        function testIVIM_NoLooseBoundsF(testCase)
            % Old loose f upper bound of 1 should be gone.
            code = ivim_code;
            % The old line was: lim = [0 0 0 0;3e-3 2*max(Y(:)) 1 1]
            testCase.verifyTrue(~contains(code, '3e-3 2*max(Y(:)) 1 1'), ...
            'Old loose f and D* limits should be removed from IVIMmodelfit.m');
        end

        function testIVIM_OptsLimSupported(testCase)
            % IVIMmodelfit.m should support opts.lim for caller overrides.
            code = ivim_code;
            testCase.verifyTrue(contains(code, '''lim''') && contains(code, 'opts.lim'), ...
            'IVIMmodelfit.m should support opts.lim override');
        end

        function testIVIM_DefaultBlim100(testCase)
            % IVIMmodelfit.m should default blim to 100.
            code = ivim_code;
            testCase.verifyTrue(contains(code, 'blim = 100'), ...
            'IVIMmodelfit.m default blim should be 100');
        end

    end
end

% Local helper functions to read source code for static analysis tests

function code = metrics_code()
    files = {'metrics_baseline.m', 'metrics_longitudinal.m', 'metrics_dosimetry.m', ...
             'metrics_stats_comparisons.m', 'metrics_stats_predictive.m', 'metrics_survival.m'};

    code = '';
    for i = 1:numel(files)
        filepath = fullfile(fileparts(mfilename('fullpath')), '..', 'core', files{i});
        fid = fopen(filepath, 'r');
        if fid ~= -1
            file_code = fread(fid, '*char')';
            fclose(fid);
            code = [code newline file_code]; %#ok<AGROW>
        end
    end
end

function code = loaddwi_code()
    files = {'load_dwi_data.m', 'discover_patient_files.m', 'compute_summary_metrics.m', 'process_single_scan.m', 'fit_models.m'};

    code = '';
    for i = 1:numel(files)
        filepath = fullfile(fileparts(mfilename('fullpath')), '..', 'core', files{i});
        fid = fopen(filepath, 'r');
        if fid ~= -1
            file_code = fread(fid, '*char')';
            fclose(fid);
            code = [code newline file_code]; %#ok<AGROW>
        end
    end
end

function code = ivim_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'dependencies', 'IVIMmodelfit.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = vis_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'core', 'visualize_results.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = param_maps_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'core', 'plot_parameter_maps.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end
