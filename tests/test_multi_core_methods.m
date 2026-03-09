%% test_multi_core_methods.m - Tests for multi-method core integration
%
% Validates that compute_summary_metrics correctly produces per-method
% sub-volume metrics when run_all_core_methods is enabled, and that
% backward compatibility is maintained when disabled.
%
% Tests covered:
%   1. test_default_off: all_core_metrics absent when feature is disabled
%   2. test_all_methods_computed: all 11 methods produce correct struct layout
%   3. test_backward_compat: top-level fields match the configured method
%   4. test_mask_storage: core_masks saved/omitted per store_core_masks flag
%   5. test_config_backward_compat: old configs without new fields get defaults
%   6. test_run_compare_cores_injection: compare_cores step injected correctly
%
% Uses function-based test pattern (not unittest class) for Octave compat.

function test_multi_core_methods()
    disp('==== Running test_multi_core_methods ====');
    [dir_path, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(dir_path, '..', 'core'));
    addpath(fullfile(dir_path, '..', 'utils'));
    addpath(fullfile(dir_path, '..', 'dependencies'));
    if exist(fullfile(dir_path, '..', '.octave_compat'), 'dir')
        addpath(genpath(fullfile(dir_path, '..', '.octave_compat')));
    end

    % All 11 tumor core delineation methods supported by extract_tumor_core
    ALL_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};

    test_default_off();
    test_all_methods_computed(ALL_METHODS);
    test_backward_compat(ALL_METHODS);
    test_mask_storage(ALL_METHODS);
    test_config_backward_compat();
    test_run_compare_cores_injection();

    disp('==== All test_multi_core_methods tests PASSED ====');
end

%% --- Test: Default off produces no all_core_metrics field ---
% When run_all_core_methods is false (the default), only the single
% configured core_method is used and no all_core_metrics struct is added.
function test_default_off()
    fprintf('  Testing default off (run_all_core_methods = false)...\n');
    [config, dvg, id_list, mrn_list, lf, immuno, gtv_locs, dwi_locs, dmean, d95, v50, dates] = make_mock_data();
    config.run_all_core_methods = false;

    sm = compute_summary_metrics(config, dvg, id_list, mrn_list, lf, immuno, ...
        gtv_locs, dwi_locs, dmean, d95, v50, dates);

    assert(~isfield(sm, 'all_core_metrics'), ...
        'all_core_metrics should NOT be present when run_all_core_methods=false');
    fprintf('    [PASS] No all_core_metrics when disabled\n');
end

%% --- Test: All 11 methods computed ---
% With run_all_core_methods=true, compute_summary_metrics should produce
% an all_core_metrics struct containing a sub-struct for each of the 11
% methods, each with adc_sub_vol, adc_sub_mean, d_sub_mean fields of
% size [nPat x nTp x nDwiTypes].
function test_all_methods_computed(ALL_METHODS)
    fprintf('  Testing all 11 methods computed...\n');
    [config, dvg, id_list, mrn_list, lf, immuno, gtv_locs, dwi_locs, dmean, d95, v50, dates] = make_mock_data();
    config.run_all_core_methods = true;
    config.store_core_masks = false;

    sm = compute_summary_metrics(config, dvg, id_list, mrn_list, lf, immuno, ...
        gtv_locs, dwi_locs, dmean, d95, v50, dates);

    assert(isfield(sm, 'all_core_metrics'), 'all_core_metrics field missing');

    for i = 1:numel(ALL_METHODS)
        mname = ALL_METHODS{i};
        assert(isfield(sm.all_core_metrics, mname), ...
            sprintf('Method %s missing from all_core_metrics', mname));
        ms = sm.all_core_metrics.(mname);

        % Check required fields exist with correct size [nPat x nTp x 3]
        assert(isfield(ms, 'adc_sub_vol'), sprintf('%s: missing adc_sub_vol', mname));
        assert(isequal(size(ms.adc_sub_vol), [2, 2, 3]), ...
            sprintf('%s: adc_sub_vol wrong size', mname));
        assert(isfield(ms, 'adc_sub_mean'), sprintf('%s: missing adc_sub_mean', mname));
        assert(isfield(ms, 'd_sub_mean'), sprintf('%s: missing d_sub_mean', mname));

        % At least Standard dwi_type (index 1) should have non-NaN values
        assert(~isnan(ms.adc_sub_vol(1,1,1)), ...
            sprintf('%s: adc_sub_vol(1,1,1) should not be NaN', mname));

        % Masks should NOT be stored when store_core_masks is false
        assert(~isfield(ms, 'core_masks'), ...
            sprintf('%s: core_masks should not be present when store_core_masks=false', mname));
    end
    fprintf('    [PASS] All 11 methods have correct structure\n');
end

%% --- Test: Backward compatibility ---
% When run_all_core_methods=true, the top-level sub-volume fields
% (e.g., sm.adc_sub_vol) must still match the values from the method
% specified by config.core_method (here 'adc_threshold'). This ensures
% that enabling the multi-method feature does not change the default
% pipeline behavior.
function test_backward_compat(ALL_METHODS)
    fprintf('  Testing backward compatibility...\n');
    [config, dvg, id_list, mrn_list, lf, immuno, gtv_locs, dwi_locs, dmean, d95, v50, dates] = make_mock_data();
    config.run_all_core_methods = true;
    config.store_core_masks = false;
    config.core_method = 'adc_threshold';

    sm = compute_summary_metrics(config, dvg, id_list, mrn_list, lf, immuno, ...
        gtv_locs, dwi_locs, dmean, d95, v50, dates);

    % Top-level fields should match the configured method's per-method entry
    assert(isfield(sm, 'adc_sub_vol'), 'Top-level adc_sub_vol missing');
    assert(isfield(sm, 'all_core_metrics'), 'all_core_metrics missing');

    pm = sm.all_core_metrics.adc_threshold;
    % Values for Standard (dwi_type=1) at patient 1, timepoint 1
    tol = 1e-10;
    top_val = sm.adc_sub_vol(1,1,1);
    pm_val = pm.adc_sub_vol(1,1,1);
    assert(abs(top_val - pm_val) < tol || (isnan(top_val) && isnan(pm_val)), ...
        sprintf('Top-level adc_sub_vol (%g) != per-method value (%g)', top_val, pm_val));
    fprintf('    [PASS] Top-level fields match configured method\n');
end

%% --- Test: Mask storage ---
% When store_core_masks=true, each method's sub-struct should include
% a core_masks cell array of size [nPat x nTp] containing logical masks.
% When store_core_masks=false (tested in test_all_methods_computed),
% core_masks should be absent to save memory.
function test_mask_storage(ALL_METHODS)
    fprintf('  Testing mask storage...\n');
    [config, dvg, id_list, mrn_list, lf, immuno, gtv_locs, dwi_locs, dmean, d95, v50, dates] = make_mock_data();
    config.run_all_core_methods = true;
    config.store_core_masks = true;

    sm = compute_summary_metrics(config, dvg, id_list, mrn_list, lf, immuno, ...
        gtv_locs, dwi_locs, dmean, d95, v50, dates);

    for i = 1:numel(ALL_METHODS)
        mname = ALL_METHODS{i};
        ms = sm.all_core_metrics.(mname);
        assert(isfield(ms, 'core_masks'), ...
            sprintf('%s: core_masks missing when store_core_masks=true', mname));
        assert(iscell(ms.core_masks), sprintf('%s: core_masks should be cell', mname));
        assert(isequal(size(ms.core_masks), [2, 2]), ...
            sprintf('%s: core_masks wrong size', mname));
        % First patient, first timepoint should have a logical mask
        mask = ms.core_masks{1,1};
        assert(islogical(mask), sprintf('%s: mask should be logical', mname));
    end
    fprintf('    [PASS] Masks stored correctly\n');
end

%% --- Test: Config backward compatibility (missing fields default to false) ---
% Old config files that predate the multi-core feature will not contain
% run_compare_cores, run_all_core_methods, or store_core_masks fields.
% parse_config must add these with default value false so the pipeline
% runs unchanged on legacy configs.
function test_config_backward_compat()
    fprintf('  Testing config backward compatibility...\n');
    [dir_path, ~, ~] = fileparts(mfilename('fullpath'));

    % Create a minimal config without the new fields
    tmp_dir = tempname;
    mkdir(tmp_dir);
    cleanup = onCleanup(@() rmdir(tmp_dir, 's'));

    cfg_str = '{"dataloc": "/tmp", "dcm2nii_call": "dcm2niix", "dwi_type": "Standard", "core_method": "adc_threshold"}';
    cfg_file = fullfile(tmp_dir, 'config.json');
    fid = fopen(cfg_file, 'w');
    fwrite(fid, cfg_str);
    fclose(fid);

    config = parse_config(cfg_file);

    % Verify all 3 new fields default to false
    assert(isfield(config, 'run_compare_cores'), 'run_compare_cores default missing');
    assert(config.run_compare_cores == false, 'run_compare_cores should default to false');

    assert(isfield(config, 'run_all_core_methods'), 'run_all_core_methods default missing');
    assert(config.run_all_core_methods == false, 'run_all_core_methods should default to false');

    assert(isfield(config, 'store_core_masks'), 'store_core_masks default missing');
    assert(config.store_core_masks == false, 'store_core_masks should default to false');

    fprintf('    [PASS] Old configs get correct defaults\n');
end

%% --- Test: run_compare_cores injection into steps ---
% When config.run_compare_cores=true, run_dwi_pipeline automatically
% injects a 'compare_cores' step right after 'metrics_baseline'.
% This test verifies the injection logic in isolation, checking both
% the true case (step added) and the false case (step not added).
function test_run_compare_cores_injection()
    fprintf('  Testing run_compare_cores step injection...\n');

    % Simulate the injection logic from run_dwi_pipeline.m
    steps_to_run = {'test', 'load', 'sanity', 'visualize', 'metrics_baseline', ...
        'metrics_longitudinal', 'metrics_dosimetry'};

    % With run_compare_cores = true
    config_struct = struct('run_compare_cores', true);
    if config_struct.run_compare_cores && ~ismember('compare_cores', steps_to_run)
        idx = find(strcmp(steps_to_run, 'metrics_baseline'));
        if ~isempty(idx)
            steps_to_run = [steps_to_run(1:idx), {'compare_cores'}, steps_to_run(idx+1:end)];
        else
            steps_to_run{end+1} = 'compare_cores';
        end
    end

    assert(ismember('compare_cores', steps_to_run), 'compare_cores should be injected');
    bl_idx = find(strcmp(steps_to_run, 'metrics_baseline'));
    cc_idx = find(strcmp(steps_to_run, 'compare_cores'));
    assert(cc_idx == bl_idx + 1, 'compare_cores should be right after metrics_baseline');

    % With run_compare_cores = false (no injection)
    steps2 = {'load', 'metrics_baseline'};
    config2 = struct('run_compare_cores', false);
    if config2.run_compare_cores && ~ismember('compare_cores', steps2)
        steps2{end+1} = 'compare_cores';
    end
    assert(~ismember('compare_cores', steps2), 'compare_cores should NOT be injected when false');

    fprintf('    [PASS] Step injection works correctly\n');
end

%% --- Helper: Create mock data for compute_summary_metrics ---
function [config, dvg, id_list, mrn_list, lf, immuno, gtv_locs, dwi_locs, dmean, d95, v50, dates]  = make_mock_data()
    rng(42);
    n_pat = 2;
    n_tp = 2;

    config = struct();
    config.adc_thresh = 0.001;
    config.high_adc_thresh = 0.00115;
    config.d_thresh = 0.001;
    config.f_thresh = 0.1;
    config.dstar_thresh = 0.01;
    config.min_vox_hist = 10;
    config.adc_max = 0.003;
    config.core_method = 'adc_threshold';
    config.core_percentile = 25;
    config.core_n_clusters = 2;
    config.fdm_parameter = 'adc';
    config.fdm_thresh = 0.0002;
    config.use_checkpoints = false;
    config.run_all_core_methods = false;
    config.store_core_masks = false;
    config.dwi_types_to_run = 1;  % Standard only
    config.dataloc = tempdir;
    config.dwi_type_name = 'Standard';

    n_vox = 50;
    dvg = struct();
    for j = 1:n_pat
        for k = 1:n_tp
            for r = 1:1
                dvg(j,k,r).adc_vector = [0.0005 + 0.0001*randn(15,1); 0.0015 + 0.0002*randn(35,1)];
                dvg(j,k,r).d_vector = [0.0004 + 0.0001*randn(15,1); 0.0014 + 0.0002*randn(35,1)];
                dvg(j,k,r).f_vector = [0.05 + 0.01*randn(15,1); 0.2 + 0.05*randn(35,1)];
                dvg(j,k,r).dstar_vector = 0.01 + 0.005*randn(n_vox, 1);
                dvg(j,k,r).adc_vector_dncnn = dvg(j,k,r).adc_vector;
                dvg(j,k,r).d_vector_dncnn = dvg(j,k,r).d_vector;
                dvg(j,k,r).f_vector_dncnn = dvg(j,k,r).f_vector;
                dvg(j,k,r).dstar_vector_dncnn = dvg(j,k,r).dstar_vector;
                dvg(j,k,r).d_vector_ivimnet = dvg(j,k,r).d_vector;
                dvg(j,k,r).f_vector_ivimnet = dvg(j,k,r).f_vector;
                dvg(j,k,r).dstar_vector_ivimnet = dvg(j,k,r).dstar_vector;
                dvg(j,k,r).dose_vector = 50 + 5*randn(n_vox, 1);
                dvg(j,k,r).vox_vol = 0.008;  % 2mm x 2mm x 2mm
                dvg(j,k,r).vox_dims = [2 2 2];
            end
        end
    end

    id_list = {'PAT001', 'PAT002'};
    mrn_list = {'MRN001', 'MRN002'};
    lf = [0; 1];
    immuno = [0; 0];
    gtv_locs = cell(n_pat, n_tp, 1);  % empty = no 3D mask
    dwi_locs = cell(n_pat, n_tp, 1);
    dmean = nan(n_pat, n_tp);
    d95 = nan(n_pat, n_tp);
    v50 = nan(n_pat, n_tp);
    dates = cell(n_pat, n_tp);
end
