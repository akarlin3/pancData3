classdef test_metrics_dosimetry < matlab.unittest.TestCase
    % TEST_METRICS_DOSIMETRY Unit tests for metrics_dosimetry.
    %
    % Verifies output dimensions, NaN-passthrough when dose/ADC vectors are
    % empty, correct dispatch of dtype_idx to the right vector field, and
    % that the GTV mask cache path does not crash on absent files.

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    % ------------------------------------------------------------------ %
    %  Helpers                                                            %
    % ------------------------------------------------------------------ %
    methods(Access = private)
        function [m_id_list, id_list, nTp, config_struct, m_data_vectors_gtvp, gtv_locations] = ...
                buildEmptyVectorInputs(~, n_valid, n_total, nTp_val)
            % Builds inputs where every dose_vector / adc_vector is empty,
            % so the inner `if ~isempty(dose_vec) && ~isempty(adc_vec)` guard
            % is never entered and all outputs stay NaN.
            m_id_list = arrayfun(@(x) sprintf('Pt%02d', x), 1:n_valid, 'UniformOutput', false)';
            id_list   = arrayfun(@(x) sprintf('Pt%02d', x), 1:n_total, 'UniformOutput', false)';
            nTp = nTp_val;

            config_struct = struct( ...
                'adc_thresh',          1.15e-3, ...
                'd_thresh',            0.9e-3, ...
                'f_thresh',            0.15, ...
                'dstar_thresh',        0.05, ...
                'dwi_types_to_run',    1);

            % Build an n_valid × nTp × 1 struct array with empty vectors.
            empty_entry = struct( ...
                'adc_vector',       [], ...
                'd_vector',         [], ...
                'f_vector',         [], ...
                'dstar_vector',     [], ...
                'adc_vector_dncnn', [], ...
                'd_vector_dncnn',   [], ...
                'f_vector_dncnn',   [], ...
                'dstar_vector_dncnn',[], ...
                'd_vector_ivimnet', [], ...
                'f_vector_ivimnet', [], ...
                'dstar_vector_ivimnet', [], ...
                'dose_vector',      []);
            m_data_vectors_gtvp = repmat(empty_entry, n_valid, nTp, 1);

            % No GTV mask files on disk.
            gtv_locations = repmat({''}, n_total, nTp, 1);
        end

        function entry = makeFilledEntry(~, n_vox)
            % Creates a struct entry with synthetic non-empty vectors.
            adc  = rand(n_vox, 1) * 2e-3;
            d    = rand(n_vox, 1) * 1e-3;
            f    = rand(n_vox, 1) * 0.3;
            ds   = rand(n_vox, 1) * 0.05;
            dose = rand(n_vox, 1) * 60 + 10;

            entry = struct( ...
                'adc_vector',           adc,  ...
                'd_vector',             d,    ...
                'f_vector',             f,    ...
                'dstar_vector',         ds,   ...
                'adc_vector_dncnn',     adc,  ...
                'd_vector_dncnn',       d,    ...
                'f_vector_dncnn',       f,    ...
                'dstar_vector_dncnn',   ds,   ...
                'd_vector_ivimnet',     d,    ...
                'f_vector_ivimnet',     f,    ...
                'dstar_vector_ivimnet', ds,   ...
                'dose_vector',          dose);
        end
    end

    methods(Test)

        function testOutputDimensionsMatchInputs(testCase)
            % With 3 valid patients and 4 timepoints the 8 output arrays must
            % each be 3 × 4.
            n_valid  = 3;
            n_total  = 3;
            nTp_val  = 4;
            [m_id, id, nTp, cfg, dv, gtv] = ...
                testCase.buildEmptyVectorInputs(n_valid, n_total, nTp_val);

            [d95a, v50a, d95d, v50d, d95f, v50f, d95ds, v50ds] = ...
                metrics_dosimetry(m_id, id, nTp, cfg, dv, gtv);

            expected_size = [n_valid, nTp_val];
            testCase.verifyEqual(size(d95a),  expected_size, 'd95_adc_sub size');
            testCase.verifyEqual(size(v50a),  expected_size, 'v50_adc_sub size');
            testCase.verifyEqual(size(d95d),  expected_size, 'd95_d_sub size');
            testCase.verifyEqual(size(v50d),  expected_size, 'v50_d_sub size');
            testCase.verifyEqual(size(d95f),  expected_size, 'd95_f_sub size');
            testCase.verifyEqual(size(v50f),  expected_size, 'v50_f_sub size');
            testCase.verifyEqual(size(d95ds), expected_size, 'd95_dstar_sub size');
            testCase.verifyEqual(size(v50ds), expected_size, 'v50_dstar_sub size');
        end

        function testEmptyVectorsYieldAllNaN(testCase)
            % When dose_vector / adc_vector are empty the guard is not entered
            % and all outputs stay at their NaN defaults.
            [m_id, id, nTp, cfg, dv, gtv] = ...
                testCase.buildEmptyVectorInputs(2, 2, 3);

            [d95a, v50a, ~, ~, ~, ~, ~, ~] = ...
                metrics_dosimetry(m_id, id, nTp, cfg, dv, gtv);

            testCase.verifyTrue(all(isnan(d95a(:))), ...
                'd95_adc_sub should be all-NaN when vectors are empty.');
            testCase.verifyTrue(all(isnan(v50a(:))), ...
                'v50_adc_sub should be all-NaN when vectors are empty.');
        end

        function testDtype1UsesStandardVectors(testCase)
            % dtype_idx = 1 → adc_vector and d_vector used.
            % With 200 voxels all below threshold and dose all above 50 Gy,
            % the output values should be non-NaN.
            n_vox = 200;
            entry = testCase.makeFilledEntry(n_vox);
            % Ensure all values are below threshold so the sub-volume fills.
            entry.adc_vector   = zeros(n_vox, 1);
            entry.d_vector     = zeros(n_vox, 1);
            entry.f_vector     = zeros(n_vox, 1);
            entry.dstar_vector = zeros(n_vox, 1);
            entry.dose_vector  = 60 * ones(n_vox, 1);

            m_id  = {'Pt01'};
            id    = {'Pt01'};
            nTp   = 1;
            cfg   = struct( ...
                'adc_thresh', 1.0, 'd_thresh', 1.0, ...
                'f_thresh', 1.0, 'dstar_thresh', 1.0, ...
                'dwi_types_to_run', 1);
            dv(1, 1, 1) = entry;
            gtv = {''};

            [d95a, v50a, ~, ~, ~, ~, ~, ~] = metrics_dosimetry(m_id, id, nTp, cfg, dv, gtv);

            testCase.verifyFalse(isnan(d95a(1, 1)), ...
                'd95_adc_sub should be computed when vectors have enough voxels.');
            testCase.verifyFalse(isnan(v50a(1, 1)), ...
                'v50_adc_sub should be computed when vectors have enough voxels.');
        end

        function testDtype2UsesDnCNNVectors(testCase)
            % dtype_idx = 2 → adc_vector_dncnn, d_vector_dncnn, etc.
            n_vox = 200;
            entry = testCase.makeFilledEntry(n_vox);
            entry.adc_vector_dncnn   = zeros(n_vox, 1);
            entry.d_vector_dncnn     = zeros(n_vox, 1);
            entry.f_vector_dncnn     = zeros(n_vox, 1);
            entry.dstar_vector_dncnn = zeros(n_vox, 1);
            entry.dose_vector        = 60 * ones(n_vox, 1);

            cfg = struct( ...
                'adc_thresh', 1.0, 'd_thresh', 1.0, ...
                'f_thresh', 1.0, 'dstar_thresh', 1.0, ...
                'dwi_types_to_run', 2);
            dv(1, 1, 1) = entry;

            [d95a, v50a, ~, ~, ~, ~, ~, ~] = ...
                metrics_dosimetry({'Pt01'}, {'Pt01'}, 1, cfg, dv, {''});

            testCase.verifyFalse(isnan(d95a(1, 1)), ...
                'dtype=2: d95_adc_sub should use dncnn vectors.');
            testCase.verifyFalse(isnan(v50a(1, 1)), ...
                'dtype=2: v50_adc_sub should use dncnn vectors.');
        end

        function testDtype3UsesIVIMnetVectors(testCase)
            % dtype_idx = 3 → adc_vector (same as Standard) but d/f/dstar use ivimnet.
            n_vox = 200;
            entry = testCase.makeFilledEntry(n_vox);
            entry.adc_vector         = zeros(n_vox, 1);
            entry.d_vector_ivimnet   = zeros(n_vox, 1);
            entry.f_vector_ivimnet   = zeros(n_vox, 1);
            entry.dstar_vector_ivimnet = zeros(n_vox, 1);
            entry.dose_vector        = 60 * ones(n_vox, 1);

            cfg = struct( ...
                'adc_thresh', 1.0, 'd_thresh', 1.0, ...
                'f_thresh', 1.0, 'dstar_thresh', 1.0, ...
                'dwi_types_to_run', 3);
            dv(1, 1, 1) = entry;

            [d95a, ~, d95d, ~, ~, ~, ~, ~] = ...
                metrics_dosimetry({'Pt01'}, {'Pt01'}, 1, cfg, dv, {''});

            testCase.verifyFalse(isnan(d95a(1, 1)), ...
                'dtype=3: d95_adc_sub should use adc_vector.');
            testCase.verifyFalse(isnan(d95d(1, 1)), ...
                'dtype=3: d95_d_sub should use d_vector_ivimnet.');
        end

        function testAbsentGtvMatFileDoesNotCrash(testCase)
            % gtv_locations contains a non-existent path.  The function should
            % fall back gracefully (has_3d = false) without erroring.
            n_vox = 200;
            entry = testCase.makeFilledEntry(n_vox);
            entry.adc_vector   = zeros(n_vox, 1);
            entry.d_vector     = zeros(n_vox, 1);
            entry.f_vector     = zeros(n_vox, 1);
            entry.dstar_vector = zeros(n_vox, 1);
            entry.dose_vector  = 60 * ones(n_vox, 1);

            cfg = struct( ...
                'adc_thresh', 1.0, 'd_thresh', 1.0, ...
                'f_thresh', 1.0, 'dstar_thresh', 1.0, ...
                'dwi_types_to_run', 1);
            dv(1, 1, 1) = entry;
            gtv = {'/nonexistent/path/mask.mat'};

            [d95a, ~, ~, ~, ~, ~, ~, ~] = ...
                metrics_dosimetry({'Pt01'}, {'Pt01'}, 1, cfg, dv, gtv);

            testCase.verifyFalse(isnan(d95a(1, 1)), ...
                'Absent GTV mat file should fall back to 1D path; output should still be computed.');
        end

        function testNTpCellWrappingHandled(testCase)
            % The function accepts nTp as a scalar cell array {nTp_val}
            % in addition to a plain scalar.
            [m_id, id, nTp, cfg, dv, gtv] = ...
                testCase.buildEmptyVectorInputs(2, 2, 2);

            % Wrap nTp in a cell as the code says: `if iscell(nTp), nTp = nTp{1}; end`
            [d95a, ~, ~, ~, ~, ~, ~, ~] = ...
                metrics_dosimetry(m_id, id, {nTp}, cfg, dv, gtv);

            testCase.verifyEqual(size(d95a), [2, nTp], ...
                'Cell-wrapped nTp should be unwrapped to a scalar before use.');
        end

    end
end
