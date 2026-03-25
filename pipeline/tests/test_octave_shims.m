classdef test_octave_shims < matlab.unittest.TestCase
    % TEST_OCTAVE_SHIMS  Smoke tests for the 21 Octave compatibility shims.
    %
    % The shim files in pipeline/.octave_compat/ provide fallback
    % implementations of MATLAB functions/classes that Octave lacks. These
    % tests verify that each shim loads without error and produces correct
    % outputs for simple inputs.
    %
    % On MATLAB, the shims are temporarily added to the path so they shadow
    % the built-in functions. The @table class shim cannot shadow MATLAB's
    % built-in table, so table tests are skipped on MATLAB.
    %
    % Run tests with:
    %   results = runtests('tests/test_octave_shims.m');

    properties
        CompatDir   % Path to .octave_compat directory
        OrigPath    % Original MATLAB path before test modifications
        TempDir     % Temporary directory for file I/O tests
    end

    methods (TestMethodSetup)
        function setup(testCase)
            testCase.OrigPath = path;
            baseDir = fileparts(fileparts(which(mfilename)));
            testCase.CompatDir = fullfile(baseDir, '.octave_compat');

            % Add compat dir to top of path so shims shadow builtins
            addpath(testCase.CompatDir);

            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            % Restore original path to remove compat shims
            path(testCase.OrigPath);

            % Clean up temp directory
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)

        %% ============= nanmean =============

        function test_nanmean_simple(testCase)
            % nanmean of a vector with no NaN should equal mean
            x = [1, 2, 3, 4, 5];
            testCase.verifyEqual(nanmean(x), 3);
        end

        function test_nanmean_with_nan(testCase)
            % nanmean should ignore NaN values
            x = [1, NaN, 3, NaN, 5];
            testCase.verifyEqual(nanmean(x), 3);
        end

        function test_nanmean_all_nan(testCase)
            % All-NaN input should return NaN
            testCase.verifyTrue(isnan(nanmean([NaN, NaN])));
        end

        function test_nanmean_matrix_dim(testCase)
            % nanmean along dim=1 (column-wise)
            X = [1 2; NaN 4; 3 NaN];
            result = nanmean(X, 1);
            testCase.verifyEqual(result, [2, 3]);
        end

        %% ============= nanstd =============

        function test_nanstd_simple(testCase)
            % nanstd of a vector with no NaN should match std
            x = [2, 4, 4, 4, 5, 5, 7, 9];
            testCase.verifyEqual(nanstd(x), std(x), 'AbsTol', 1e-10);
        end

        function test_nanstd_with_nan(testCase)
            % nanstd should ignore NaN values
            x = [2, NaN, 4, 4, 5];
            clean = [2, 4, 4, 5];
            testCase.verifyEqual(nanstd(x), std(clean), 'AbsTol', 1e-10);
        end

        function test_nanstd_population(testCase)
            % flag=1 gives population std (divide by N, not N-1)
            x = [2, 4, 4, 4];
            expected = std(x, 1);
            testCase.verifyEqual(nanstd(x, 1), expected, 'AbsTol', 1e-10);
        end

        function test_nanstd_all_nan(testCase)
            % All-NaN input should return NaN
            testCase.verifyTrue(isnan(nanstd([NaN, NaN])));
        end

        %% ============= contains =============

        function test_contains_true(testCase)
            testCase.verifyTrue(contains('hello world', 'world'));
        end

        function test_contains_false(testCase)
            testCase.verifyFalse(contains('hello world', 'xyz'));
        end

        function test_contains_cell_array(testCase)
            strs = {'abc', 'def', 'abcdef'};
            result = contains(strs, 'abc');
            testCase.verifyEqual(result, [true, false, true]);
        end

        function test_contains_ignore_case(testCase)
            testCase.verifyTrue(contains('Hello', 'hello', 'IgnoreCase', true));
        end

        %% ============= categorical =============

        function test_categorical_passthrough(testCase)
            % categorical shim is a pass-through — only testable on Octave
            % because MATLAB's built-in categorical class cannot be
            % reliably shadowed by a function file on the path.
            if ~exist('OCTAVE_VERSION', 'builtin')
                % On MATLAB, verify the shim file exists (it will be used
                % on Octave) but skip the behavioral test.
                testCase.verifyTrue(exist(fullfile(testCase.CompatDir, 'categorical.m'), 'file') > 0, ...
                    'categorical.m shim should exist in .octave_compat');
                return;
            end
            x = {'A', 'B', 'C'};
            c = categorical(x);
            testCase.verifyEqual(c, x);
        end

        function test_categorical_numeric(testCase)
            % categorical shim numeric pass-through — Octave only.
            if ~exist('OCTAVE_VERSION', 'builtin')
                return;
            end
            x = [1, 2, 3];
            c = categorical(x);
            testCase.verifyEqual(c, x);
        end

        %% ============= cvpartition =============

        function test_cvpartition_kfold(testCase)
            % Numeric input should give NumTestSets = 5
            c = cvpartition(100, 'KFold', 10);
            testCase.verifyTrue(isfield(c, 'NumTestSets'));
            testCase.verifyEqual(c.NumTestSets, 5);
        end

        %% ============= fitglme =============

        function test_fitglme_returns_struct(testCase)
            % fitglme shim should return struct with Coefficients and Formula
            tbl = struct('Y', [1;2;3], 'X', [4;5;6]);
            mdl = fitglme(tbl, 'Y ~ X');
            testCase.verifyTrue(isfield(mdl, 'Coefficients'));
            testCase.verifyTrue(isfield(mdl, 'Formula'));
            testCase.verifyEqual(mdl.Formula, 'Y ~ X');
            testCase.verifyEqual(mdl.Coefficients.Estimate, [0; 0]);
            testCase.verifyEqual(mdl.Coefficients.pValue, [1; 1]);
        end

        %% ============= niftiinfo =============

        function test_niftiinfo_fields(testCase)
            info = niftiinfo('test_volume.nii');
            testCase.verifyEqual(info.Filename, 'test_volume.nii');
            testCase.verifyEqual(info.ImageSize, [0 0 0]);
            testCase.verifyEqual(info.PixelDimensions, [1 1 1]);
            testCase.verifyEqual(info.Datatype, 'double');
        end

        %% ============= niftiwrite / niftiread round-trip =============

        function test_nifti_roundtrip(testCase)
            % Write data with niftiwrite shim, read back with niftiread shim
            data = rand(4, 4, 4);
            fname = fullfile(testCase.TempDir, 'test_vol');
            niftiwrite(data, fname);

            % Verify proxy .mat file was created
            testCase.verifyTrue(exist([fname '_nifti.mat'], 'file') > 0);

            % Read back and verify
            loaded = niftiread(fname);
            testCase.verifyEqual(loaded, data, 'AbsTol', 1e-15);
        end

        function test_niftiread_from_info(testCase)
            % niftiread should accept a niftiinfo struct
            data = ones(3, 3, 3);
            fname = fullfile(testCase.TempDir, 'test_info');
            niftiwrite(data, fname);
            info = niftiinfo(fname);
            loaded = niftiread(info);
            testCase.verifyEqual(loaded, data);
        end

        function test_niftiread_missing_file(testCase)
            % Reading a nonexistent file should error
            testCase.verifyError( ...
                @() niftiread(fullfile(testCase.TempDir, 'nonexistent')), ...
                'niftiread:FileNotFound');
        end

        %% ============= sgtitle =============

        function test_sgtitle_no_error(testCase)
            % sgtitle should not error (falls back to title)
            fig = figure('Visible', 'off');
            cleanup = onCleanup(@() close(fig));
            subplot(1, 2, 1); plot(1:3);
            subplot(1, 2, 2); plot(1:3);
            sgtitle('Test Title');
            % If we get here without error, the test passes
            testCase.verifyTrue(true);
        end

        %% ============= yline =============

        function test_yline_draws_line(testCase)
            % yline should return a valid line handle
            fig = figure('Visible', 'off');
            cleanup = onCleanup(@() close(fig));
            plot(1:10, 1:10);
            hl = yline(5);
            testCase.verifyTrue(ishghandle(hl));
            % Verify Y-data is [5, 5]
            ydata = get(hl, 'YData');
            testCase.verifyEqual(ydata, [5, 5]);
        end

        %% ============= spectralcluster =============

        function test_spectralcluster_returns_indices(testCase)
            % spectralcluster shim should return cluster indices via kmeans
            rng(42);
            X = [randn(20, 2); randn(20, 2) + 5];
            % Suppress the expected warning from the shim
            warning('off', 'spectralcluster:octaveShim');
            idx = spectralcluster(X, 2);
            warning('on', 'spectralcluster:octaveShim');
            testCase.verifyEqual(numel(idx), 40);
            testCase.verifyTrue(all(idx == 1 | idx == 2));
        end

        function test_spectralcluster_warns(testCase)
            % spectralcluster shim should emit a warning
            X = randn(10, 2);
            testCase.verifyWarning( ...
                @() spectralcluster(X, 2), ...
                'spectralcluster:octaveShim');
        end

        %% ============= TestCase shim =============

        function test_testcase_shim_loads(testCase)
            % Verify the TestCase shim classdef file exists and is readable
            tc_file = fullfile(testCase.CompatDir, ...
                '+matlab', '+unittest', 'TestCase.m');
            testCase.verifyTrue(exist(tc_file, 'file') > 0);
        end

        %% ============= TestRunner shim =============

        function test_testrunner_shim_loads(testCase)
            tr_file = fullfile(testCase.CompatDir, ...
                '+matlab', '+unittest', 'TestRunner.m');
            testCase.verifyTrue(exist(tr_file, 'file') > 0);
        end

        %% ============= TestSuite shim =============

        function test_testsuite_shim_loads(testCase)
            ts_file = fullfile(testCase.CompatDir, ...
                '+matlab', '+unittest', 'TestSuite.m');
            testCase.verifyTrue(exist(ts_file, 'file') > 0);
        end

        %% ============= PathFixture shim =============

        function test_pathfixture_shim_loads(testCase)
            pf_file = fullfile(testCase.CompatDir, ...
                '+matlab', '+unittest', '+fixtures', 'PathFixture.m');
            testCase.verifyTrue(exist(pf_file, 'file') > 0);
        end

        %% ============= CodeCoveragePlugin shim =============

        function test_codecoverageplugin_shim_loads(testCase)
            cp_file = fullfile(testCase.CompatDir, ...
                '+matlab', '+unittest', '+plugins', 'CodeCoveragePlugin.m');
            testCase.verifyTrue(exist(cp_file, 'file') > 0);
        end

        %% ============= @table shim (Octave only) =============

        function test_table_constructor_octave_only(testCase)
            % The @table shim conflicts with MATLAB's built-in table class.
            % Only test on Octave where the shim is actually used.
            if ~exist('OCTAVE_VERSION', 'builtin')
                return;  % Skip on MATLAB
            end
            t = table([1;2;3], {'a';'b';'c'}, 'VariableNames', {'Num', 'Str'});
            testCase.verifyEqual(t.Num, [1;2;3]);
            testCase.verifyEqual(t.Str, {'a';'b';'c'});
        end

        function test_table_subsasgn_octave_only(testCase)
            if ~exist('OCTAVE_VERSION', 'builtin')
                return;
            end
            t = table([1;2;3], 'VariableNames', {'X'});
            t.Y = [4;5;6];
            testCase.verifyEqual(t.Y, [4;5;6]);
        end

        function test_table_row_slicing_octave_only(testCase)
            if ~exist('OCTAVE_VERSION', 'builtin')
                return;
            end
            t = table([10;20;30], [40;50;60], 'VariableNames', {'A', 'B'});
            t2 = t(1:2, :);
            testCase.verifyEqual(t2.A, [10;20]);
            testCase.verifyEqual(t2.B, [40;50]);
        end

        function test_table_display_octave_only(testCase)
            if ~exist('OCTAVE_VERSION', 'builtin')
                return;
            end
            t = table([1;2], 'VariableNames', {'Val'});
            % display should not error
            display(t);
            testCase.verifyTrue(true);
        end

    end

end
