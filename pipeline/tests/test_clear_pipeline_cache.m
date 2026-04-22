classdef test_clear_pipeline_cache < matlab.unittest.TestCase
    % TEST_CLEAR_PIPELINE_CACHE Unit tests for clear_pipeline_cache.
    %
    % Validates cache clearing: deletion of derived .mat files
    % (summary_metrics*, adc_vectors), full protection of every
    % dwi_vectors*.mat file, sentinel-verified checkpoint removal,
    % once-per-session guard via persistent variable, and no-op when
    % clear_cache is false.

    properties
        TmpDir
    end

    methods(TestMethodSetup)
        function setupTempDir(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            % Reset the persistent variable by clearing the function
            clear clear_pipeline_cache;
        end
    end

    methods(TestMethodTeardown)
        function cleanupTempDir(testCase)
            % Reset persistent state after each test
            clear clear_pipeline_cache;
            if isfolder(testCase.TmpDir)
                rmdir(testCase.TmpDir, 's');
            end
        end
    end

    methods(Test)
        function test_deletes_derived_cache_files(testCase)
            % Derived caches (summary_metrics*, adc_vectors) should be
            % deleted, but dwi_vectors*.mat must always be preserved.
            f_dwi = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            f_sum = fullfile(testCase.TmpDir, 'summary_metrics_Standard.mat');
            f_adc = fullfile(testCase.TmpDir, 'adc_vectors.mat');
            fclose(fopen(f_dwi, 'w'));
            fclose(fopen(f_sum, 'w'));
            fclose(fopen(f_adc, 'w'));

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);

            testCase.verifyTrue(cleared);
            testCase.verifyTrue(exist(f_dwi, 'file') == 2, ...
                'dwi_vectors*.mat must never be deleted by clear_pipeline_cache');
            testCase.verifyFalse(exist(f_sum, 'file') == 2);
            testCase.verifyFalse(exist(f_adc, 'file') == 2);
        end

        function test_protects_all_dwi_vectors_files(testCase)
            % Every dwi_vectors*.mat variant must survive clear_cache.
            variants = {'dwi_vectors.mat', 'dwi_vectors_Standard.mat', ...
                'dwi_vectors_dnCNN.mat', 'dwi_vectors_IVIMnet.mat', ...
                'dwi_vectors_ea.mat', 'dwi_vectors_2026_Apr_22.mat'};
            for vi = 1:numel(variants)
                fclose(fopen(fullfile(testCase.TmpDir, variants{vi}), 'w'));
            end

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            clear_pipeline_cache(cfg);

            for vi = 1:numel(variants)
                p = fullfile(testCase.TmpDir, variants{vi});
                testCase.verifyTrue(exist(p, 'file') == 2, ...
                    sprintf('%s was deleted but must be preserved', variants{vi}));
            end
        end

        function test_removes_checkpoint_with_sentinel(testCase)
            % Checkpoint dir with .pipeline_created sentinel should be removed.
            cp_dir = fullfile(testCase.TmpDir, 'processed_patients');
            mkdir(cp_dir);
            fid = fopen(fullfile(cp_dir, '.pipeline_created'), 'w');
            fprintf(fid, 'test sentinel');
            fclose(fid);

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            clear_pipeline_cache(cfg);

            testCase.verifyFalse(isfolder(cp_dir));
        end

        function test_skips_checkpoint_without_sentinel(testCase)
            % Checkpoint dir without sentinel should be preserved.
            cp_dir = fullfile(testCase.TmpDir, 'processed_patients');
            mkdir(cp_dir);

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            clear_pipeline_cache(cfg);

            testCase.verifyTrue(isfolder(cp_dir));
        end

        function test_noop_when_clear_cache_false(testCase)
            % When clear_cache is false, nothing should be deleted.
            % Uses summary_metrics (deletable) so the assertion would
            % actually fail if the early-return path were broken.
            f1 = fullfile(testCase.TmpDir, 'summary_metrics_Standard.mat');
            fclose(fopen(f1, 'w'));

            cfg = struct('clear_cache', false, 'dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);

            testCase.verifyFalse(cleared);
            testCase.verifyTrue(exist(f1, 'file') == 2);
        end

        function test_noop_when_field_missing(testCase)
            % When clear_cache field is absent, nothing should be deleted.
            f1 = fullfile(testCase.TmpDir, 'summary_metrics_Standard.mat');
            fclose(fopen(f1, 'w'));

            cfg = struct('dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);

            testCase.verifyFalse(cleared);
            testCase.verifyTrue(exist(f1, 'file') == 2);
        end

        function test_once_per_session_guard(testCase)
            % Second call should be a no-op (returns false).  Use a
            % deletable file (summary_metrics_*) to exercise the guard:
            % dwi_vectors*.mat is protected regardless of the guard.
            f1 = fullfile(testCase.TmpDir, 'summary_metrics_Standard.mat');
            fclose(fopen(f1, 'w'));

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            cleared1 = clear_pipeline_cache(cfg);
            testCase.verifyTrue(cleared1);

            % Create another file and call again — should NOT be deleted
            f2 = fullfile(testCase.TmpDir, 'summary_metrics_dnCNN.mat');
            fclose(fopen(f2, 'w'));
            cleared2 = clear_pipeline_cache(cfg);
            testCase.verifyFalse(cleared2);
            testCase.verifyTrue(exist(f2, 'file') == 2);
        end

        function test_no_files_to_clear(testCase)
            % Empty directory should still succeed without error.
            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);
            testCase.verifyTrue(cleared);
        end
    end
end
