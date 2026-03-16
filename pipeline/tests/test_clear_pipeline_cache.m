classdef test_clear_pipeline_cache < matlab.unittest.TestCase
    % TEST_CLEAR_PIPELINE_CACHE Unit tests for clear_pipeline_cache.
    %
    % Validates cache clearing: deletion of pipeline-generated .mat files,
    % protection of curated files, sentinel-verified checkpoint removal,
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
        function test_deletes_cache_files(testCase)
            % Pipeline-generated .mat files should be deleted.
            f1 = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            f2 = fullfile(testCase.TmpDir, 'summary_metrics_Standard.mat');
            f3 = fullfile(testCase.TmpDir, 'adc_vectors.mat');
            fclose(fopen(f1, 'w'));
            fclose(fopen(f2, 'w'));
            fclose(fopen(f3, 'w'));

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);

            testCase.verifyTrue(cleared);
            testCase.verifyFalse(exist(f1, 'file') == 2);
            testCase.verifyFalse(exist(f2, 'file') == 2);
            testCase.verifyFalse(exist(f3, 'file') == 2);
        end

        function test_protects_curated_files(testCase)
            % dwi_vectors_ea.mat should NOT be deleted.
            protected = fullfile(testCase.TmpDir, 'dwi_vectors_ea.mat');
            fclose(fopen(protected, 'w'));

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            clear_pipeline_cache(cfg);

            testCase.verifyTrue(exist(protected, 'file') == 2);
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
            f1 = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            fclose(fopen(f1, 'w'));

            cfg = struct('clear_cache', false, 'dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);

            testCase.verifyFalse(cleared);
            testCase.verifyTrue(exist(f1, 'file') == 2);
        end

        function test_noop_when_field_missing(testCase)
            % When clear_cache field is absent, nothing should be deleted.
            f1 = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            fclose(fopen(f1, 'w'));

            cfg = struct('dataloc', testCase.TmpDir);
            cleared = clear_pipeline_cache(cfg);

            testCase.verifyFalse(cleared);
            testCase.verifyTrue(exist(f1, 'file') == 2);
        end

        function test_once_per_session_guard(testCase)
            % Second call should be a no-op (returns false).
            f1 = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            fclose(fopen(f1, 'w'));

            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            cleared1 = clear_pipeline_cache(cfg);
            testCase.verifyTrue(cleared1);

            % Create another file and call again — should NOT be deleted
            f2 = fullfile(testCase.TmpDir, 'dwi_vectors_dnCNN.mat');
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
