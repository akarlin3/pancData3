classdef test_backfill_checkpoint_sentinel < matlab.unittest.TestCase
% TEST_BACKFILL_CHECKPOINT_SENTINEL — Unit tests for backfill_checkpoint_sentinel.m
%
% Validates that the sentinel is written only when the directory looks
% pipeline-owned (contains pipeline-style files and nothing foreign),
% and that all no-op paths (missing dir, existing sentinel, empty dir,
% foreign files) leave the directory untouched.

    properties
        TmpDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            testCase.OriginalPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            if exist(testCase.TmpDir, 'dir')
                rmdir(testCase.TmpDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods(Access = private)
        function touch(~, p)
            fid = fopen(p, 'w'); fclose(fid);
        end

        function p = checkpoint_subdir(testCase)
            p = fullfile(testCase.TmpDir, 'processed_patients');
            mkdir(p);
        end
    end

    methods(Test)
        function test_backfills_when_pipeline_files_present(testCase)
            % Directory contains valid patient_NNN_*.mat files, no sentinel.
            cp = testCase.checkpoint_subdir();
            testCase.touch(fullfile(cp, 'patient_001_P01.mat'));
            testCase.touch(fullfile(cp, 'patient_002_P02.mat'));

            backfilled = backfill_checkpoint_sentinel(cp);

            testCase.verifyTrue(backfilled);
            testCase.verifyTrue(exist(fullfile(cp, '.pipeline_created'), 'file') == 2);
        end

        function test_backfills_with_lock_file_present(testCase)
            % Lock files also count as pipeline-owned; sentinel still backfills.
            cp = testCase.checkpoint_subdir();
            testCase.touch(fullfile(cp, 'patient_001_P01.mat'));
            testCase.touch(fullfile(cp, 'patient_002_P02.lock'));

            backfilled = backfill_checkpoint_sentinel(cp);

            testCase.verifyTrue(backfilled);
            testCase.verifyTrue(exist(fullfile(cp, '.pipeline_created'), 'file') == 2);
        end

        function test_skips_when_sentinel_already_exists(testCase)
            % Pre-existing sentinel: function should be a no-op.
            cp = testCase.checkpoint_subdir();
            testCase.touch(fullfile(cp, 'patient_001_P01.mat'));
            sentinel_path = fullfile(cp, '.pipeline_created');
            fid = fopen(sentinel_path, 'w');
            fprintf(fid, 'pre-existing'); fclose(fid);

            backfilled = backfill_checkpoint_sentinel(cp);

            testCase.verifyFalse(backfilled);
            % And the original content is preserved
            content = fileread(sentinel_path);
            testCase.verifySubstring(content, 'pre-existing');
        end

        function test_skips_empty_directory(testCase)
            % Empty dir: no evidence it's pipeline-owned, do not backfill.
            cp = testCase.checkpoint_subdir();

            backfilled = backfill_checkpoint_sentinel(cp);

            testCase.verifyFalse(backfilled);
            testCase.verifyFalse(exist(fullfile(cp, '.pipeline_created'), 'file') == 2);
        end

        function test_skips_when_foreign_file_present(testCase)
            % Foreign file mixed in: refuse to claim the directory.
            cp = testCase.checkpoint_subdir();
            testCase.touch(fullfile(cp, 'patient_001_P01.mat'));
            testCase.touch(fullfile(cp, 'my_notes.txt'));

            backfilled = backfill_checkpoint_sentinel(cp);

            testCase.verifyFalse(backfilled);
            testCase.verifyFalse(exist(fullfile(cp, '.pipeline_created'), 'file') == 2);
        end

        function test_skips_when_only_lock_files(testCase)
            % No patient_*.mat means no completed checkpoints — refuse.
            cp = testCase.checkpoint_subdir();
            testCase.touch(fullfile(cp, 'patient_001_P01.lock'));

            backfilled = backfill_checkpoint_sentinel(cp);

            testCase.verifyFalse(backfilled);
            testCase.verifyFalse(exist(fullfile(cp, '.pipeline_created'), 'file') == 2);
        end

        function test_skips_when_directory_missing(testCase)
            % No directory at all — clean no-op, no error.
            missing = fullfile(testCase.TmpDir, 'does_not_exist');

            backfilled = backfill_checkpoint_sentinel(missing);

            testCase.verifyFalse(backfilled);
        end

        function test_clear_cache_succeeds_after_backfill(testCase)
            % End-to-end: backfill, then clear_pipeline_cache should now
            % sweep the directory because the sentinel is present.
            cp = testCase.checkpoint_subdir();
            testCase.touch(fullfile(cp, 'patient_001_P01.mat'));
            backfill_checkpoint_sentinel(cp);

            % Reset clear_pipeline_cache's persistent guard
            clear clear_pipeline_cache;
            cfg = struct('clear_cache', true, 'dataloc', testCase.TmpDir);
            clear_pipeline_cache(cfg);
            clear clear_pipeline_cache;

            testCase.verifyFalse(isfolder(cp), ...
                'processed_patients should be cleared after sentinel backfill');
        end
    end
end
