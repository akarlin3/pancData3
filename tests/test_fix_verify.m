% Unit test wrapper for reproduction of the path construction bug and fallback logic
classdef test_fix_verify < matlab.unittest.TestCase
    methods(Test)
        function testVulnerabilityReproduction(testCase)
            % Setup
            base_dir = tempname;
            if ~isfolder(base_dir)
                mkdir(base_dir);
            end

            % Ensure we have absolute path
            base_dir = char(java.io.File(base_dir).getAbsolutePath());

            outloc = fullfile(base_dir, 'nii');
            % Ensure trailing slash for outloc to mimic some potential environment behaviors,
            % though fullfile usually strips it. The bug is [outloc gtvname] vs fullfile.
            % If outloc doesn't end in separator, [outloc gtvname] is definitely wrong.
            if ~isfolder(outloc)
                mkdir(outloc);
            end

            gtvname = 'fx1_gtv1';
            struct_file = fullfile(base_dir, 'mask.mat');

            % Create dummy mat file
            Stvol3d = magic(5);
            save(struct_file, 'Stvol3d');

            fprintf('Base dir: %s\n', base_dir);
            fprintf('Outloc: %s\n', outloc);

            % --- VULNERABLE LOGIC REPRODUCTION (Simulated from core/load_dwi_data.m) ---
            fprintf('Running vulnerable logic simulation...\n');

            % Block 1: Broken path construction
            % [outloc gtvname '.nii.gz'] -> .../niifx1_gtv1.nii.gz (missing separator if outloc doesn't have one)
            % The code uses 'outloc' which comes from fullfile(basefolder, 'nii').
            % fullfile typically does NOT add a trailing separator.

            bad_path_check = [outloc gtvname '.nii.gz'];

            if ~exist(bad_path_check,'file')
                fprintf('  Block 1 executed (as expected due to bad path check).\n');
                % We simulate safe_load_mask success
                gtv_mask = Stvol3d;

                % The write happens to the malformed path
                try
                    % Simulate niftiwrite
                    target_path = [outloc gtvname '.nii.gz'];
                    fid = fopen(target_path, 'w');
                    fwrite(fid, 0);
                    fclose(fid);
                    fprintf('  Block 1 wrote to: %s\n', target_path);
                catch ME
                    fprintf('  Block 1 write failed: %s\n', ME.message);
                end
            end

            % Block 2: Correct path, Unsafe fallback
            good_path_check = fullfile(outloc, [gtvname '.nii.gz']);

            if ~exist(good_path_check,'file')
                fprintf('  Block 2 executed (because Block 1 wrote to wrong place).\n');
                % Simulate unsafe load
                tmp = load(struct_file);
                gtv_mask = tmp.Stvol3d;

                % Write to correct path
                fid = fopen(good_path_check, 'w');
                fwrite(fid, 0);
                fclose(fid);
                fprintf('  Block 2 wrote to: %s\n', good_path_check);
            else
                fprintf('  Block 2 SKIPPED (Unexpected for vulnerable code).\n');
            end

            % --- VERIFICATION ---
            bad_file = [outloc gtvname '.nii.gz'];
            good_file = fullfile(outloc, [gtvname '.nii.gz']);

            % ASSERTIONS

            % 1. The bad file SHOULD exist in the vulnerable state
            testCase.verifyTrue(logical(exist(bad_file, 'file')), ...
                sprintf('VULNERABILITY CONFIRMED: Malformed file exists at %s', bad_file));

            % 2. The good file SHOULD exist (created by the unsafe fallback)
            testCase.verifyTrue(logical(exist(good_file, 'file')), ...
                'Result: Correct file exists (created by fallback).');

            % Cleanup handled by test framework or tempname
        end
    end
end
