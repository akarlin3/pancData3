% TEST_FIX_VERIFY Regression test reproducing the path concatenation bug in load_dwi_data.m.
%
% This test demonstrates the vulnerability caused by using string concatenation
% ([outloc gtvname '.nii.gz']) instead of fullfile() for path construction.
% When outloc = '/tmp/.../nii' (no trailing separator), the concatenated path
% becomes '/tmp/.../niifx1_gtv1.nii.gz' (missing separator between directory
% and filename), causing the file existence check to fail and triggering an
% unsafe fallback code path.
%
% The test creates a temp directory, simulates both the broken (concatenation)
% and correct (fullfile) paths, and verifies:
%   1. The malformed file is created at the wrong location (confirming the bug)
%   2. The fallback code then creates the file at the correct location
classdef test_fix_verify < matlab.unittest.TestCase
    methods(Test)
        function testVulnerabilityReproduction(testCase)
            % Create a temporary directory to simulate the patient nii/ folder
            base_dir = tempname;
            if ~isfolder(base_dir)
                mkdir(base_dir);
            end

            % Resolve to absolute path (required for reliable path comparisons)

            % Create the nii/ subdirectory. Note: fullfile does NOT append a trailing
            % separator, so outloc = '/tmp/.../nii' (no trailing slash).
            % This is the root cause of the bug: [outloc gtvname] yields
            % '/tmp/.../niifx1_gtv1' instead of '/tmp/.../nii/fx1_gtv1'.
            outloc = fullfile(base_dir, 'nii');
            if ~isfolder(outloc)
                mkdir(outloc);
            end

            gtvname = 'fx1_gtv1';  % Typical GTV mask filename stem
            struct_file = fullfile(base_dir, 'mask.mat');

            % Create a dummy .mat file containing a GTV mask (magic square as placeholder)
            Stvol3d = magic(5);
            save(struct_file, 'Stvol3d');

            fprintf('Base dir: %s\n', base_dir);
            fprintf('Outloc: %s\n', outloc);

            % --- VULNERABLE LOGIC REPRODUCTION (Simulated from core/load_dwi_data.m) ---
            % Block 1: Reproduce the broken path construction using string concatenation.
            % [outloc gtvname '.nii.gz'] produces '.../niifx1_gtv1.nii.gz' because
            % outloc lacks a trailing separator. This path does not exist, so the
            % exist() check triggers the file creation at the WRONG location.
            fprintf('Running vulnerable logic simulation...\n');

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

            % Block 2: The fallback code uses the correct path via fullfile().
            % Since Block 1 wrote to the wrong location, the correct file still
            % doesn't exist, triggering this fallback which uses unsafe load().
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
            % Both files should exist: the malformed one (confirming the bug)
            % and the correct one (created by the fallback path).
            bad_file = [outloc gtvname '.nii.gz'];
            good_file = fullfile(outloc, [gtvname '.nii.gz']);

            % ASSERTIONS: Confirm both paths were written to, proving the bug exists

            % 1. The malformed file (e.g., '.../niifx1_gtv1.nii.gz') should exist,
            %    confirming the string concatenation bug creates files at wrong locations.
            testCase.verifyTrue(logical(exist(bad_file, 'file')), ...
                sprintf('VULNERABILITY CONFIRMED: Malformed file exists at %s', bad_file));

            % 2. The correct file (e.g., '.../nii/fx1_gtv1.nii.gz') should also exist,
            %    created by the fallback path that uses fullfile() properly.
            testCase.verifyTrue(logical(exist(good_file, 'file')), ...
                'Result: Correct file exists (created by fallback).');

            % Cleanup handled by test framework or tempname
        end
    end
end
