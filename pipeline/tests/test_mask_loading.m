% TEST_MASK_LOADING  Diagnostic spot-check for 3D GTV mask consistency.
%
% This script (not a unittest) verifies that the voxel count stored in
% data_vectors_gtvp matches the number of nonzero voxels in the
% corresponding 3D GTV mask file (Stvol3d).  A mismatch indicates that
% the mask and the extracted feature vectors have diverged, which would
% invalidate all downstream per-voxel analyses.
%
% The script requires local patient data (pipeline_voxels.mat) and the
% original GTV .mat mask files on disk.  If neither is available it
% exits gracefully with a skip message.
%
% Output is captured to test_mask_loading_output.txt via MATLAB diary.

test_output_dir = fullfile(tempdir, 'test_mask_loading_output');
if ~exist(test_output_dir, 'dir'), mkdir(test_output_dir); end
diary(fullfile(test_output_dir, 'test_mask_loading_output.txt'));
diary on;

% Locate the data directory from config.json (if present) or fall back
% to the repository root.
config_file = fullfile(fileparts(mfilename('fullpath')), '..', 'config.json');
if isfile(config_file)
    config_struct = parse_config(config_file);
    dataloc = config_struct.dataloc;
else
    dataloc = fullfile(fileparts(mfilename('fullpath')), '..', filesep);
end

% Load the pre-computed DWI data vectors and GTV file paths.
data_file = fullfile(dataloc, 'pipeline_voxels.mat');
if ~exist(data_file, 'file')
    fprintf('Skipping test_mask_loading: no local data_vectors.mat found.\n');
    diary off;
    return;
end

load(data_file, 'id_list', 'data_vectors_gtvp', 'gtv_locations');

% Iterate over patients and fractions to find the first available GTV mask
% file on disk, then compare the mask voxel count against the stored
% feature vector length.
found = false;
for j = 1:size(gtv_locations,1)
    for k = 1:size(gtv_locations,2)
        gtv_mat = gtv_locations{j,k,1};
        if ~isempty(gtv_mat)
            % Normalize path separators for the current OS (the stored
            % paths may use Windows or Unix separators).
            path_parts = strsplit(gtv_mat, {'/', '\'});
            gtv_mat = fullfile(path_parts{:});
            % Restore leading '/' on Unix if strsplit consumed it
            if isunix && (startsWith(gtv_mat, filesep) == 0) && isempty(path_parts{1})
                gtv_mat = [filesep gtv_mat];
            end

            if exist(gtv_mat, 'file')
                % Load the 3D binary mask into a struct to avoid workspace
                % injection (security best practice for untrusted .mat files).
                loaded_data = load(gtv_mat, 'Stvol3d');
                gtv_mask = loaded_data.Stvol3d;

                % Compare: vector length should equal the number of
                % mask-positive voxels.
                vec_len = length(data_vectors_gtvp(j,k,1).adc_vector);
                mask_sum = sum(gtv_mask(:) == 1);

                fprintf('Patient %d, Fraction %d\n', j, k);
                fprintf('Vector length: %d\n', vec_len);
                fprintf('Mask sum: %d\n', mask_sum);
                if vec_len == mask_sum
                    fprintf('MATCH!\n');
                else
                    fprintf('MISMATCH!\n');
                end
                found = true;
                break;
            end
        end
    end
    if found, break; end
end

diary off;
