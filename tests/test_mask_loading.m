% Test 3D mask loading
if ~exist('saved_figures', 'dir'), mkdir('saved_figures'); end
diary(fullfile('saved_figures', 'test_mask_loading_output.txt'));
diary on;

config_file = fullfile(fileparts(mfilename('fullpath')), '..', 'config.json');
if isfile(config_file)
    config_struct = parse_config(config_file);
    dataloc = config_struct.dataloc;
else
    dataloc = fullfile(fileparts(mfilename('fullpath')), '..', filesep);
end

tmp_dwi = load(fullfile(dataloc, 'dwi_vectors.mat'), 'id_list', 'data_vectors_gtvp', 'gtv_locations');
id_list = tmp_dwi.id_list;
data_vectors_gtvp = tmp_dwi.data_vectors_gtvp;
gtv_locations = tmp_dwi.gtv_locations;

found = false;
for j = 1:size(gtv_locations,1)
    for k = 1:size(gtv_locations,2)
        gtv_mat = gtv_locations{j,k,1};
        if ~isempty(gtv_mat)
            % Ensure correct file separator for current platform
            path_parts = strsplit(gtv_mat, {'/', '\'});
            gtv_mat = fullfile(path_parts{:});
            % Handle potential leading separator (absolute path on Unix) lost by strsplit
            if isunix && (startsWith(gtv_mat, filesep) == 0) && isempty(path_parts{1})
                gtv_mat = [filesep gtv_mat];
            end

            if exist(gtv_mat, 'file')
                tmp_gtv = load(gtv_mat, 'Stvol3d');
                gtv_mask = tmp_gtv.Stvol3d;
                
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
