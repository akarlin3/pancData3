% Test 3D mask loading
dataloc = '\\pensmph6\mpcsresearch1\aliottae\pancreas_dwi\';
load([dataloc 'adc_vectors.mat'], 'id_list', 'data_vectors_gtvp', 'gtv_locations');

found = false;
for j = 1:size(gtv_locations,1)
    for k = 1:size(gtv_locations,2)
        gtv_mat = gtv_locations{j,k,1};
        if ~isempty(gtv_mat)
            % convert MacOS path to Windows dataloc
            gtv_mat = strrep(gtv_mat, '/Volumes/aliottae/pancreas_dwi/', dataloc);
            % replace forward slashes with backward slashes
            gtv_mat = strrep(gtv_mat, '/', '\');
            if exist(gtv_mat, 'file')
                load(gtv_mat, 'Stvol3d');
                gtv_mask = Stvol3d;
                
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
