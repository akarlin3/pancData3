function [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub] = metrics_dosimetry(m_id_list, id_list, nTp, config_struct, m_data_vectors_gtvp, gtv_locations)
% METRICS_DOSIMETRY — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 3/5 of the metrics step.

fprintf('  --- SECTION 6: Target Coverage (Sub-Volume Dose Metrics) ---\n');

adc_thresh   = config_struct.adc_thresh;    
d_thresh     = config_struct.d_thresh;      
f_thresh     = config_struct.f_thresh;      
dstar_thresh = config_struct.dstar_thresh;  

d95_adc_sub = nan(length(m_id_list), nTp);
v50_adc_sub = nan(length(m_id_list), nTp);

d95_d_sub = nan(length(m_id_list), nTp);
v50_d_sub = nan(length(m_id_list), nTp);

d95_f_sub = nan(length(m_id_list), nTp);
v50_f_sub = nan(length(m_id_list), nTp);

d95_dstar_sub = nan(length(m_id_list), nTp);
v50_dstar_sub = nan(length(m_id_list), nTp);

for j = 1:length(m_id_list)
    j_orig = find(strcmp(id_list, m_id_list{j}));
    for k = 1
        if isfield(config_struct, 'dwi_types_to_run') && isscalar(config_struct.dwi_types_to_run)
            dtype_idx = config_struct.dwi_types_to_run;
        else
            dtype_idx = 1;
        end
        
        switch dtype_idx
            case 1
                adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector;
                d_vec   = m_data_vectors_gtvp(j,k,1).d_vector;
                f_vec   = m_data_vectors_gtvp(j,k,1).f_vector;
                dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector;
            case 2
                adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector_dncnn;
                d_vec   = m_data_vectors_gtvp(j,k,1).d_vector_dncnn;
                f_vec   = m_data_vectors_gtvp(j,k,1).f_vector_dncnn;
                dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector_dncnn;
            case 3
                adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector;
                d_vec   = m_data_vectors_gtvp(j,k,1).d_vector_ivimnet;
                f_vec   = m_data_vectors_gtvp(j,k,1).f_vector_ivimnet;
                dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector_ivimnet;
        end
        dose_vec  = m_data_vectors_gtvp(j,k,1).dose_vector;
        
        if ~isempty(dose_vec) && ~isempty(adc_vec)
            gtv_mat = gtv_locations{j_orig, k, 1};
            has_3d = false;
            if ~isempty(gtv_mat)
                path_parts = strsplit(gtv_mat, {'/', '\'});
                gtv_mat = fullfile(path_parts{:});
                if isunix && (startsWith(gtv_mat, filesep) == 0) && isempty(path_parts{1})
                    gtv_mat = [filesep gtv_mat];
                end

                if exist(gtv_mat, 'file')
                    tmp = load(gtv_mat, 'Stvol3d');
                    gtv_mask_3d = tmp.Stvol3d;
                    if sum(gtv_mask_3d(:) == 1) == length(adc_vec)
                        has_3d = true;
                    end
                end
            end
            
            se = strel('sphere', 1);
            min_cc_voxels = 10;
            min_subvol_voxels = 100;
            
            adc_mask_1d = adc_vec < adc_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = adc_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                adc_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_adc_sub = dose_vec(adc_mask_1d);
            if ~isempty(dose_adc_sub) && sum(adc_mask_1d) >= min_subvol_voxels
                d95_adc_sub(j,k) = prctile(dose_adc_sub, 5); 
                v50_adc_sub(j,k) = sum(dose_adc_sub >= 50) / length(dose_adc_sub) * 100;
            end
            
            d_mask_1d = d_vec < d_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = d_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                d_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_d_sub = dose_vec(d_mask_1d);
            if ~isempty(dose_d_sub) && sum(d_mask_1d) >= min_subvol_voxels
                d95_d_sub(j,k) = prctile(dose_d_sub, 5);
                v50_d_sub(j,k) = sum(dose_d_sub >= 50) / length(dose_d_sub) * 100;
            end
            
            f_mask_1d = f_vec < f_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = f_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                f_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_f_sub = dose_vec(f_mask_1d);
            if ~isempty(dose_f_sub) && sum(f_mask_1d) >= min_subvol_voxels
                d95_f_sub(j,k) = prctile(dose_f_sub, 5);
                v50_f_sub(j,k) = sum(dose_f_sub >= 50) / length(dose_f_sub) * 100;
            end
            
            dstar_mask_1d = dstar_vec < dstar_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = dstar_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                dstar_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_dstar_sub = dose_vec(dstar_mask_1d);
            if ~isempty(dose_dstar_sub) && sum(dstar_mask_1d) >= min_subvol_voxels
                d95_dstar_sub(j,k) = prctile(dose_dstar_sub, 5);
                v50_dstar_sub(j,k) = sum(dose_dstar_sub >= 50) / length(dose_dstar_sub) * 100;
            end
            
        end
    end
end
end
