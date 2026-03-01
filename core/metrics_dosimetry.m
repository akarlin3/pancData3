function [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub] = metrics_dosimetry(m_id_list, id_list, nTp, config_struct, m_data_vectors_gtvp, gtv_locations)
% METRICS_DOSIMETRY — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 3/5 of the metrics step. Computes dose metrics (D95, V50) for resistant sub-volumes.
%
% Inputs:
%   m_id_list         - Cell array of valid patient identifiers
%   id_list           - Original cell array of all patient identifiers
%   nTp               - Number of timepoints
%   config_struct     - Configuration settings (thresholds)
%   m_data_vectors_gtvp- Struct array containing parametric vectors and dose vectors
%   gtv_locations     - Cell array showing paths to 3D mask arrays
%
% Outputs:
%   d95_*_sub         - D95 dose (Gy) delivered to the resistant sub-volume defined by *
%   v50_*_sub         - V50Gy dose coverage (%) delivered to the resistant sub-volume defined by *
%

fprintf('  --- SECTION 6: Target Coverage (Sub-Volume Dose Metrics) ---\n');

adc_thresh   = config_struct.adc_thresh;    
d_thresh     = config_struct.d_thresh;      
f_thresh     = config_struct.f_thresh;      
dstar_thresh = config_struct.dstar_thresh;  

if iscell(nTp)
    nTp = nTp{1};
end

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
            gtv_mask_3d = [];
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
            
            [d95_adc_sub(j,k), v50_adc_sub(j,k)] = calculate_subvolume_metrics(adc_vec, adc_thresh, dose_vec, has_3d, gtv_mask_3d);
            [d95_d_sub(j,k), v50_d_sub(j,k)]     = calculate_subvolume_metrics(d_vec, d_thresh, dose_vec, has_3d, gtv_mask_3d);
            [d95_f_sub(j,k), v50_f_sub(j,k)]     = calculate_subvolume_metrics(f_vec, f_thresh, dose_vec, has_3d, gtv_mask_3d);
            [d95_dstar_sub(j,k), v50_dstar_sub(j,k)] = calculate_subvolume_metrics(dstar_vec, dstar_thresh, dose_vec, has_3d, gtv_mask_3d);
            
        end
    end
end
end
