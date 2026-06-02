function [data_vectors_gtvp, data_vectors_gtvn, dmean_gtvp, dmean_gtvn, ...
    d95_gtvp, d95_gtvn, v50gy_gtvp, v50gy_gtvn, adc_mean, adc_kurtosis, ...
    d_mean, d_kurtosis, d_mean_dncnn, d_mean_ivimnet, lf, immuno, ...
    bad_dwi_locations, bad_dwi_count] = load_dwi_data_reconstruct( ...
    mrn_list, id_list, checkpoint_dir, dwi_locations, ...
    data_vectors_gtvp, data_vectors_gtvn, dmean_gtvp, dmean_gtvn, ...
    d95_gtvp, d95_gtvn, v50gy_gtvp, v50gy_gtvn, adc_mean, adc_kurtosis, ...
    d_mean, d_kurtosis, d_mean_dncnn, d_mean_ivimnet, lf, immuno, ...
    bad_dwi_locations_per_patient)
    % LOAD_DWI_DATA_RECONSTRUCT Reconstruct cohort-wide arrays from per-patient checkpoints.
    %   Extracted verbatim from load_dwi_data.m Section 2 (post-parfor reconstruction
    %   loop). Loads each patient's checkpoint .mat, validates it, and assigns the
    %   results back into the pre-initialized cohort-level arrays passed in. Returns
    %   the populated arrays plus the flattened bad-DWI location list and count.

% Reconstruct global arrays from per-patient checkpoints.
% Each checkpoint contains one patient's complete processing results
% (all fractions, all repeats). This reconstruction step assembles the
% individual patient results into cohort-wide arrays that downstream
% analysis expects. The separation between per-patient checkpointing
% (during parfor) and global reconstruction (sequential) is necessary
% because parfor does not support direct assignment to shared struct arrays
% with dynamic field sets.
n_reconstruct = length(mrn_list);
for j = 1:n_reconstruct
    text_progress_bar(j, n_reconstruct, 'Reconstructing checkpoints');
    mrn = mrn_list{j};
    patient_id = id_list{j};
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, patient_id));

    if exist(checkpoint_file, 'file')
        % Load checkpoint with basic corruption detection
        loaded_data = load(checkpoint_file);

        required_fields = {'data_vectors_gtvp', 'data_vectors_gtvn', ...
            'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', ...
            'v50gy_gtvp', 'v50gy_gtvn', 'adc_mean', 'd_mean', ...
            'd_mean_dncnn', 'd_mean_ivimnet', 'lf', 'immuno', 'bad_dwi_list'};
        missing_fields = setdiff(required_fields, fieldnames(loaded_data));
        if ~isempty(missing_fields)
            warning('load_dwi_data:corruptCheckpoint', ...
                'Checkpoint for patient %d (%s) is missing fields: %s. Skipping.', ...
                j, patient_id, strjoin(missing_fields, ', '));
            continue;
        end

        % Validate dimensions of checkpoint arrays against expected sizes
        expected_n_fx = size(dwi_locations, 2);
        expected_n_rp = size(dwi_locations, 3);
        cp_adc_size = size(loaded_data.adc_mean);
        if length(cp_adc_size) < 2 || cp_adc_size(2) ~= expected_n_fx
            warning('load_dwi_data:checkpointSizeMismatch', ...
                'Checkpoint for patient %d (%s): adc_mean has %d fractions, expected %d. Skipping.', ...
                j, patient_id, cp_adc_size(2), expected_n_fx);
            continue;
        end
        if length(cp_adc_size) >= 3 && cp_adc_size(3) ~= expected_n_rp
            warning('load_dwi_data:checkpointSizeMismatch', ...
                'Checkpoint for patient %d (%s): adc_mean has %d repeats, expected %d. Skipping.', ...
                j, patient_id, cp_adc_size(3), expected_n_rp);
            continue;
        end

        % Assign back to global arrays
        % Struct arrays
        data_vectors_gtvp = align_and_assign_struct(data_vectors_gtvp, loaded_data.data_vectors_gtvp, j);
        data_vectors_gtvn = align_and_assign_struct(data_vectors_gtvn, loaded_data.data_vectors_gtvn, j);

        % Scalar/Vector arrays (patient x fraction)
        dmean_gtvp(j,:) = loaded_data.dmean_gtvp;
        dmean_gtvn(j,:) = loaded_data.dmean_gtvn;
        d95_gtvp(j,:) = loaded_data.d95_gtvp;
        d95_gtvn(j,:) = loaded_data.d95_gtvn;
        v50gy_gtvp(j,:) = loaded_data.v50gy_gtvp;
        v50gy_gtvn(j,:) = loaded_data.v50gy_gtvn;

        % Summary metrics (patient x fraction x repeat)
        adc_mean(j,:,:) = loaded_data.adc_mean;
        if isfield(loaded_data, 'adc_kurtosis')
            adc_kurtosis(j,:,:) = loaded_data.adc_kurtosis;
        end
        d_mean(j,:,:) = loaded_data.d_mean;
        if isfield(loaded_data, 'd_kurtosis')
            d_kurtosis(j,:,:) = loaded_data.d_kurtosis;
        end
        d_mean_dncnn(j,:,:) = loaded_data.d_mean_dncnn;
        d_mean_ivimnet(j,:,:) = loaded_data.d_mean_ivimnet;

        % Clinical data and tracking
        lf(j) = loaded_data.lf;
        immuno(j) = loaded_data.immuno;
        bad_dwi_locations_per_patient{j} = loaded_data.bad_dwi_list;
    else
        fprintf('Warning: No checkpoint found for patient %d (Patient ID %s) during reconstruction.\n', j, patient_id);
    end
end

% Flatten bad_dwi_locations from per-patient cell arrays into a single
% cohort-wide list.  These flagged acquisitions are reported in the
% pipeline log for the physicist to review and decide on exclusion.
bad_dwi_locations = [bad_dwi_locations_per_patient{:}];
bad_dwi_count = length(bad_dwi_locations);
end
