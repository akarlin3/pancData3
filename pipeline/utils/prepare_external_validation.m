function prepare_external_validation(trained_model_struct, config_struct, output_path)
% PREPARE_EXTERNAL_VALIDATION  Export trained model for external validation.
%
%   Saves a self-contained .mat file with all information needed to apply
%   the trained elastic net model to an external dataset: coefficients,
%   feature names, scaling parameters, imputation reference data, and
%   the risk threshold from Youden index.
%
% Inputs:
%   trained_model_struct - Struct with fields:
%       .coefficients     - Elastic net coefficients
%       .selected_features - Indices of selected features
%       .feature_names    - Cell array of feature names
%       .scaling_mu       - Mean per feature (training set)
%       .scaling_sigma    - Std per feature (training set)
%       .imputation_ref   - Reference data for KNN imputation
%       .risk_threshold   - Optimal threshold from Youden index
%       .auc              - Training AUC
%       .n_patients       - Number of training patients
%       .event_rate       - Event rate in training cohort
%   config_struct        - Pipeline configuration struct
%   output_path          - Full path to save the .mat file

    fprintf('  --- Exporting Validation Model ---\n');

    % Validate required fields
    required_fields = {'coefficients', 'selected_features', 'feature_names', ...
        'scaling_mu', 'scaling_sigma', 'risk_threshold'};
    for i = 1:numel(required_fields)
        if ~isfield(trained_model_struct, required_fields{i})
            error('prepare_external_validation:missingField', ...
                'trained_model_struct missing required field: %s', required_fields{i});
        end
    end

    % Build export struct
    validation_model = struct();
    validation_model.model_coefficients = trained_model_struct.coefficients;
    validation_model.feature_names = trained_model_struct.feature_names;
    validation_model.selected_features = trained_model_struct.selected_features;

    validation_model.feature_scaling = struct();
    validation_model.feature_scaling.mu = trained_model_struct.scaling_mu;
    validation_model.feature_scaling.sigma = trained_model_struct.scaling_sigma;

    validation_model.risk_threshold = trained_model_struct.risk_threshold;

    % Optional fields
    if isfield(trained_model_struct, 'imputation_ref')
        validation_model.imputation_reference = trained_model_struct.imputation_ref;
    end

    % Training cohort summary
    validation_model.training_cohort_summary = struct();
    if isfield(trained_model_struct, 'n_patients')
        validation_model.training_cohort_summary.n = trained_model_struct.n_patients;
    end
    if isfield(trained_model_struct, 'event_rate')
        validation_model.training_cohort_summary.event_rate = trained_model_struct.event_rate;
    end
    if isfield(trained_model_struct, 'auc')
        validation_model.training_cohort_summary.auc = trained_model_struct.auc;
    end

    % Version and metadata
    validation_model.version = '1.0';
    validation_model.export_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    if isfield(config_struct, 'dwi_type')
        validation_model.dwi_type = config_struct.dwi_type;
    end

    % Save
    save(output_path, 'validation_model', '-v7.3');
    fprintf('  📁 Validation model exported: %s\n', output_path);
    fprintf('    Features: %d, Threshold: %.3f\n', ...
        numel(validation_model.feature_names), validation_model.risk_threshold);
end
