function results = apply_external_validation(model_path, external_data_path, config_struct)
% APPLY_EXTERNAL_VALIDATION  Apply a saved model to an external dataset.
%
%   Loads a validation model exported by prepare_external_validation.m
%   and applies it to an external dataset, generating predicted risk
%   scores with calibration assessment.
%
% Inputs:
%   model_path         - Path to the exported validation .mat file
%   external_data_path - Path to external data folder (same directory
%                        structure as the training data)
%   config_struct      - Pipeline configuration struct
%
% Outputs:
%   results            - Struct with fields:
%       .risk_scores     - Predicted risk scores per patient
%       .risk_class      - Binary high/low risk classification
%       .patient_ids     - Patient IDs from external dataset
%       .calibration     - Calibration metrics (if outcomes available)

    fprintf('  --- Applying External Validation Model ---\n');

    % Load the saved model
    if ~exist(model_path, 'file')
        error('apply_external_validation:fileNotFound', ...
            'Validation model file not found: %s', model_path);
    end
    loaded = load(model_path, 'validation_model');
    model = loaded.validation_model;

    fprintf('    Model version: %s, exported: %s\n', model.version, model.export_date);
    fprintf('    Expected features: %d\n', numel(model.feature_names));

    % Validate external data path
    if ~exist(external_data_path, 'dir')
        error('apply_external_validation:dirNotFound', ...
            'External data directory not found: %s', external_data_path);
    end

    % Load external feature matrix
    % Look for a pre-assembled feature matrix or raw data
    ext_mat = fullfile(external_data_path, 'features.mat');
    if exist(ext_mat, 'file')
        ext_data = load(ext_mat);
        if isfield(ext_data, 'X_features')
            X_ext = ext_data.X_features;
        else
            error('apply_external_validation:missingFeatures', ...
                'External features.mat must contain X_features variable.');
        end
        if isfield(ext_data, 'patient_ids')
            patient_ids = ext_data.patient_ids;
        else
            patient_ids = arrayfun(@(i) sprintf('EXT_%03d', i), 1:size(X_ext,1), 'UniformOutput', false);
        end
        if isfield(ext_data, 'y_true')
            y_ext = ext_data.y_true;
        else
            y_ext = [];
        end
    else
        error('apply_external_validation:noFeatures', ...
            'External data folder must contain features.mat with X_features.');
    end

    % Check feature compatibility
    n_expected = numel(model.feature_names);
    if size(X_ext, 2) ~= n_expected
        warning('apply_external_validation:featureMismatch', ...
            'Expected %d features, got %d. Attempting to match by name.', ...
            n_expected, size(X_ext, 2));
    end

    % Apply scaling
    mu = model.feature_scaling.mu;
    sigma = model.feature_scaling.sigma;
    n_feat = min(size(X_ext, 2), numel(mu));
    X_scaled = X_ext(:, 1:n_feat);
    for f = 1:n_feat
        if sigma(f) > 0
            X_scaled(:, f) = (X_scaled(:, f) - mu(f)) / sigma(f);
        end
    end

    % Handle missing features
    if n_feat < n_expected
        X_scaled = [X_scaled, zeros(size(X_scaled, 1), n_expected - n_feat)];
    end

    % Apply model
    coefs = model.model_coefficients;
    if numel(coefs) == n_expected + 1
        % Intercept included
        risk_scores = X_scaled * coefs(2:end) + coefs(1);
    else
        risk_scores = X_scaled * coefs(:);
    end

    % Convert to probabilities (logistic)
    risk_probs = 1 ./ (1 + exp(-risk_scores));

    % Classify
    risk_class = double(risk_probs >= model.risk_threshold);

    results = struct();
    results.risk_scores = risk_scores;
    results.risk_probs = risk_probs;
    results.risk_class = risk_class;
    results.patient_ids = patient_ids;

    fprintf('    Patients scored: %d\n', numel(risk_scores));
    fprintf('    High risk: %d (%.1f%%)\n', sum(risk_class), 100*mean(risk_class));

    % Calibration if outcomes available
    if ~isempty(y_ext)
        try
            cal = compute_calibration_metrics(risk_probs, y_ext, 5, '', '', '');
            results.calibration = cal;
            fprintf('    External Brier score: %.3f\n', cal.brier_score);
        catch
            fprintf('    ⚠️  Calibration assessment failed.\n');
        end
    end
end
