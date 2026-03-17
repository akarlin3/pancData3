function [X_lasso_all, feat_names_lasso, original_feature_indices, feat_names_lasso_full] = assemble_predictive_features( ...
    valid_pts, target_fx, nTp, fx_label, output_folder, ...
    ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
    m_d95_gtvp, m_v50gy_gtvp, ...
    d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
    d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, ...
    adc_kurt, adc_skew, auxiliary_features, auxiliary_feature_names)
% ASSEMBLE_PREDICTIVE_FEATURES  Build the feature matrix for elastic net.
%
%   Assembles candidate features from baseline covariates (Fx1 absolute
%   values), absolute values at the target fraction, change metrics, and
%   dose features.  Post-treatment timepoints exclude dose columns.
%   When kurtosis and skewness arrays are provided, they are appended as
%   columns 23-24 (yielding a 24-column matrix before NaN removal).
%   Optional auxiliary features (e.g., biomarkers) are appended after
%   the radiomics columns when provided.
%   All-NaN columns are removed.
%
% Inputs:
%   valid_pts              - Logical mask of valid patients
%   target_fx              - Target fraction index (column in parameter arrays)
%   nTp                    - Total number of timepoints
%   fx_label               - Label for the current fraction (for error messages)
%   output_folder          - Output folder path (for debug file on error)
%   ADC_abs ... v50_dstar_sub - Parameter arrays (n_patients x nTp)
%   adc_kurt               - (Optional, arg 23) ADC kurtosis (n_patients x nTp)
%   adc_skew               - (Optional, arg 24) ADC skewness (n_patients x nTp)
%   auxiliary_features     - (Optional, arg 25) [n_patients x n_aux] auxiliary features
%   auxiliary_feature_names - (Optional, arg 26) {1 x n_aux} cell array of names
%
% Outputs:
%   X_lasso_all             - Feature matrix with all-NaN columns removed
%   feat_names_lasso        - Feature names after column removal
%   original_feature_indices - 1-based indices into the original features
%   feat_names_lasso_full   - Unfiltered feature names (before NaN column removal)

    % Assemble the feature matrix with 22 candidate predictors.
    % Include baseline (Fx1) absolute values as covariates to adjust for
    % pre-existing differences in tumor burden and diffusion properties.
    % Without baseline adjustment, observed differences in change metrics
    % could be confounded by baseline tumour characteristics (e.g., larger
    % tumours may show smaller percent changes due to dilution effects).
    X_lasso_all = [ADC_abs(valid_pts, 1), D_abs(valid_pts, 1), ...
                   f_abs(valid_pts, 1),   Dstar_abs(valid_pts, 1), ...
                   ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
                   f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
                   ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
                   f_delta(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx), ...
                   m_d95_gtvp(valid_pts, target_fx), m_v50gy_gtvp(valid_pts, target_fx), ...
                   d95_adc_sub(valid_pts, target_fx), v50_adc_sub(valid_pts, target_fx), ...
                   d95_d_sub(valid_pts, target_fx),   v50_d_sub(valid_pts, target_fx), ...
                   d95_f_sub(valid_pts, target_fx),   v50_f_sub(valid_pts, target_fx), ...
                   d95_dstar_sub(valid_pts, target_fx), v50_dstar_sub(valid_pts, target_fx)];


    % Feature name labels corresponding to the columns above:
    %   Cols  1-4:  Baseline (Fx1) absolute diffusion values (_BL)
    %   Cols  5-8:  Absolute diffusion values at the target fraction (_Abs)
    %   Cols  9-12: Percent/delta change from baseline (_Pct, _Delta)
    %   Cols 13-14: Dose metrics for the whole GTVp (D95 = dose covering 95%
    %               of volume; V50 = fraction receiving >= 50 Gy)
    %   Cols 15-22: Dose metrics for parameter-specific sub-volumes
    %               (e.g., D95_Sub_ADC = D95 within the ADC-defined core)
    %   Cols 23-24: ADC kurtosis and skewness (optional, from summary_metrics)
    feat_names_lasso = {'ADC_BL', 'D_BL', 'f_BL', 'Dstar_BL', ...
                        'ADC_Abs', 'D_Abs', 'f_Abs', 'Dstar_Abs', ...
                        'ADC_Pct', 'D_Pct', 'f_Delta', 'Dstar_Pct', ...
                        'D95_GTVp', 'V50_GTVp', ...
                        'D95_Sub_ADC', 'V50_Sub_ADC', ...
                        'D95_Sub_D', 'V50_Sub_D', ...
                        'D95_Sub_f', 'V50_Sub_f', ...
                        'D95_Sub_Dstar', 'V50_Sub_Dstar'};

    n_base_features = 22;

    % Append ADC kurtosis and skewness if provided (columns 23-24)
    if nargin >= 24 && ~isempty(adc_kurt) && ~isempty(adc_skew)
        X_lasso_all = [X_lasso_all, ...
                       adc_kurt(valid_pts, target_fx), ...
                       adc_skew(valid_pts, target_fx)];
        feat_names_lasso = [feat_names_lasso, 'ADC_Kurt', 'ADC_Skew'];
        n_base_features = 24;
    end

    % Append auxiliary features if provided (e.g., ctDNA, PET SUV, genomic)
    if nargin >= 26 && ~isempty(auxiliary_features) && ~isempty(auxiliary_feature_names)
        X_lasso_all = [X_lasso_all, auxiliary_features(valid_pts, :)];
        feat_names_lasso = [feat_names_lasso, auxiliary_feature_names(:)'];
        n_base_features = n_base_features + size(auxiliary_features, 2);
    end

    % Track original column positions through filtering so we can map
    % surviving features back to their canonical names after NaN removal
    original_feature_indices = 1:n_base_features;

    if target_fx == nTp || target_fx == 6
        % Post-treatment timepoint: exclude dose features (columns 13-22).
        % Dose metrics are only meaningful during active treatment when the
        % dose is being delivered.  At the post-treatment scan (typically
        % 3 months after RT), the full dose has been delivered and there is
        % no additional dose to correlate with — the therapeutic window for
        % dose-response analysis has closed.
        % Keep columns 1-12, skip 13-22, keep 23-24 (kurtosis/skewness) if present
        if n_base_features > 22
            keep_cols = [1:12, 23:n_base_features];
        else
            keep_cols = 1:12;
        end
        X_lasso_all = X_lasso_all(:, keep_cols);
        feat_names_lasso = feat_names_lasso(keep_cols);
        original_feature_indices = original_feature_indices(keep_cols);
    end

    % Defensive check: if any input array has incompatible types (e.g.,
    % cell instead of numeric), MATLAB's horzcat produces a cell array
    % instead of a numeric matrix. Catch this early and dump diagnostic
    % info for debugging rather than failing cryptically downstream.
    if iscell(X_lasso_all)
        vars = {'ADC_abs_BL', 'D_abs_BL', 'f_abs_BL', 'Dstar_abs_BL', 'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', 'm_d95_gtvp', 'm_v50gy_gtvp', 'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub'};
        vars_vals = {ADC_abs, D_abs, f_abs, Dstar_abs, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub};
        fp = fopen(fullfile(output_folder, 'debug_concat_error.txt'), 'a');
        if fp >= 0
            fprintf(fp, '\n--- X_lasso_all is a cell array at fx_label=%s ---\n', fx_label);
            for i_v = 1:length(vars)
                tmp_v = vars_vals{i_v};
                fprintf(fp, '%s -> Size: %s, Class: %s\n', vars{i_v}, mat2str(size(tmp_v)), class(tmp_v));
            end
            fclose(fp);
        end
        error('Invalid data type: cell array detected in X_lasso_all concatenation. See debug_concat_error.txt');
    end

    % Remove all-NaN columns (features with no valid data for any patient).
    % This typically occurs when dose data is unavailable or when a DWI type
    % lacks IVIM parameters (e.g., D/f/D* may be all-NaN for Standard DWI
    % if only ADC was fitted).
    valid_cols = ~all(isnan(X_lasso_all), 1);
    feat_names_lasso_full = feat_names_lasso;  % keep unfiltered copy for display
    X_lasso_all = X_lasso_all(:, valid_cols);
    feat_names_lasso = feat_names_lasso(valid_cols);
    original_feature_indices = original_feature_indices(valid_cols);
end
