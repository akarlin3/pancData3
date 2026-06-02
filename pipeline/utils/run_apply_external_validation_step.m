function run_apply_external_validation_step(session, ext_data_path, validation_model_path)
    fprintf('\xf0\x9f\x92\xa1 External validation: applying saved model to external dataset\n');
    config_struct = session.config_struct;
    if ~exist(validation_model_path, 'file')
        fprintf('      \xe2\x9a\xa0\xef\xb8\x8f Validation model not found at: %s. Skipping.\n', validation_model_path);
        return;
    end
    ext_results = apply_external_validation(validation_model_path, ext_data_path, config_struct);
    % Store results in the session results file
    if exist(session.results_file, 'file')
        tmp = load(session.results_file, 'calculated_results');
        cr = tmp.calculated_results;
        cr.external_validation = ext_results;
        calculated_results = cr; %#ok<NASGU>
        save(session.results_file, 'calculated_results');
        fprintf('      \xf0\x9f\x93\x81 External validation results saved to %s\n', session.results_file);
    end
    fprintf('      \xe2\x9c\x85 Done.\n');
end
