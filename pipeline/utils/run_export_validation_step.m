function run_export_validation_step(session, validation_model_path)
    fprintf('\xf0\x9f\x92\xa1 External validation: exporting trained model...\n');
    config_struct = session.config_struct;
    if exist(session.results_file, 'file')
        tmp = load(session.results_file, 'calculated_results');
        cr = tmp.calculated_results;
        % Build trained_model_struct from calculated_results if it has model info
        if isfield(cr, 'trained_model')
            prepare_external_validation(cr.trained_model, config_struct, validation_model_path);
        else
            fprintf('      \xe2\x9a\xa0\xef\xb8\x8f calculated_results does not contain trained_model. Export skipped.\n');
        end
    else
        fprintf('      \xe2\x9a\xa0\xef\xb8\x8f No calculated_results file found. Export skipped.\n');
    end
    fprintf('      \xe2\x9c\x85 Done.\n');
end
