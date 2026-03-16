function baseline = load_baseline_from_disk(baseline_results_file)
% LOAD_BASELINE_FROM_DISK  Load persisted metrics_baseline outputs.
%
%   baseline = load_baseline_from_disk(baseline_results_file)
%
%   Loads the .mat file saved by the metrics_baseline step and returns
%   all variables in a single struct.  This is used when the pipeline
%   skips metrics_baseline but downstream steps (longitudinal, dosimetry,
%   stats, survival) need the baseline outputs.
%
%   Inputs:
%     baseline_results_file - Full path to metrics_baseline_results_*.mat
%
%   Outputs:
%     baseline - Struct with all baseline output fields:
%                .m_lf, .m_total_time, .m_total_follow_up_time,
%                .m_gtv_vol, .m_adc_mean, .m_d_mean, .m_f_mean,
%                .m_dstar_mean, .m_id_list, .m_mrn_list, .m_d95_gtvp,
%                .m_v50gy_gtvp, .m_data_vectors_gtvp, .lf_group,
%                .valid_pts, .ADC_abs, .D_abs, .f_abs, .Dstar_abs,
%                .ADC_pct, .D_pct, .f_delta, .Dstar_pct, .nTp,
%                .metric_sets, .set_names, .time_labels, .dtype_label,
%                .dl_provenance
%
%   Errors:
%     Throws 'LoadBaseline:NotFound' if the file does not exist.

    if ~exist(baseline_results_file, 'file')
        error('LoadBaseline:NotFound', ...
            'metrics_baseline results not found at: %s', baseline_results_file);
    end

    baseline = load(baseline_results_file);
end
