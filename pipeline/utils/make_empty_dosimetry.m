function dr = make_empty_dosimetry(baseline_results)
% MAKE_EMPTY_DOSIMETRY  Create NaN-filled dosimetry struct when results are missing.
    nPt = length(baseline_results.m_id_list);
    nTp_val = baseline_results.nTp;
    dr.d95_adc_sub = nan(nPt, nTp_val);
    dr.v50_adc_sub = nan(nPt, nTp_val);
    dr.d95_d_sub = nan(nPt, nTp_val);
    dr.v50_d_sub = nan(nPt, nTp_val);
    dr.d95_f_sub = nan(nPt, nTp_val);
    dr.v50_f_sub = nan(nPt, nTp_val);
    dr.d95_dstar_sub = nan(nPt, nTp_val);
    dr.v50_dstar_sub = nan(nPt, nTp_val);
end
