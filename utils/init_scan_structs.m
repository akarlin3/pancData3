function [gtvp_structs, gtvn_structs] = init_scan_structs(n_fx, n_rp)
% INIT_SCAN_STRUCTS Create empty GTVp/GTVn struct arrays with consistent field ordering.
%   [gtvp_structs, gtvn_structs] = init_scan_structs(n_fx, n_rp) returns two
%   n_fx-by-n_rp struct arrays pre-populated with all biomarker, dosimetry,
%   clinical, and deep-learning fields set to []. This ensures a uniform
%   field order across the pipeline and prevents struct concatenation errors
%   when fields are added in varying order.
%
%   Used by load_dwi_data and process_single_scan to initialise storage
%   before populating per-scan results.
%
%   Inputs:
%       n_fx  — number of fractions (rows)
%       n_rp  — number of repeatability indices (columns)
%
%   Outputs:
%       gtvp_structs — n_fx x n_rp struct array for primary GTV data
%       gtvn_structs — n_fx x n_rp struct array for nodal GTV data

    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'dose_vector', [], 'dvh', [], 'd95', [], 'v50gy', [], ...
        'd_vector_dncnn', [], 'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'adc_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], 'dstar_vector_ivimnet', [], ...
        'ID', [], 'MRN', [], 'LF', [], 'Immuno', [], ...
        'Fraction', [], 'Repeatability_index', [], 'vox_vol', []);
    gtvp_structs = repmat(empty_entry, n_fx, n_rp);
    gtvn_structs = repmat(empty_entry, n_fx, n_rp);
end
