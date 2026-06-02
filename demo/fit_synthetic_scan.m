function rec = fit_synthetic_scan(dwi, bvalues, mask, bthr)
% FIT_SYNTHETIC_SCAN  Run the REAL pipeline fitter on one synthetic scan.
%
% This is the whole point of the demo: nothing here re-implements the
% physics fit. It calls pipeline/core/fit_models.m unmodified — the exact
% segmented IVIM + weighted-ADC estimator the clinical pipeline uses — on
% phantom data whose ground truth we know. The verbose progress prints from
% fit_models are captured so the demo log stays readable.
%
% Returns per-voxel recovered parameter vectors (over the mask), in the same
% column-major voxel order as truth maps indexed by the same logical mask.
    opts = struct('bthr', bthr);
    % evalc captures fit_models' fprintf chatter; the assignment still lands
    % in this workspace (evalc runs in the caller scope).
    evalc('[d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask, opts);'); %#ok<NASGU>
    m = (mask == 1);
    rec = struct('D', d_map(m), 'f', f_map(m), 'Dstar', dstar_map(m), 'ADC', adc_map(m));
end
