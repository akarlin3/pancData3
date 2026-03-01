function [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts)
% FIT_MODELS Fits IVIM and ADC diffusion models to DWI data
%
%   [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts)
%   Flattens the 3D masked volume to 1D, calls IVIMmodelfit, and calculates
%   monoexponential ADC using vectorized OLS. Reconstructs outputs to 3D.

    % [MODULARIZATION STAGE 3]: Masked 1D Flattening + `parfor`
    % Flatten 3D volume to a 1D array of strictly non-zero mask voxels
    % to bypass thousands of empty space computations within the dependency.
    sz3 = [size(dwi,1), size(dwi,2), size(dwi,3)];
    valid_voxels_idx = find(mask_ivim);
    n_valid = length(valid_voxels_idx);

    % Preallocate output 1D arrays
    d_vec = nan(n_valid, 1);
    f_vec = nan(n_valid, 1);
    dstar_vec = nan(n_valid, 1);

    has_sufficient_bvalues = sum(bvalues >= opts.bthr) >= 2;
    if has_sufficient_bvalues && n_valid > 0
        fprintf('  [Stage 3 Opt] Flattening %d valid voxels for accelerated IVIM fit...\n', n_valid);

        % Extract 1D signal decay curves for valid voxels
        % Reshape DWI to (voxels x bvalues)
        dwi_flat = reshape(dwi, [prod(sz3), length(bvalues)]);
        dwi_valid = dwi_flat(valid_voxels_idx, :);

        % [MODULARIZATION STAGE 3]: Masked 1D Flattening
        % Pad to even number of elements to ensure 3rd dimension > 1
        % so MATLAB doesn't drop it in size() checks inside dependency.
        pad_len = mod(2 - mod(n_valid, 2), 2);
        dwi_valid_padded = [dwi_valid; zeros(pad_len, length(bvalues))];
        n_padded = n_valid + pad_len;

        % We must pass a 3D volume to the dependency, [N, 1, 2, bval]
        dwi_1d_vol = reshape(dwi_valid_padded, [n_padded/2, 1, 2, length(bvalues)]);
        mask_1d_vol = true(n_padded/2, 1, 2);

        % Execute the untouched dependency on the flattened array
        ivim_fit_1d = IVIMmodelfit(dwi_1d_vol, bvalues, "seg", mask_1d_vol, opts);

        % Restructure output back to strictly 1D and snip padding
        ivim_out_flat = reshape(ivim_fit_1d, [n_padded, 4]);
        d_vec = squeeze(ivim_out_flat(1:n_valid, 1));
        f_vec = squeeze(ivim_out_flat(1:n_valid, 3));
        dstar_vec = squeeze(ivim_out_flat(1:n_valid, 4));

        % Replace zero-fit voxels with NaN (failed fits)
        zero_mask = (d_vec == 0);
        d_vec(zero_mask) = nan;
        f_vec(zero_mask) = nan;
        dstar_vec(zero_mask) = nan;
    elseif ~has_sufficient_bvalues
        fprintf('Insufficient b-values >= %d for IVIM fit; skipping IVIM (maps set to NaN)\n', opts.bthr);
    end

    % Safely reconstruct 1D outputs back into 3D volume geometry
    d_map = nan(sz3);
    f_map = nan(sz3);
    dstar_map = nan(sz3);

    if n_valid > 0
        d_map(valid_voxels_idx) = d_vec;
        f_map(valid_voxels_idx) = f_vec;
        dstar_map(valid_voxels_idx) = dstar_vec;
    end

    % Monoexponential ADC fit — log-linear OLS on active voxels only.
    % [PHYSICS EXPLANATION]:
    % Apparent Diffusion Coefficient (ADC) assumes a simple mono-exponential decay:
    % S = S0 * exp(-b * ADC)
    % Linearized form: ln(S/S0) = -b * ADC
    % Solved via Ordinary Least Squares (OLS) on the log-signal.

    adc_sz  = [size(dwi,1), size(dwi,2), size(dwi,3)];
    adc_map = nan(adc_sz);

    if n_valid > 0
        % Extract 1D signal decay curves if not already computed
        if ~exist('dwi_valid', 'var')
            dwi_flat = reshape(dwi, [prod(sz3), length(bvalues)]);
            dwi_valid = dwi_flat(valid_voxels_idx, :);
        end

        % Filter to voxels with all-positive signal to prevent log() issues
        adc_valid_idx = all(dwi_valid > 0, 2);

        if any(adc_valid_idx)
            S_a = dwi_valid(adc_valid_idx, :);

            % Vectorized OLS on valid masked voxels
            adc_vals = (-bvalues(2:end) \ ...
                permute(log(S_a(:,2:end) ./ S_a(:,1)), [2 1]))';

            adc_vals(adc_vals < 0) = nan;  % clamp noise-driven negative ADC estimates

            % Prepare a temporary 1D vector for all valid mask voxels
            adc_vec_out = nan(n_valid, 1);
            adc_vec_out(adc_valid_idx) = adc_vals;

            % Reconstruct into 3D geometry
            adc_map(valid_voxels_idx) = adc_vec_out;
        end
    end
end
