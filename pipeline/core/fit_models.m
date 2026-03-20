function [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts)
% FIT_MODELS Fits IVIM and ADC diffusion models to DWI data
%
%   [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts)
%   Flattens the 3D masked volume to 1D, calls IVIMmodelfit, and calculates
%   monoexponential ADC using vectorized OLS. Reconstructs outputs to 3D.
%
%   opts.use_gpu (optional, default false): When true and a GPU is available,
%   the ADC weighted least-squares computation is offloaded to the GPU via
%   gpuArray for faster vectorized matrix operations.
%
% ANALYTICAL RATIONALE — TWO-MODEL APPROACH
%   This function fits two complementary diffusion models to multi-b-value
%   DWI data within the GTV mask:
%
%   1. IVIM biexponential model (segmented fit):
%      S(b) = S0 * [f * exp(-b*D*) + (1-f) * exp(-b*D)]
%      Separates true tissue diffusion (D) from pseudo-diffusion (D*) and
%      perfusion fraction (f). In pancreatic tumors, D reflects cellular
%      density (restricted diffusion = high cellularity = low D), while f
%      and D* reflect microvasculature and capillary blood flow. Changes in
%      these parameters during RT can indicate tumor response before
%      morphological changes are visible.
%
%   2. ADC monoexponential model (weighted least-squares):
%      S(b) = S0 * exp(-b*ADC)
%      A simpler model that combines diffusion and perfusion into a single
%      "apparent" coefficient. While less physically specific than IVIM, ADC
%      is more robust and widely used clinically. ADC < ~1.0e-3 mm^2/s in
%      pancreatic tissue suggests restricted diffusion (high cellularity).
%
%   Both models are fit only within the GTV mask to avoid wasting compute
%   on background voxels and to prevent normal tissue from contaminating
%   tumor-specific parameter distributions.

    %% ---- Input validation ----
    % Validate that bvalues(1)==0 up-front before any model fitting.
    % The ADC WLS formula divides by S(:,1) assuming it is S(b=0), and the
    % IVIM segmented fit also relies on S0 being the b=0 signal. A non-zero
    % first b-value produces silently wrong estimates for both models.
    if bvalues(1) ~= 0
        error('fit_models:noB0', ...
            'First b-value is %g, not 0. Both IVIM and ADC fits require S(b=0) as reference. Reorder b-values so b=0 comes first.', bvalues(1));
    end

    % Handle repeated b-values by averaging signal across duplicates.
    % Clinical DWI protocols legitimately include repeated b-values for
    % signal averaging (e.g., b=[0,0,100,100,500,500]) where the scanner
    % acquires multiple averages at each b-value. We detect duplicates,
    % average the corresponding DWI volumes, and proceed with unique b-values.
    [bvalues_unique, ~, ic] = unique(bvalues);
    if length(bvalues_unique) < length(bvalues)
        warning('fit_models:repeatedBValues', ...
            'Detected %d repeated b-values (e.g., multiple averages). Averaging signal across %d unique b-values before fitting.', ...
            length(bvalues) - length(bvalues_unique), length(bvalues_unique));
        
        % Average DWI signal across repeated b-values
        sz_dwi = size(dwi);
        dwi_flat_all = reshape(dwi, [], length(bvalues)); % [n_spatial x n_bvalues]
        n_unique = length(bvalues_unique);
        dwi_averaged = zeros(size(dwi_flat_all, 1), n_unique);
        for ub = 1:n_unique
            idx_this_b = (ic == ub);
            dwi_averaged(:, ub) = mean(dwi_flat_all(:, idx_this_b), 2);
        end
        dwi = reshape(dwi_averaged, [sz_dwi(1), sz_dwi(2), sz_dwi(3), n_unique]);
        bvalues = bvalues_unique;
        clear dwi_flat_all dwi_averaged;
    end

    % Validate that (unique) b-values are monotonically increasing
    if any(diff(bvalues) <= 0)
        error('fit_models:nonMonotonicBValues', ...
            'B-values must be monotonically increasing after averaging duplicates. Non-monotonic ordering can cause convergence failures in IVIM model fitting.');
    end

    % Validate sufficient dynamic range for IVIM parameter estimation
    % Need adequate separation between low-b and high-b values to distinguish
    % perfusion and diffusion components. Minimum b-value range of 200 s/mm²
    % provides sufficient signal decay for reliable parameter estimation.
    bvalue_range = bvalues(end) - bvalues(1);
    min_bvalue_range = 200; % s/mm²
    if bvalue_range < min_bvalue_range
        warning('fit_models:insufficientBValueRange', ...
            'B-value range (%.1f s/mm²) is less than recommended minimum (%.1f s/mm²). IVIM parameter estimation may be unreliable.', ...
            bvalue_range, min_bvalue_range);
    end

    % Validate minimum number of b-values for model fitting
    if length(bvalues) < 4
        warning('fit_models:fewBValues', ...
            'Only %d b-values provided. IVIM fitting requires at least 4 b-values for reliable parameter estimation.', ...
            length(bvalues));
    end

    %% ---- Optimized solver configuration for IVIM fitting ----
    % Configure optimized solver options for lsqnonlin used in IVIM fitting
    % Improved settings for better convergence speed and robustness
    optimized_opts = optimoptions('lsqnonlin', ...
        'UseParallel', true, ...                 % Enable parallel gradient computation
        'OptimalityTolerance', 1e-8, ...        % Tighter optimality tolerance for better convergence
        'FunctionTolerance', 1e-8, ...          % Tighter function tolerance
        'StepTolerance', 1e-10, ...             % Smaller step tolerance for precision
        'MaxIterations', 400, ...               % Increased max iterations
        'MaxFunctionEvaluations', 800, ...      % Increased max function evaluations
        'Display', 'off');                      % Suppress individual fit output
    
    % Store original options to restore later if needed
    if isfield(opts, 'optimoptions')
        original_opts = opts.optimoptions;
    else
        original_opts = [];
    end
    
    % Apply optimized options
    opts.optimoptions = optimized_opts;

    %% ---- IVIM MODEL FITTING (Segmented Biexponential) ----
    % [MODULARIZATION STAGE 3]: Masked 1D Flattening + `parfor`
    % Flatten 3D volume to a 1D array of strictly non-zero mask voxels
    % to bypass thousands of empty space computations within the dependency.
    % The GTV mask typically covers only ~100-2000 voxels out of a
    % ~256x256x20 volume (~1.3M voxels), so masking provides a ~1000x speedup.
    sz3 = [size(dwi,1), size(dwi,2), size(dwi,3)];  % spatial dimensions of DWI volume
    valid_voxels_idx = find(mask_ivim);  % linear indices of tumor voxels within the GTV mask
    n_valid = length(valid_voxels_idx);  % number of voxels to fit (typically ~100-2000)

    % Preallocate output 1D arrays with NaN so that voxels where fitting
    % fails (convergence failure, non-physical result) naturally propagate
    % as missing data through nanmean/nanstd in downstream summary metrics.
    d_vec = nan(n_valid, 1);
    f_vec = nan(n_valid, 1);
    dstar_vec = nan(n_valid, 1);

    % The segmented IVIM fit requires:
    %   - At least 2 b-values >= bthr (typically 100 s/mm^2) for the high-b
    %     linear fit that estimates D (Stage 1 of segmented approach)
    %   - At least 1 b-value < bthr to separate the perfusion component
    %     (Stage 2 estimates f and D*)
    % Without these, the model is underdetermined. Typical pancreatic DWI
    % protocols use b = [0, 30, 100, 150, 550] which gives 3 high-b and
    % 2 low-b values — sufficient for segmented fitting.
    has_sufficient_bvalues = sum(bvalues >= opts.bthr) >= 2 && sum(bvalues < opts.bthr) >= 1;
    
    if ~has_sufficient_bvalues
        fprintf('Insufficient b-values >= %d for IVIM fit; skipping IVIM (maps set to NaN)\n', opts.bthr);
        % Early return for IVIM parameters - ADC can still be computed
        d_map = nan(sz3);
        f_map = nan(sz3);
        dstar_map = nan(sz3);
    else
        % Pre-compute ADC-derived initial guess for D parameter to improve convergence
        adc_initial_guess = [];
        if n_valid > 0
            fprintf('  [Stage 3 Opt] Computing ADC warm start for IVIM D parameter...\n');
            
            % Extract 1D signal decay curves for valid voxels for initial ADC estimation
            % Avoid creating a full reshape of the 4D volume; instead index
            % directly into each 3D sub-volume to save memory.
            n_bvals = length(bvalues);
            dwi_valid = zeros(n_valid, n_bvals);
            for b = 1:n_bvals
                vol_b = dwi(:,:,:,b);
                dwi_valid(:,b) = vol_b(valid_voxels_idx);
            end
            
            % Quick ADC estimation for warm start using high-b values only
            % This provides better initial guess for D parameter than default values
            high_b_idx = bvalues >= opts.bthr;
            if sum(high_b_idx) >= 2
                % Use only high b-values for more accurate D initial guess
                b_high = bvalues(high_b_idx);
                s_high = dwi_valid(:, high_b_idx);
                s0_valid = dwi_valid(:, 1); % b=0 signal
                
                % Simple linear regression on log data for initial D estimate
                valid_signal_idx = all(s_high > 0, 2) & s0_valid > 0;
                if any(valid_signal_idx)
                    s_ratio = s_high(valid_signal_idx, :) ./ s0_valid(valid_signal_idx);
                    log_s = log(s_ratio);
                    
                    % Linear regression: log(S/S0) = -b*D
                    A_matrix = [-b_high(:), ones(length(b_high), 1)];
                    ATA = A_matrix' * A_matrix;
                    if rcond(ATA) < eps
                        % Singular or near-singular (e.g., identical b-values); skip warm start
                        adc_initial_guess = [];
                    else
                        adc_coeffs = ATA \ (A_matrix' * log_s');
                        adc_initial_guess = max(adc_coeffs(1, :), 1e-5); % Ensure positive values
                    end
                    adc_initial_guess = min(adc_initial_guess, 3e-3); % Cap at reasonable upper limit
                end
            end
            
            if ~isempty(adc_initial_guess)
                fprintf('  [Warm Start] Using ADC-derived initial D values (mean: %.2e mm²/s)\n', mean(adc_initial_guess));
                % Store initial guess in opts for IVIMmodelfit to use
                opts.d_initial = adc_initial_guess;
            end
        end
        
        if n_valid > 0
            fprintf('  [Stage 3 Opt] Flattening %d valid voxels for accelerated IVIM fit...\n', n_valid);

            % Extract 1D signal decay curves for valid voxels (reuse if already computed)
            if ~exist('dwi_valid', 'var')
                n_bvals = length(bvalues);
                dwi_valid = zeros(n_valid, n_bvals);
                for b = 1:n_bvals
                    vol_b = dwi(:,:,:,b);
                    dwi_valid(:,b) = vol_b(valid_voxels_idx);
                end
            end

            % [MODULARIZATION STAGE 3]: Masked 1D Flattening
            % The IVIMmodelfit dependency expects a 4D volume (x,y,z,b) as input.
            % We reshape our 1D valid-voxel array into a minimal 3D volume with
            % shape [N/2, 1, 2] to satisfy the dependency's dimensionality checks.
            % Padding to even length ensures the 3rd dimension is exactly 2;
            % MATLAB's size() would drop a trailing singleton dimension of 1,
            % causing the dependency to misinterpret the data layout.
            pad_len = mod(2 - mod(n_valid, 2), 2);
            dwi_valid_padded = [dwi_valid; zeros(pad_len, length(bvalues))];
            n_padded = n_valid + pad_len;

            % Reshape into a pseudo-3D volume for the dependency interface
            dwi_1d_vol = reshape(dwi_valid_padded, [n_padded/2, 1, 2, length(bvalues)]);
            mask_1d_vol = true(n_padded/2, 1, 2);
            if pad_len > 0
                mask_1d_vol(n_padded/2, 1, 2) = false;  % exclude padded voxel from fit
            end

            % Execute the segmented IVIM fit on the flattened masked array.
            % "seg" = segmented fitting strategy:
            %   Stage 1: High-b log-linear fit to isolate D (now with ADC warm start)
            %   Stage 2: Full-model fit with D fixed to estimate f and D*
            % This avoids the ill-conditioning of simultaneous 3-parameter
            % nonlinear fitting, which is especially problematic for D* because
            % the pseudo-diffusion signal decays rapidly and has low SNR.
            ivim_fit_1d = IVIMmodelfit(dwi_1d_vol, bvalues, "seg", mask_1d_vol, opts);

            % Restructure output back to strictly 1D and snip padding.
            % IVIMmodelfit returns a 4-column output: [D, S0, f, D*].
            % Column 2 (S0) is discarded — it is the signal at b=0, which is
            % already available from the raw DWI data and not a diffusion parameter.
            ivim_out_flat = reshape(ivim_fit_1d, [n_padded, 4]);
            d_vec = reshape(ivim_out_flat(1:n_valid, 1), [n_valid, 1]);
            f_vec = reshape(ivim_out_flat(1:n_valid, 3), [n_valid, 1]);
            dstar_vec = reshape(ivim_out_flat(1:n_valid, 4), [n_valid, 1]);

            % Release large intermediate arrays no longer needed to reduce
            % peak memory during parfor (each worker holds its own copy).
            clear dwi_valid dwi_valid_padded dwi_1d_vol mask_1d_vol ivim_fit_1d ivim_out_flat;

            % Replace zero-fit voxels with NaN (failed fits).
            % The segmented IVIM fitter returns D=0 when the log-linear
            % regression yields a non-positive slope (physically impossible for
            % diffusion). These represent voxels where noise dominates the
            % signal decay — common in low-SNR regions of pancreatic DWI.
            % Setting to NaN ensures they are excluded from nanmean-based
            % summary statistics rather than biasing the mean downward.
            zero_mask = (d_vec == 0);
            d_vec(zero_mask) = nan;
            f_vec(zero_mask) = nan;
            dstar_vec(zero_mask) = nan;
        end

        % Reconstruct 1D fitted parameters back into 3D volume geometry.
        % Background voxels (outside GTV mask) remain NaN, creating parameter
        % maps where only tumor voxels carry fitted values. This spatial
        % correspondence between mask and parameter maps is essential for
        % subsequent overlay visualization and voxel-level dose-response analysis.
        d_map = nan(sz3);
        f_map = nan(sz3);
        dstar_map = nan(sz3);

        if n_valid > 0
            % Place fitted 1D parameter values back into their original 3D
            % spatial positions using the linear indices from the GTV mask.
            d_map(valid_voxels_idx) = d_vec;
            f_map(valid_voxels_idx) = f_vec;
            dstar_map(valid_voxels_idx) = dstar_vec;
        end
    end

    % Restore original optimization options if they existed
    if ~isempty(original_opts)
        opts.optimoptions = original_opts;
    elseif isfield(opts, 'optimoptions')
        opts = rmfield(opts, 'optimoptions');
    end

    % Monoexponential ADC fit — weighted log-linear regression on active voxels.
    % [PHYSICS EXPLANATION]:
    % Apparent Diffusion Coefficient (ADC) assumes a simple mono-exponential decay:
    % S = S0 * exp(-b * ADC)
    % Linearized form: ln(S/S0) = -b * ADC
    % The log transform amplifies noise at low signal intensities (variance of
    % ln(S) ~ 1/S^2 for Gaussian noise).  Weighted Least Squares (WLS) with
    % weights = S^2 corrects for this heteroscedasticity, giving more weight
    % to high-SNR measurements.

    adc_sz  = [size(dwi,1), size(dwi,2), size(dwi,3)];  % same as sz3 above
    adc_map = nan(adc_sz);  % NaN background for non-tumor voxels

    % Determine whether to use GPU for the ADC WLS computation.
    % First verify a GPU actually exists via gpu_available, then check that
    % the device has enough free memory for the data transfer.
    use_gpu_adc = false;
    if isfield(opts, 'use_gpu') && opts.use_gpu
        gpu_device_idx = 1;
        if isfield(opts, 'gpu_device')
            gpu_device_idx = opts.gpu_device;
        end
        [gpu_ok, gpu_dev] = gpu_available(gpu_device_idx);
        if gpu_ok
            % Estimate GPU memory needed: signal matrix (n_valid x n_bvalues)
            % stored as double (8 bytes per element), plus b-value vector.
            estimated_bytes = n_valid * length(bvalues) * 8 + length(bvalues) * 8;
            if gpu_dev.AvailableMemory >= estimated_bytes
                use_gpu_adc = true;
            else
                fprintf('  [GPU] Insufficient GPU memory (need %.1f MB, available %.1f MB) — falling back to CPU.\n', ...
                    estimated_bytes / 1e6, gpu_dev.AvailableMemory / 1e6);
            end
        end
    end

    if n_valid > 0
        % Extract 1D signal decay curves if not already computed
        if ~exist('dwi_valid', 'var')
            n_bvals = length(bvalues);
            dwi_valid = zeros(n_valid, n_bvals);
            for b = 1:n_bvals
                vol_b = dwi(:,:,:,b);
                dwi_valid(:,b) = vol_b(valid_voxels_idx);
            end
        end

        % Filter to voxels with all-positive signal across ALL b-values.
        % The log-linear fit requires ln(S), so any zero or negative signal
        % intensity would produce -Inf or complex values. Zero signal at
        % high b-values occurs when tissue diffusivity is very high (signal
        % fully attenuated) or due to noise floor effects in magnitude DWI
        % images. These voxels cannot produce meaningful ADC estimates.
        adc_valid_idx = all(dwi_valid > 0, 2);

        if any(adc_valid_idx)
            S_a = dwi_valid(adc_valid_idx, :);

            % bvalues(1)==0 already validated at function entry.

            % Vectorized WLS: weights = S^2 at each b-value.
            %
            % MATHEMATICAL DERIVATION:
            % The monoexponential model S(b) = S0 * exp(-b * ADC) becomes
            % linear after taking logarithms: ln(S/S0) = -b * ADC.
            % However, ln() transforms Gaussian noise into heteroscedastic
            % noise with variance ~ 1/S^2 (via delta method).  Weighted
            % Least Squares with w = S^2 corrects for this heteroscedasticity,
            % giving more weight to high-SNR (high-signal) measurements and
            % reducing the influence of noisy high-b-value data points.
            %
            % The closed-form WLS solution for a single predictor is:
            %   ADC = sum(w * A * y) / sum(w * A^2)
            % where A = -b (design matrix), y = ln(S/S0), w = S^2.
            %
            % b=0 is excluded from the regression (used only as S0 reference)
            % because ln(S0/S0) = 0 provides no information about ADC.
            A_b = -bvalues(2:end);                             % [n_b x 1]

            if use_gpu_adc
                % GPU-accelerated path: transfer signal matrix and b-values
                % to GPU memory for vectorized WLS computation. gpuArray
                % operations use cuBLAS for element-wise and reduction ops,
                % which is beneficial when n_vox is large.
                try
                    S_a_gpu = gpuArray(S_a);
                    A_b_gpu = gpuArray(A_b);

                    Y_gpu = log(S_a_gpu(:,2:end) ./ S_a_gpu(:,1));
                    W_gpu = S_a_gpu(:,2:end).^2;

                    numer_gpu = sum(W_gpu .* Y_gpu .* A_b_gpu', 2);
                    denom_gpu = sum(W_gpu .* (A_b_gpu'.^2), 2);
                    adc_vals_gpu = numer_gpu ./ denom_gpu;

                    % Ensure all GPU operations complete before gathering
                    gpuDevice().wait();
                    adc_vals = gather(adc_vals_gpu);
                catch gpu_err
                    fprintf('  [GPU] ADC fitting GPU error: %s — falling back to CPU.\n', gpu_err.message);
                    use_gpu_adc = false;
                end
            end

            if ~use_gpu_adc
                % CPU path: standard vectorized WLS
                Y = log(S_a(:,2:end) ./ S_a(:,1));                % [n_vox x n_b]
                W = S_a(:,2:end).^2;                               % [n_vox x n_b]

                % Compute weighted numerator and denominator per voxel
                numer = sum(W .* Y .* A_b', 2);                   % [n_vox x 1]
                denom = sum(W .* (A_b'.^2), 2);                   % [n_vox x 1]
                adc_vals = numer ./ denom;
            end

            % Negative ADC values are physically impossible (diffusion cannot
            % be negative). They arise from noise-dominated signal where
            % signal increases with b-value (e.g., due to fat contamination,
            % T2 shine-through, or severe motion artifacts). Setting to NaN
            % excludes these from downstream analysis.
            adc_vals(adc_vals < 0 | ~isfinite(adc_vals)) = nan;

            % Prepare a temporary 1D vector for all valid mask voxels
            adc_vec_out = nan(n_valid, 1);
            adc_vec_out(adc_valid_idx) = adc_vals;

            % Reconstruct into 3D geometry
            adc_map(valid_voxels_idx) = adc_vec_out;
        end
    end
end