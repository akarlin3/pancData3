function [dwi_dncnn, havedenoised] = compute_dncnn_fallback(dwi, i_sort, gtv_mask, gtvn_mask, use_gpu)
    % COMPUTE_DNCNN_FALLBACK On-the-fly DnCNN denoising when cache is missing
    %   Extracted local helper from process_single_scan.m. Loads the DnCNN
    %   model and denoises each b-value volume of DWI (optionally on GPU),
    %   returning the normalised denoised volume and a success flag. Failures
    %   are caught and reported, leaving HAVEDENOISED at 0.
    if nargin < 5, use_gpu = false; end
    havedenoised = 0;
    dwi_dncnn = [];
    gpu_active = false;
    if use_gpu
        [gpu_ok, ~] = gpu_available();
        if gpu_ok
            gpu_active = true;
            fprintf('  [DnCNN] Cache missing. Executing deep learning denoising on GPU...\n');
        else
            fprintf('  [DnCNN] GPU requested but unavailable. Falling back to CPU.\n');
        end
    end
    if ~gpu_active
        fprintf('  [DnCNN] Cache missing. Executing deep learning denoising on CPU...\n');
    end
    try
        loaded_model = load(fullfile(fileparts(mfilename('fullpath')), '..', 'dependencies', 'dncnn_model.mat'), 'net');
        dncnn_net = loaded_model.net;
        if gpu_active && isa(dncnn_net, 'dlnetwork')
            try
                dncnn_net = dlupdate(@gpuArray, dncnn_net);
                fprintf('  [DnCNN] Network moved to GPU.\n');
            catch
                fprintf('  [DnCNN] Could not move network to GPU. Using CPU.\n');
                gpu_active = false;
            end
        end
        dwi_input = single(dwi);
        if ~isempty(gtv_mask)
            mask_input = single(gtv_mask);
        elseif ~isempty(gtvn_mask)
            mask_input = single(gtvn_mask);
        else
            mask_input = ones(size(dwi, 1), size(dwi, 2), size(dwi, 3), 'single');
        end
        dwi_dncnn_out = zeros(size(dwi_input), 'single');
        n_bvals_dncnn = size(dwi_input, 4);
        for b_idx = 1:n_bvals_dncnn
            text_progress_bar(b_idx, n_bvals_dncnn, 'DnCNN denoising b-values');
            dwi_dncnn_out(:,:,:,b_idx) = apply_dncnn_symmetric(dwi_input(:,:,:,b_idx), mask_input, dncnn_net, 15);
        end
        dwi_dncnn = double(mat2gray(dwi_dncnn_out));
        havedenoised = 1;
        device_label = 'GPU';
        if ~gpu_active, device_label = 'CPU'; end
        fprintf('  [DnCNN] Deep learning denoising completed on %s.\n', device_label);
    catch ME
        fprintf('  [DnCNN] Computation failed: %s\n', ME.message);
        fprintf('  [DnCNN] Proceeding without denoised data.\n');
    end
end
