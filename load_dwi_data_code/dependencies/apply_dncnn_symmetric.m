function img_denoised = apply_dncnn_symmetric(img, net, pad_size)
% apply_dncnn_symmetric - Applies spatial dnCNN denoising with symmetric padding 
% to mitigate edge-artifacts at the tissue margins (e.g., GTV boundaries).
%
% Zero-padding pulls background values into the convolutional kernels,
% artificially suppressing signals at the boundaries and corrupting dose 
% metrics like D95. Symmetric reflection mitigates this.
%
% Inputs:
%   img      - The original 2D or 3D image array before GTV masking
%   net      - The trained DnCNN dlnetwork or SeriesNetwork
%   pad_size - (Optional) The padding size. Should match or exceed 
%              the largest receptive field of the network. For a standard 
%              17-layer DnCNN, the receptive field is 35x35 (pad=17).
%              Default is 17 if not provided.
%
% Outputs:
%   img_denoised - The denoised image, perfectly cropped back to original 
%                  matrix dimensions.

    if nargin < 3
        % Default pad size for a standard 17-layer DnCNN
        pad_size = 17; 
    end

    % 1. Get original dimensions for cropping
    orig_size = size(img);
    
    % Handle 2D vs 3D cases for padding
    if ndims(img) == 2
        pad_dims = [pad_size, pad_size];
    elseif ndims(img) == 3
        % Assuming 2D spatial convolution is applied slice-by-slice,
        % or 3D spatial convolution. If 3D conv, we pad in all 3 dims.
        % We will conservatively pad all 3 spatial dimensions here, 
        % or you can adjust to [pad_size, pad_size, 0] for 2D slice-by-slice.
        pad_dims = [pad_size, pad_size, pad_size];
    else
        error('Input image must be 2D or 3D.');
    end

    % 2. Pad the array using symmetric reflection
    % This reflects the edge pixels outward, preventing boundary zeros.
    img_padded = padarray(img, pad_dims, 'symmetric');
    
    % 3. Format for the network (Deep Learning Toolbox usually expects 
    % specific formatting like 'SSCB' or numeric arrays for predict)
    % Depending on how your pipeline calls the network, you might use predict()
    % or forward(). For standard predict with a formatted dlarray:
    
    % (Note: If you use the standard predict(net, img_padded), uncomment below)
    % img_denoised_padded = predict(net, img_padded);
    
    % If using MATLAB's built-in denoiseImage (which uses a pretrained DnCNN):
    % img_denoised_padded = denoiseImage(img_padded, net);
    
    % Placeholder for inference (replace with your exact inference call):
    if isa(net, 'dlnetwork') || isa(net, 'SeriesNetwork') || isa(net, 'DAGNetwork')
        % Extract underlying data if it returns a dlarray
        img_denoised_padded = predict(net, single(img_padded));
        if isa(img_denoised_padded, 'dlarray')
            img_denoised_padded = extractdata(img_denoised_padded);
        end
    else
        % Fallback for built-in or custom wrappers
        img_denoised_padded = denoiseImage(img_padded, net);
    end
    
    % 4. Precisely crop the output back to original dimensions
    % Crop syntax: img( start : end, start : end, ... )
    if ndims(img) == 2
        img_denoised = img_denoised_padded(...
            pad_dims(1) + 1 : pad_dims(1) + orig_size(1), ...
            pad_dims(2) + 1 : pad_dims(2) + orig_size(2));
    elseif ndims(img) == 3
        img_denoised = img_denoised_padded(...
            pad_dims(1) + 1 : pad_dims(1) + orig_size(1), ...
            pad_dims(2) + 1 : pad_dims(2) + orig_size(2), ...
            pad_dims(3) + 1 : pad_dims(3) + orig_size(3));
    end

    % The resulting img_denoised can now be safely masked by the GTV 
    % without suppressing boundary diffusion metrics.
end
