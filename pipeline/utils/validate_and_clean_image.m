function img_clean = validate_and_clean_image(img, mask)
%VALIDATE_AND_CLEAN_IMAGE  Remove NaN/Inf values and handle negative values for 2D image.
    img_clean = double(img);

    % Replace NaN/Inf values outside mask with 0
    invalid_mask = ~isfinite(img_clean);
    img_clean(invalid_mask & ~mask) = 0;

    % For values inside mask, replace NaN/Inf with median of valid values
    if any(invalid_mask & mask)
        valid_vals = img_clean(mask & isfinite(img_clean));
        if ~isempty(valid_vals)
            replacement_val = median(valid_vals);
        else
            replacement_val = 0;
        end
        img_clean(invalid_mask & mask) = replacement_val;
    end

    % Handle negative values by taking absolute value for texture calculations
    % (preserves original values for first-order features computed earlier)
    img_clean = abs(img_clean);
end
