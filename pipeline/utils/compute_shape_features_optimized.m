function [volume, surface_area, sphericity, elongation, compactness] = compute_shape_features_optimized(mask, voxel_spacing)
%COMPUTE_SHAPE_FEATURES_OPTIMIZED  Optimized shape features using regionprops3 and gradient-based surface area.
%
%   Uses built-in regionprops3 for volume calculation and vectorized gradient-based
%   method for surface area calculation to avoid O(n³) pixel iteration.

    % Initialize default values
    volume = NaN;
    surface_area = NaN;
    sphericity = NaN;
    elongation = NaN;
    compactness = NaN;

    % Check if mask is empty or all zeros
    if ~any(mask(:))
        return;
    end

    % Ensure mask is logical
    mask = logical(mask);

    % Handle 2D vs 3D cases
    is_3d = ndims(mask) == 3 && size(mask, 3) > 1;

    if is_3d
        dx = voxel_spacing(1);
        dy = voxel_spacing(2);
        dz = voxel_spacing(3);
        voxel_volume = dx * dy * dz;
    else
        dx = voxel_spacing(1);
        dy = voxel_spacing(2);
        voxel_area = dx * dy;
    end

    % Try using regionprops3 for 3D volumes or regionprops for 2D
    if is_3d && license('test', 'image_toolbox') && exist('regionprops3', 'file') == 2
        try
            % Use regionprops3 for efficient volume calculation
            props = regionprops3(mask, 'Volume');
            if ~isempty(props) && height(props) >= 1
                volume = props.Volume(1) * voxel_volume;
            else
                volume = sum(mask(:)) * voxel_volume;
            end

            % Compute surface area using gradient-based edge detection
            surface_area = compute_surface_area_gradient(mask, voxel_spacing);

        catch ME_regionprops3
            warning('compute_texture_features:regionprops3Failed', ...
                'regionprops3 failed (%s), using manual calculation', ME_regionprops3.message);
            % Fallback to manual calculation
            volume = sum(mask(:)) * voxel_volume;
            surface_area = compute_surface_area_gradient(mask, voxel_spacing);
        end
    elseif ~is_3d && license('test', 'image_toolbox') && exist('regionprops', 'file') == 2
        try
            % Use regionprops for 2D area calculation
            props = regionprops(mask, 'Area');
            if ~isempty(props)
                volume = props(1).Area * voxel_area; % 2D "volume" is area
            else
                volume = sum(mask(:)) * voxel_area;
            end

            % For 2D, "surface area" is perimeter
            perim_props = regionprops(mask, 'Perimeter');
            if ~isempty(perim_props)
                surface_area = perim_props(1).Perimeter * sqrt(dx * dy);
            else
                surface_area = compute_perimeter_2d(mask, dx, dy);
            end

        catch ME_regionprops
            warning('compute_texture_features:regionpropsFailed', ...
                'regionprops failed (%s), using manual calculation', ME_regionprops.message);
            % Fallback to manual calculation
            volume = sum(mask(:)) * voxel_area;
            surface_area = compute_perimeter_2d(mask, dx, dy);
        end
    else
        % Manual calculation for systems without Image Processing Toolbox
        if is_3d
            volume = sum(mask(:)) * voxel_volume;
            surface_area = compute_surface_area_gradient(mask, voxel_spacing);
        else
            volume = sum(mask(:)) * voxel_area;
            surface_area = compute_perimeter_2d(mask, dx, dy);
        end
    end

    % Calculate sphericity and compactness
    if ~isnan(volume) && ~isnan(surface_area) && volume > 0 && surface_area > 0
        if is_3d
            % 3D sphericity: (36π * V²)^(1/3) / A
            sphericity = ((36 * pi * volume^2)^(1/3)) / surface_area;
            % 3D compactness: A³ / (36π * V²)
            compactness = (surface_area^3) / (36 * pi * volume^2);
        else
            % 2D sphericity: 2√(π * A) / P
            sphericity = (2 * sqrt(pi * volume)) / surface_area;
            % 2D compactness: P² / (4π * A)
            compactness = (surface_area^2) / (4 * pi * volume);
        end
    end

    % Calculate elongation using covariance matrix eigenvalues
    elongation = compute_elongation_optimized(mask, voxel_spacing);
end
