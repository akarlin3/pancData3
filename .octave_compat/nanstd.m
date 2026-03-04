function y = nanstd(x, flag, dim)
    if nargin < 2 || isempty(flag)
        flag = 0;
    end
    if nargin < 3
        dim = find(size(x) ~= 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end

    mask = ~isnan(x);
    n = sum(mask, dim);

    x(isnan(x)) = 0;
    mean_x = sum(x, dim) ./ n;

    % Broadcast mean_x back to size of x for element-wise subtraction
    sz = size(x);
    rep_dims = ones(1, ndims(x));
    rep_dims(dim) = sz(dim);
    mean_x_rep = repmat(mean_x, rep_dims);

    diff = x - mean_x_rep;
    diff(~mask) = 0;

    if flag == 0
        denom = max(n - 1, 1);
    else
        denom = n;
    end

    y = sqrt(sum(diff.^2, dim) ./ denom);
    y(n == 0) = NaN;
end
