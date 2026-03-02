function y = nanmean(x, dim)
    if nargin == 1
        dim = find(size(x) ~= 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end
    mask = ~isnan(x);
    n = sum(mask, dim);
    x(isnan(x)) = 0;
    y = sum(x, dim) ./ n;
    y(n == 0) = NaN;
end
