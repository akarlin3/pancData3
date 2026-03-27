% NANMEDIAN  NaN-ignoring median (Octave-compatible shim).
%   M = nanmedian(X) returns the median of X, ignoring NaN values.
%   M = nanmedian(X, dim) operates along dimension dim.
function m = nanmedian(x, dim)
    if nargin < 2
        if isvector(x)
            x = x(~isnan(x));
            m = median(x);
        else
            dim = 1;
            m = nanmedian(x, dim);
        end
    else
        nd = ndims(x);
        sz = size(x);
        % Move target dimension to first position
        perm = [dim, 1:dim-1, dim+1:nd];
        x = permute(x, perm);
        sz_p = size(x);
        x = reshape(x, sz_p(1), []);
        m = zeros(1, size(x, 2));
        for i = 1:size(x, 2)
            col = x(:, i);
            col = col(~isnan(col));
            if isempty(col)
                m(i) = NaN;
            else
                m(i) = median(col);
            end
        end
        % Reshape back
        out_sz = sz_p;
        out_sz(1) = 1;
        m = reshape(m, out_sz);
        % Undo permutation
        m = ipermute(m, perm);
    end
end
