function result = nanmean_safe(v, dim)
% NANMEAN_SAFE — Octave-compatible NaN-ignoring mean.
% MATLAB's nanmean is in the Statistics Toolbox; Octave may lack it.
%
% Inputs:
%   v   - numeric array
%   dim - (optional) dimension along which to take the mean.  When
%         omitted, behaves as nanmean(v) (single-arg form: collapses
%         to a scalar over all non-NaN elements).
%
% Output:
%   result - mean of non-NaN elements (NaN where the slice has no
%            non-NaN entries).
if nargin < 2
    if exist('OCTAVE_VERSION', 'builtin')
        tmp = v(~isnan(v));
        if isempty(tmp)
            result = NaN;
        else
            result = mean(tmp);
        end
    else
        result = nanmean(v);
    end
else
    if exist('OCTAVE_VERSION', 'builtin')
        result = mean(v, dim, 'omitnan');
    else
        result = nanmean(v, dim);
    end
end
end
