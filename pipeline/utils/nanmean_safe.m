function result = nanmean_safe(v)
% NANMEAN_SAFE — Octave-compatible NaN-ignoring mean.
% MATLAB's nanmean is in the Statistics Toolbox; Octave may lack it.
%
% Input:
%   v - numeric vector
%
% Output:
%   result - mean of non-NaN elements (NaN if all elements are NaN)
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
end
