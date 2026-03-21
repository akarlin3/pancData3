function result = nanstd_safe(v)
% NANSTD_SAFE — Octave-compatible NaN-ignoring standard deviation.
%
% Input:
%   v - numeric vector
%
% Output:
%   result - standard deviation of non-NaN elements (NaN if all NaN)
if exist('OCTAVE_VERSION', 'builtin')
    tmp = v(~isnan(v));
    if isempty(tmp)
        result = NaN;
    else
        result = std(tmp);
    end
else
    result = nanstd(v);
end
end
