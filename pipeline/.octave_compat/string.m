% STRING  Octave-compatible shim for MATLAB's string type.
%   Converts the input to a char array. MATLAB's string type is not
%   available in Octave, so this provides a basic workaround.
function s = string(x)
    if ischar(x)
        s = x;
    elseif iscell(x)
        s = x;
    elseif isnumeric(x)
        s = num2str(x);
    else
        s = char(x);
    end
end
