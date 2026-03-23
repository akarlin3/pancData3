% ISCATEGORICAL  Octave-compatible shim for MATLAB's iscategorical.
%   Returns false since Octave does not have a native categorical type.
function tf = iscategorical(x)
    tf = false;
end
