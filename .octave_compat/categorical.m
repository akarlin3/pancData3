% CATEGORICAL  Octave-compatible stub for MATLAB's categorical type.
%
%   MATLAB's categorical array (R2013b+) stores data as a set of discrete
%   categories with an underlying integer encoding. Octave does not have a
%   native categorical type.
%
%   This stub simply returns the input unchanged, which works when the
%   downstream code only uses the categorical for display or grouping via
%   unique(). It does NOT provide category ordering, merging, or relabeling.
function c = categorical(x)
    % Pass-through: return the input as-is (cell array of strings, numeric, etc.).
    c = x;
end
