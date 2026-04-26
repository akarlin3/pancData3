function tf = has_table_var(T, varname)
% HAS_TABLE_VAR  Test whether a MATLAB table contains a named variable.
%
%   tf = HAS_TABLE_VAR(T, varname) returns true if VARNAME is one of the
%   variable (column) names of table T, false otherwise.
%
%   This is a workaround for the fact that ``isfield(T, name)`` returns
%   false for tables on MATLAB releases earlier than R2019b — even when
%   the column genuinely exists. Production data discovery (load_dwi_data
%   reading the clinical spreadsheet via readtable) silently dropped the
%   LF / Immuno / Pat lookups on the MSK MATLAB runtime, defaulting every
%   patient to ``pat_lf = 0`` regardless of the actual spreadsheet value.
%
%   Routing the existence check through ``ismember(varname,
%   T.Properties.VariableNames)`` works on every MATLAB version back to
%   R2013b (when tables were introduced) AND in Octave with the table
%   shim, so callers don't have to gate by version.
%
%   See also: load_dwi_data, isfield, ismember.

    if isempty(T)
        tf = false;
        return;
    end
    if istable(T)
        % Canonical accessor; works on every MATLAB release that
        % supports tables (R2013b+).
        try
            tf = ismember(varname, T.Properties.VariableNames);
            return;
        catch
            % Fall through to the secondary path on the off chance
            % Properties.VariableNames is unavailable (e.g. shimmed table).
        end
        try
            tf = ismember(varname, fieldnames(T));
        catch
            tf = false;
        end
    elseif isstruct(T)
        % Octave-fallback path in load_dwi_data builds T as a struct when
        % readtable can't be used; preserve isfield semantics there.
        tf = isfield(T, varname);
    else
        tf = false;
    end
end
