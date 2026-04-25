function val = parse_clinical_value(raw)
% PARSE_CLINICAL_VALUE  Coerce a clinical-spreadsheet cell to a numeric scalar.
%
%   val = PARSE_CLINICAL_VALUE(raw) takes whatever readtable returned for a
%   single binary outcome cell (e.g. T.LF(i_pat(1))) and returns a numeric
%   scalar in {0, 1, NaN}. Handles the type variants Excel + readtable can
%   produce when a column is technically numeric but the cell format is
%   "Text", or when blanks/categoricals sneak in:
%
%     numeric (double, single, int*, logical) → cast to double
%     char array  ('0', '1', '', '   ')        → str2double after strtrim
%     string scalar ("0", "1", "", missing)    → str2double after strtrim
%     cell {<one of the above>}                → unwrap and recurse
%     categorical                              → str2double(char(...))
%     missing / NaN / empty                    → NaN
%
%   The previous direct assignment ``pat_lf = T.LF(i_pat(1))`` worked when
%   T.LF was double but silently produced 0 / NaN / errors when readtable
%   inferred a non-numeric type (e.g. all 45 patients reading as lf=0
%   despite the spreadsheet's LF column having ~7 LF=1 entries). Routing
%   through this helper makes the read robust to those type variants.

    if nargin < 1
        val = NaN;
        return;
    end

    % Unwrap one level of cell; readtable returns 1x1 cells when the column
    % is a cell-of-char.
    if iscell(raw)
        if isempty(raw)
            val = NaN;
            return;
        end
        val = parse_clinical_value(raw{1});
        return;
    end

    % missing / empty
    if isempty(raw)
        val = NaN;
        return;
    end
    if isnumeric(raw) || islogical(raw)
        val = double(raw);
        if numel(val) ~= 1
            val = double(val(1));
        end
        return;
    end

    % string / char / categorical → strip whitespace, parse as number
    if isa(raw, 'string')
        if ismissing(raw)
            val = NaN;
            return;
        end
        s = strtrim(char(raw));
    elseif iscategorical(raw)
        if any(ismissing(raw))
            val = NaN;
            return;
        end
        s = strtrim(char(raw));
    elseif ischar(raw)
        s = strtrim(raw);
    else
        % Unknown type — best-effort coerce
        try
            s = strtrim(char(raw));
        catch
            val = NaN;
            return;
        end
    end

    if isempty(s)
        val = NaN;
        return;
    end
    val = str2double(s);
end
