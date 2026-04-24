function fxtmp_idx = match_fraction_folder(names, fx_tag)
% MATCH_FRACTION_FOLDER  Logical mask of folder names that represent a given fraction.
%
%   fxtmp_idx = MATCH_FRACTION_FOLDER(names, fx_tag) returns a logical row
%   vector the same size as the cell array of folder NAMES, with true
%   entries where the folder name starts with FX_TAG (case-insensitive)
%   AND either ends there or is followed by a non-alphanumeric character.
%
%   Examples (fx_tag = 'Fx1'):
%       'Fx1'                  -> true   (exact match)
%       'Fx1 - repeatability'  -> true   (space is non-alphanumeric)
%       'fx1_old'              -> true   (underscore is non-alphanumeric)
%       'Fx10'                 -> false  ('0' is alphanumeric)
%       'Fx11'                 -> false  ('1' is alphanumeric)
%       'xFx1'                 -> false  (prefix does not start at index 1)
%
%   Examples (fx_tag = 'post'):
%       'post'                 -> true
%       'post - treatment'     -> true
%       'postprocessing'       -> false  ('p' is alphanumeric)
%
% Rationale
%   Earlier versions used regex ('Fx1(\b|$)' and later '^Fx1([^A-Za-z0-9]|$)')
%   but both silently failed against 'Fx1 - repeatability' in some MATLAB
%   runtimes on the MSK cohort, leaving gtv_locations{:, 1, :} empty for
%   ~95% of patients and cascading into all-NaN dice_rpt_*. strncmpi plus
%   an explicit next-character check removes any dependency on regex-engine
%   semantics for word boundaries / anchors.

    prefix_len = numel(fx_tag);
    % Case-insensitive prefix check across the whole cell array in one shot.
    prefix_hit = strncmpi(names, fx_tag, prefix_len);

    fxtmp_idx = false(size(prefix_hit));
    for k = find(prefix_hit)
        nm = names{k};
        if numel(nm) == prefix_len
            fxtmp_idx(k) = true;  % exact match: 'Fx1' / 'post'
        else
            next_char = nm(prefix_len + 1);
            % Accept any non-alphanumeric separator (space, hyphen,
            % underscore, dot, etc.). Reject alphanumeric continuations
            % so 'Fx10'/'Fx11' do not collide with 'Fx1' and
            % 'postprocessing' does not collide with 'post'.
            if ~(isstrprop(next_char, 'alpha') || isstrprop(next_char, 'digit'))
                fxtmp_idx(k) = true;
            end
        end
    end
end
