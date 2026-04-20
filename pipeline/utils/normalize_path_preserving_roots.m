function out = normalize_path_preserving_roots(p)
% NORMALIZE_PATH_PRESERVING_ROOTS  Normalize separators without losing roots.
%
%   Out = normalize_path_preserving_roots(p) converts '/' and '\' in P to
%   the native filesep (fullfile semantics) while preserving:
%     - Windows UNC prefix '\\server\share\...'
%     - Windows drive letter 'C:\...'
%     - Unix absolute path '/abs/path'
%
%   The naive pattern
%       fullfile(strsplit(p, {'/','\'}){:})
%   strips empty tokens, which eats the '\\' of a UNC path and the leading
%   '/' of a Unix absolute path, yielding unresolvable relative paths.

    if isempty(p)
        out = p;
        return;
    end

    % Detect Windows UNC prefix: leading '\\' or '//'.
    is_unc = (length(p) >= 2) && any(strcmp(p(1:2), {'\\', '//'}));
    % Detect Unix absolute path.
    is_unix_abs = p(1) == '/' && ~is_unc;

    path_parts = strsplit(p, {'/', '\'});
    % strsplit on leading separator produces an empty first token — drop
    % the empties we are about to re-add as an explicit root.
    non_empty = path_parts(~cellfun('isempty', path_parts));

    if isempty(non_empty)
        out = p;
        return;
    end

    if is_unc
        % Rebuild '\\server\share\rest...' with native filesep.
        out = [filesep filesep strjoin(non_empty, filesep)];
    elseif is_unix_abs
        out = [filesep strjoin(non_empty, filesep)];
    else
        out = fullfile(non_empty{:});
    end
end
