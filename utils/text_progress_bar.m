function text_progress_bar(current, total, label)
%TEXT_PROGRESS_BAR Displays an in-place text progress bar for long-running loops.
%
%   text_progress_bar(current, total, label)
%
%   Inputs:
%       current - Current iteration number (1-based)
%       total   - Total number of iterations
%       label   - Descriptive label string for the progress bar
%
%   Example:
%       for i = 1:100
%           text_progress_bar(i, 100, 'Processing');
%           % ... work ...
%       end

    if nargin < 3, label = 'Progress'; end
    if total <= 0, return; end

    bar_width = 30;
    fraction = current / total;
    filled = round(bar_width * fraction);
    empty_count = bar_width - filled;

    fill_char = char(9608);   % U+2588 FULL BLOCK
    empty_char = char(9617);  % U+2591 LIGHT SHADE

    bar_str = [repmat(fill_char, 1, filled), repmat(empty_char, 1, empty_count)];
    pct = floor(fraction * 100);

    if current >= total
        fprintf('\r  %s: |%s| %3d%% (%d/%d)\n', label, bar_str, pct, current, total);
    else
        fprintf('\r  %s: |%s| %3d%% (%d/%d)', label, bar_str, pct, current, total);
    end
end
