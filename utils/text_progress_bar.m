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
%
% --- Analytical Rationale ---
% The DWI pipeline processes patient cohorts that can take hours to complete.
% Key bottlenecks include: DICOM-to-NIfTI conversion via dcm2niix, IVIM
% bi-exponential fitting (iterative nonlinear least squares per voxel),
% DnCNN inference, and deformable image registration. A visual progress
% indicator is essential for the researcher to estimate completion time and
% detect stalled processing early.
%
% This function uses carriage return ('\r') to overwrite the current line
% in-place, avoiding log flooding. When the pipeline's diary logging is
% active, each '\r' update does NOT create a new line in the log file on
% most platforms, keeping log files compact. A newline is only emitted on
% the final iteration to cleanly terminate the progress bar.

    if nargin < 3, label = 'Progress'; end
    if total <= 0, return; end

    % --- Bar Rendering ---
    % A 30-character-wide bar provides sufficient resolution (~3.3% per
    % character) while fitting within an 80-column terminal alongside the
    % label, percentage, and count annotation.
    bar_width = 30;
    fraction = current / total;
    filled = round(bar_width * fraction);
    empty_count = bar_width - filled;

    % Unicode block characters provide a clean visual appearance in modern
    % terminals and the MATLAB command window. U+2588 (FULL BLOCK) for
    % completed portion, U+2591 (LIGHT SHADE) for remaining.
    fill_char = char(9608);   % U+2588 FULL BLOCK
    empty_char = char(9617);  % U+2591 LIGHT SHADE

    bar_str = [repmat(fill_char, 1, filled), repmat(empty_char, 1, empty_count)];
    pct = fraction * 100;

    if current >= total
        % Final iteration: emit a newline so subsequent console output
        % starts on a fresh line rather than overwriting the progress bar.
        fprintf('\r  %s: |%s| %5.1f%% (%d/%d)\n', label, bar_str, pct, current, total);
    else
        % Intermediate iterations: carriage return without newline
        % overwrites the previous bar in-place, creating an animated effect.
        fprintf('\r  %s: |%s| %5.1f%% (%d/%d)', label, bar_str, pct, current, total);
    end
end
