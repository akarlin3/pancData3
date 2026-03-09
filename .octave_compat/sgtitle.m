% SGTITLE  Octave-compatible shim for MATLAB's sgtitle (R2018b+).
%
%   sgtitle(txt) adds a title to an entire figure containing subplots.
%   MATLAB's sgtitle places text centered above all subplot axes. Octave
%   does not have sgtitle, so this shim falls back to title(), which adds
%   the text to the current axes instead of the figure as a whole.
%
%   Behavioral differences from MATLAB's sgtitle:
%   - Text is placed on the current axes, not as a figure-level annotation.
%   - Additional name-value formatting arguments (FontSize, etc.) are accepted
%     via varargin but passed to title(), which may not support all of them.
function sgtitle(txt, varargin)
    % Fallback: apply as a regular axes title.
    title(txt);
end
