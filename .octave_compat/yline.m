% YLINE  Octave-compatible shim for MATLAB's yline (R2018b+).
%
%   hl = yline(yval) draws a horizontal constant line at y = yval.
%   hl = yline(yval, linespec, ...) passes additional arguments to plot().
%
%   MATLAB's yline creates a ConstantLine object that automatically adjusts
%   when axes limits change. This shim uses a simple plot() call spanning the
%   current x-axis limits, so the line will NOT auto-update if the axes are
%   later rescaled or panned.
%
%   Behavioral differences from MATLAB's yline:
%   - Returns a standard line handle, not a ConstantLine object.
%   - The line endpoints are fixed to the x-limits at the time of creation.
%   - Label argument (3rd positional in MATLAB's yline) is not supported;
%     all extra arguments are forwarded to plot() as line properties.
function hl = yline(yval, varargin)
    ax = gca();
    xlims = get(ax, 'XLim');
    % Save and restore hold state so this function does not alter it.
    hold_state = ishold(ax);
    hold(ax, 'on');

    % Draw a horizontal line spanning the current x-axis limits.
    hl = plot(xlims, [yval, yval], varargin{:});

    % Restore the original hold state.
    if ~hold_state
        hold(ax, 'off');
    end
end
