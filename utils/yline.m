function hl = yline(yval, varargin)
    % YLINE Create a horizontal line
    ax = gca();
    xlims = get(ax, 'XLim');
    hold_state = ishold(ax);
    hold(ax, 'on');

    % Basic implementation just to plot the line
    hl = plot(xlims, [yval, yval], varargin{:});

    if ~hold_state
        hold(ax, 'off');
    end
end
