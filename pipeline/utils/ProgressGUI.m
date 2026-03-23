classdef ProgressGUI < handle
%PROGRESSGUI Professional custom-figure progress bar for MATLAB pipelines.
%   Creates a styled figure window with a colored progress bar, percentage
%   display, detail text, summary, and elapsed time. Designed to replace
%   MATLAB's built-in waitbar() with a more polished look.
%
%   The figure uses HandleVisibility='off' so it survives
%   set(0,'DefaultFigureVisible','off') used by the pipeline.
%
%   Usage:
%       gui = ProgressGUI('Running Tests', 60);
%       gui.update(0.5, struct('completed',30,'total',60,'passed',28,'failed',2), ...
%                  'test_foo/test_bar', 'running');
%       gui.close();

    properties (Access = private)
        Figure              % figure handle (HandleVisibility='off' to survive global visibility changes)
        BarAxes             % axes containing the progress bar patches
        BarBackground       % patch: full-width gray track (the "empty" portion)
        BarFill             % patch: colored fill that grows left-to-right with progress
        TitleText           % uicontrol: bold heading at top of figure
        PercentText         % uicontrol: large percentage display (e.g. "45.0%")
        CountText           % uicontrol: count display beside percentage (e.g. "(9/20)")
        DetailText          % uicontrol: current item name or step description
        SummaryText         % uicontrol: pass/fail counts or step-level summary (bottom-left)
        ElapsedText         % uicontrol: wall-clock elapsed time (bottom-right)
        StartTime           % tic token captured at construction for elapsed time calculation
    end

    properties (Constant, Access = private)
        % Material Design-inspired color palette for status indication
        COLOR_GREEN  = [0.30, 0.69, 0.31]   % success / running normally
        COLOR_AMBER  = [1.00, 0.76, 0.03]   % partial failures (some tests failed)
        COLOR_RED    = [0.90, 0.22, 0.21]   % failure / critical error
        COLOR_BG     = [0.97, 0.97, 0.97]   % figure background (near-white)
        COLOR_TRACK  = [0.88, 0.88, 0.88]   % unfilled progress bar track
        COLOR_GRAY   = [0.45, 0.45, 0.45]   % secondary text (counts, elapsed time)
    end

    methods
        function obj = ProgressGUI(titleStr, total, options)
        %PROGRESSGUI Create a professional progress bar figure.
        %   obj = ProgressGUI(title, total)
        %   obj = ProgressGUI(title, total, options)
        %
        %   options struct fields (all optional):
        %       .width   - figure width in pixels (default 480)
        %       .height  - figure height in pixels (default 200)
            if nargin < 3, options = struct(); end
            figW = 480; if isfield(options, 'width'),  figW = options.width;  end
            figH = 200; if isfield(options, 'height'), figH = options.height; end

            obj.StartTime = tic;

            % Center on screen
            screenSz = get(0, 'ScreenSize');
            figX = (screenSz(3) - figW) / 2;
            figY = (screenSz(4) - figH) / 2;

            % Create a minimal figure with no menu/toolbar.
            % HandleVisibility='off' prevents this figure from being
            % affected by set(0,'DefaultFigureVisible','off') which the
            % pipeline uses to suppress plot windows during batch runs.
            % CloseRequestFcn hides rather than deletes so re-opening
            % is possible if the user accidentally closes the window.
            obj.Figure = figure( ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off', ...
                'Name', titleStr, ...
                'Resize', 'off', ...
                'Color', ProgressGUI.COLOR_BG, ...
                'Position', [figX, figY, figW, figH], ...
                'HandleVisibility', 'off', ...
                'Visible', 'on', ...
                'CloseRequestFcn', @(~,~) set(gcbo, 'Visible', 'off'));

            % --- Title ---
            obj.TitleText = uicontrol(obj.Figure, 'Style', 'text', ...
                'String', titleStr, ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.82, 0.90, 0.14], ...
                'FontSize', 13, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'BackgroundColor', ProgressGUI.COLOR_BG);

            % --- Percentage ---
            obj.PercentText = uicontrol(obj.Figure, 'Style', 'text', ...
                'String', '0.0%', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.62, 0.50, 0.18], ...
                'FontSize', 18, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'BackgroundColor', ProgressGUI.COLOR_BG);

            % --- Count ---
            countStr = '';
            if total > 0
                countStr = sprintf('(0/%d)', total);
            end
            obj.CountText = uicontrol(obj.Figure, 'Style', 'text', ...
                'String', countStr, ...
                'Units', 'normalized', ...
                'Position', [0.55, 0.65, 0.40, 0.12], ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'center', ...
                'ForegroundColor', ProgressGUI.COLOR_GRAY, ...
                'BackgroundColor', ProgressGUI.COLOR_BG);

            % --- Progress bar (axes + patches) ---
            % Uses a hidden axes with [0,1] x [0,1] coordinate space.
            % Two overlapping patches: background track (full width) and
            % colored fill (width updated via XData in update()).
            obj.BarAxes = axes('Parent', obj.Figure, ...
                'Units', 'normalized', ...
                'Position', [0.08, 0.42, 0.84, 0.14], ...
                'XLim', [0 1], 'YLim', [0 1], ...
                'XTick', [], 'YTick', [], ...
                'Box', 'off', ...
                'Color', 'none', ...
                'XColor', 'none', 'YColor', 'none', ...
                'HandleVisibility', 'off');

            % Track background: full-width gray rectangle
            obj.BarBackground = patch(obj.BarAxes, ...
                [0, 0, 1, 1], [0, 1, 1, 0], ...
                ProgressGUI.COLOR_TRACK, ...
                'EdgeColor', [0.80, 0.80, 0.80], ...
                'LineWidth', 0.5);

            % Fill patch: starts at zero width; XData(3:4) are updated to
            % the current fraction in update() to animate the bar
            obj.BarFill = patch(obj.BarAxes, ...
                [0, 0, 0, 0], [0, 1, 1, 0], ...
                ProgressGUI.COLOR_GREEN, ...
                'EdgeColor', 'none');

            % --- Detail text ---
            obj.DetailText = uicontrol(obj.Figure, 'Style', 'text', ...
                'String', '', ...
                'Units', 'normalized', ...
                'Position', [0.08, 0.26, 0.84, 0.12], ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'center', ...
                'ForegroundColor', ProgressGUI.COLOR_GRAY, ...
                'BackgroundColor', ProgressGUI.COLOR_BG);

            % --- Summary text (bottom-left) ---
            obj.SummaryText = uicontrol(obj.Figure, 'Style', 'text', ...
                'String', '', ...
                'Units', 'normalized', ...
                'Position', [0.08, 0.04, 0.55, 0.18], ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', ProgressGUI.COLOR_BG);

            % --- Elapsed time (bottom-right) ---
            obj.ElapsedText = uicontrol(obj.Figure, 'Style', 'text', ...
                'String', 'Elapsed: 0s', ...
                'Units', 'normalized', ...
                'Position', [0.60, 0.04, 0.35, 0.18], ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'right', ...
                'ForegroundColor', ProgressGUI.COLOR_GRAY, ...
                'BackgroundColor', ProgressGUI.COLOR_BG);

            % 'limitrate' caps redraws at ~20 fps to avoid GUI overhead
            if exist('OCTAVE_VERSION', 'builtin')
                drawnow();
            else
                drawnow('limitrate');
            end
        end

        function update(obj, fraction, counts, detail, status)
        %UPDATE Refresh the progress bar display.
        %   update(fraction, counts, detail, status)
        %
        %   fraction - 0.0 to 1.0
        %   counts   - struct with optional fields: completed, total, passed, failed
        %   detail   - string for the detail line
        %   status   - 'running' | 'success' | 'failure'
            if ~obj.isValid(), return; end

            % Clamp fraction to [0, 1] to prevent patch drawing errors
            fraction = max(0, min(1, fraction));

            % Animate the bar by setting the right edge of the fill patch
            % XData = [left-bottom, left-top, right-top, right-bottom]
            set(obj.BarFill, 'XData', [0, 0, fraction, fraction]);

            % Update bar color based on status
            if nargin >= 5
                barColor = obj.getBarColor(counts, status);
                set(obj.BarFill, 'FaceColor', barColor);
            end

            % Percentage
            set(obj.PercentText, 'String', sprintf('%.1f%%', fraction * 100));

            % Count
            if isstruct(counts) && isfield(counts, 'completed') && isfield(counts, 'total')
                set(obj.CountText, 'String', sprintf('(%d/%d)', counts.completed, counts.total));
            end

            % Detail
            if nargin >= 4 && ~isempty(detail)
                % Truncate long detail strings
                if length(detail) > 70
                    detail = ['...' detail(end-66:end)];
                end
                set(obj.DetailText, 'String', detail);
            end

            % Summary
            if isstruct(counts)
                summaryStr = obj.buildSummary(counts);
                set(obj.SummaryText, 'String', summaryStr);
            end

            % Elapsed time
            elapsed = toc(obj.StartTime);
            set(obj.ElapsedText, 'String', ['Elapsed: ' ProgressGUI.formatTime(elapsed)]);

            if exist('OCTAVE_VERSION', 'builtin')
                drawnow();
            else
                drawnow('limitrate');
            end
        end

        function setDetail(obj, text)
        %SETDETAIL Update only the detail line.
            if ~obj.isValid(), return; end
            if length(text) > 70
                text = ['...' text(end-66:end)];
            end
            set(obj.DetailText, 'String', text);
            if exist('OCTAVE_VERSION', 'builtin')
                drawnow();
            else
                drawnow('limitrate');
            end
        end

        function close(obj)
        %CLOSE Close the progress figure if it exists.
            try
                if ~isempty(obj.Figure) && isvalid(obj.Figure)
                    delete(obj.Figure);
                end
            catch
                % Silently ignore
            end
        end

        function valid = isValid(obj)
        %ISVALID Returns true if the figure still exists and is valid.
            valid = ~isempty(obj.Figure) && isvalid(obj.Figure);
        end
    end

    methods (Static)
        function available = isDisplayAvailable()
        %ISDISPLAYAVAILABLE Check if a GUI display is available.
        %   Returns false in Octave (no Java/figure support in headless mode),
        %   in headless MATLAB (no DISPLAY on Linux), or when JVM is absent.
            available = false;
            if exist('OCTAVE_VERSION', 'builtin'), return; end
            try
                % usejava('desktop') is true for interactive MATLAB sessions;
                % for -nodisplay/-batch on Linux, check DISPLAY + JVM instead
                available = usejava('desktop') || ...
                            (~isempty(getenv('DISPLAY')) && usejava('jvm'));
            catch
                available = false;
            end
        end
    end

    methods (Access = private)
        function color = getBarColor(~, counts, status)
        %GETBARCOLOR Determine bar fill color from status and failure count.
        %   Red = explicit failure, Amber = partial failures (some tests
        %   failed but run continues), Green = success or running normally.
            if strcmp(status, 'failure')
                color = ProgressGUI.COLOR_RED;
            elseif strcmp(status, 'success')
                color = ProgressGUI.COLOR_GREEN;
            elseif isstruct(counts) && isfield(counts, 'failed') && counts.failed > 0
                color = ProgressGUI.COLOR_AMBER;
            else
                color = ProgressGUI.COLOR_GREEN;
            end
        end

        function str = buildSummary(~, counts)
        %BUILDSUMMARY Build summary text from counts struct.
            str = '';
            if isfield(counts, 'passed') && isfield(counts, 'failed')
                if counts.failed > 0
                    str = sprintf('%d passed, %d failed', counts.passed, counts.failed);
                else
                    str = sprintf('%d passed', counts.passed);
                end
            elseif isfield(counts, 'completed') && isfield(counts, 'total')
                str = sprintf('%d / %d complete', counts.completed, counts.total);
            end
            if isfield(counts, 'stepName') && ~isempty(counts.stepName)
                str = counts.stepName;
            end
        end
    end

    methods (Static, Access = private)
        function str = formatTime(seconds)
        %FORMATTIME Format elapsed seconds as human-readable string.
            if seconds < 60
                str = sprintf('%.0fs', seconds);
            elseif seconds < 3600
                m = floor(seconds / 60);
                s = round(mod(seconds, 60));
                str = sprintf('%dm %ds', m, s);
            else
                h = floor(seconds / 3600);
                m = floor(mod(seconds, 3600) / 60);
                s = round(mod(seconds, 60));
                str = sprintf('%dh %dm %ds', h, m, s);
            end
        end
    end
end
