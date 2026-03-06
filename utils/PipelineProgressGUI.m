classdef PipelineProgressGUI < handle
%PIPELINEPROGRESSGUI Pipeline-aware progress bar wrapper around ProgressGUI.
%   Maps pipeline step keys (e.g., 'load', 'sanity') to human-readable
%   display names and tracks step-level progress through the DWI pipeline.
%
%   Usage:
%       gui = PipelineProgressGUI({'load','sanity','metrics_baseline'}, 'Standard');
%       gui.startStep('load');
%       % ... step runs ...
%       gui.completeStep('load', 'success');
%       gui.startStep('sanity');
%       % ...
%       gui.close();

    properties (Access = private)
        GUI                 % ProgressGUI instance
        StepKeys            % cell array of step keys being tracked
        StepDisplayNames    % cell array of display names (same order)
        TotalSteps          % number of tracked steps
        CompletedSteps      % number of steps finished
        DWITypeName         % e.g. 'Standard'
    end

    properties (Constant, Access = private)
        % Master mapping of step keys to display names
        STEP_MAP = { ...
            'load',                      'Load DWI Data'; ...
            'sanity',                    'Sanity Checks'; ...
            'metrics_baseline',          'Metrics: Baseline'; ...
            'metrics_longitudinal',      'Metrics: Longitudinal'; ...
            'metrics_dosimetry',         'Metrics: Dosimetry'; ...
            'metrics_stats_comparisons', 'Metrics: Comparisons'; ...
            'metrics_stats_predictive',  'Metrics: Predictive'; ...
            'visualize',                 'Visualization'; ...
            'metrics_survival',          'Metrics: Survival' ...
        };
    end

    methods
        function obj = PipelineProgressGUI(stepsToRun, dwiTypeName)
        %PIPELINEPROGRESSGUI Create a pipeline progress bar.
        %   stepsToRun  - cell array of step key strings
        %   dwiTypeName - string identifying the DWI type (e.g. 'Standard')
            obj.DWITypeName = dwiTypeName;
            obj.CompletedSteps = 0;

            % Filter to only the steps present in the master map
            % (excludes 'test' and config parse which are not tracked)
            allKeys = PipelineProgressGUI.STEP_MAP(:, 1);
            allNames = PipelineProgressGUI.STEP_MAP(:, 2);

            obj.StepKeys = {};
            obj.StepDisplayNames = {};
            for k = 1:numel(stepsToRun)
                idx = find(strcmp(stepsToRun{k}, allKeys), 1);
                if ~isempty(idx)
                    obj.StepKeys{end+1} = allKeys{idx};
                    obj.StepDisplayNames{end+1} = allNames{idx};
                end
            end
            obj.TotalSteps = numel(obj.StepKeys);

            if obj.TotalSteps == 0
                obj.GUI = [];
                return;
            end

            titleStr = sprintf('DWI Pipeline — %s', dwiTypeName);
            obj.GUI = ProgressGUI(titleStr, obj.TotalSteps);

            % Initial update
            counts = struct('completed', 0, 'total', obj.TotalSteps, ...
                            'stepName', sprintf('Step 0/%d', obj.TotalSteps));
            obj.GUI.update(0, counts, 'Initializing...', 'running');
        end

        function startStep(obj, stepKey)
        %STARTSTEP Signal the beginning of a pipeline step.
            if ~obj.isValid(), return; end

            idx = find(strcmp(stepKey, obj.StepKeys), 1);
            if isempty(idx), return; end

            displayName = obj.StepDisplayNames{idx};
            fraction = obj.CompletedSteps / max(obj.TotalSteps, 1);
            counts = struct('completed', obj.CompletedSteps, 'total', obj.TotalSteps, ...
                            'stepName', sprintf('Step %d/%d: %s', obj.CompletedSteps + 1, obj.TotalSteps, displayName));
            obj.GUI.update(fraction, counts, displayName, 'running');
        end

        function completeStep(obj, stepKey, status)
        %COMPLETESTEP Signal that a pipeline step has finished.
        %   status: 'success', 'warning', or 'skipped'
            if ~obj.isValid(), return; end

            idx = find(strcmp(stepKey, obj.StepKeys), 1);
            if isempty(idx), return; end

            obj.CompletedSteps = obj.CompletedSteps + 1;
            displayName = obj.StepDisplayNames{idx};
            fraction = obj.CompletedSteps / max(obj.TotalSteps, 1);

            if obj.CompletedSteps >= obj.TotalSteps
                stepLabel = 'All steps complete';
                guiStatus = 'success';
            else
                nextIdx = min(obj.CompletedSteps + 1, numel(obj.StepDisplayNames));
                nextName = obj.StepDisplayNames{nextIdx};
                stepLabel = sprintf('Step %d/%d: %s', obj.CompletedSteps + 1, obj.TotalSteps, nextName);
                if strcmp(status, 'warning')
                    guiStatus = 'failure';
                else
                    guiStatus = 'running';
                end
            end

            statusSuffix = '';
            if strcmp(status, 'warning')
                statusSuffix = ' (warning)';
            elseif strcmp(status, 'skipped')
                statusSuffix = ' (skipped)';
            end

            counts = struct('completed', obj.CompletedSteps, 'total', obj.TotalSteps, ...
                            'stepName', stepLabel);
            obj.GUI.update(fraction, counts, [displayName statusSuffix], guiStatus);
        end

        function close(obj)
        %CLOSE Close the pipeline progress GUI.
            if ~isempty(obj.GUI)
                obj.GUI.close();
            end
        end

        function valid = isValid(obj)
        %ISVALID Returns true if the underlying GUI is still valid.
            valid = ~isempty(obj.GUI) && obj.GUI.isValid();
        end
    end
end
