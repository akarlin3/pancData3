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
        GUI                 % ProgressGUI instance (the underlying figure-based progress bar)
        StepKeys            % cell array of step keys being tracked (e.g. {'load','sanity'})
        StepDisplayNames    % cell array of human-readable display names (same order as StepKeys)
        TotalSteps          % number of tracked steps (excludes unrecognized keys)
        CompletedSteps      % running count of steps finished so far
        DWITypeName         % DWI processing type label, e.g. 'Standard', 'dnCNN', 'IVIMnet'
    end

    properties (Constant, Access = private)
        % Master mapping of step keys to display names
        STEP_MAP = { ...
            'load',                      'Load DWI Data'; ...
            'sanity',                    'Sanity Checks'; ...
            'metrics_baseline',          'Metrics: Baseline'; ...
            'compare_cores',             'Compare Core Methods'; ...
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

            % Filter stepsToRun to only those present in the master STEP_MAP.
            % Internal steps like 'test' and config parsing are not tracked
            % in the progress bar because they run before the main pipeline.
            allKeys = PipelineProgressGUI.STEP_MAP(:, 1);
            allNames = PipelineProgressGUI.STEP_MAP(:, 2);

            obj.StepKeys = {};
            obj.StepDisplayNames = {};
            for k = 1:numel(stepsToRun)
                % Look up the step key in the master map to get its display name
                idx = find(strcmp(stepsToRun{k}, allKeys), 1);
                if ~isempty(idx)
                    obj.StepKeys{end+1} = allKeys{idx};
                    obj.StepDisplayNames{end+1} = allNames{idx};
                end
            end
            obj.TotalSteps = numel(obj.StepKeys);

            if obj.TotalSteps == 0
                % No recognized steps to track; skip GUI creation entirely
                obj.GUI = [];
                return;
            end

            % Create the underlying ProgressGUI figure with a DWI-type-specific title
            titleStr = sprintf('DWI Pipeline — %s', dwiTypeName);
            obj.GUI = ProgressGUI(titleStr, obj.TotalSteps);

            % Show the initial "0%" state before any step begins
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
