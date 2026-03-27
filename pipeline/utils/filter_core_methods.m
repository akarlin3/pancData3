function [active_methods, pruned_info] = filter_core_methods(all_methods, failure_table, config_struct)
%FILTER_CORE_METHODS Remove core methods exceeding failure rate cutoff.
%
%   Reads the failure rate table from compute_core_failure_rates and removes
%   methods whose total failure rate (any_failure_rate) exceeds the configured
%   threshold. Also removes methods listed in config.excluded_core_methods.
%
%   Inputs:
%       all_methods     - cell array of all method name strings (e.g. 11)
%       failure_table   - struct from compute_core_failure_rates containing
%                         .any_failure_rate [N x 3], .method_names,
%                         .median_core_voxels [N x 3]
%       config_struct   - pipeline config with fields:
%           .max_core_failure_rate   - scalar 0-1, prune methods above this
%           .excluded_core_methods   - cell array of method names to always exclude
%           .min_core_voxels         - minimum median core voxel count (0 = disabled)
%
%   Outputs:
%       active_methods  - cell array of method names that passed the filter
%       pruned_info     - struct array (one per pruned method) with fields:
%           .name, .reason, .failure_rate, .pipeline

    active_methods = all_methods(:)';
    pruned_info = struct('name', {}, 'reason', {}, 'failure_rate', {}, 'pipeline', {});

    ft_methods = failure_table.method_names;
    pipeline_names = failure_table.pipeline_names;

    % --- 1. Manual exclusions ---
    if isfield(config_struct, 'excluded_core_methods') && ~isempty(config_struct.excluded_core_methods)
        manual_list = config_struct.excluded_core_methods;
        for i = 1:numel(manual_list)
            mname = manual_list{i};
            % Safety: never prune adc_threshold
            if strcmp(mname, 'adc_threshold')
                warning('filter_core_methods:adcProtected', ...
                    'adc_threshold cannot be manually excluded (it is the fallback method).');
                continue;
            end
            if ismember(mname, active_methods)
                active_methods(strcmp(active_methods, mname)) = [];
                % Look up failure rate for logging
                ft_idx = find(strcmp(ft_methods, mname), 1);
                if ~isempty(ft_idx)
                    [max_fr, max_p] = max(failure_table.any_failure_rate(ft_idx, :));
                else
                    max_fr = NaN;
                    max_p = 1;
                end
                entry = struct('name', mname, 'reason', 'manual_exclusion', ...
                    'failure_rate', max_fr, 'pipeline', pipeline_names{max_p});
                pruned_info(end+1) = entry; %#ok<AGROW>
            end
        end
    end

    % --- 2. Failure rate filter ---
    max_rate = 1.0;
    if isfield(config_struct, 'max_core_failure_rate')
        max_rate = config_struct.max_core_failure_rate;
    end
    if max_rate < 1.0
        to_remove = {};
        for i = 1:numel(active_methods)
            mname = active_methods{i};
            ft_idx = find(strcmp(ft_methods, mname), 1);
            if isempty(ft_idx), continue; end

            [max_fr, max_p] = max(failure_table.any_failure_rate(ft_idx, :));
            if isnan(max_fr), continue; end

            if max_fr > max_rate
                % Safety: never prune adc_threshold
                if strcmp(mname, 'adc_threshold')
                    warning('filter_core_methods:adcProtected', ...
                        'adc_threshold exceeds failure threshold (%.1f%%) but is protected.', max_fr * 100);
                    continue;
                end
                to_remove{end+1} = mname; %#ok<AGROW>
                entry = struct('name', mname, 'reason', 'failure_rate', ...
                    'failure_rate', max_fr, 'pipeline', pipeline_names{max_p});
                pruned_info(end+1) = entry; %#ok<AGROW>
            end
        end
        active_methods = setdiff(active_methods, to_remove, 'stable');
    end

    % --- 3. Minimum voxel filter ---
    min_vox = 0;
    if isfield(config_struct, 'min_core_voxels')
        min_vox = config_struct.min_core_voxels;
    end
    if min_vox > 0 && isfield(failure_table, 'median_core_voxels')
        to_remove = {};
        for i = 1:numel(active_methods)
            mname = active_methods{i};
            ft_idx = find(strcmp(ft_methods, mname), 1);
            if isempty(ft_idx), continue; end

            med_vox = failure_table.median_core_voxels(ft_idx, :);
            % Use nanmedian across pipelines
            overall_med = nanmedian(med_vox(~isnan(med_vox)));
            if isnan(overall_med), continue; end

            if overall_med < min_vox
                if strcmp(mname, 'adc_threshold')
                    warning('filter_core_methods:adcProtected', ...
                        'adc_threshold has low median voxels (%.0f) but is protected.', overall_med);
                    continue;
                end
                to_remove{end+1} = mname; %#ok<AGROW>
                [max_fr, max_p] = max(failure_table.any_failure_rate(ft_idx, :));
                entry = struct('name', mname, 'reason', 'insufficient_voxels', ...
                    'failure_rate', max_fr, 'pipeline', pipeline_names{max_p});
                pruned_info(end+1) = entry; %#ok<AGROW>
            end
        end
        active_methods = setdiff(active_methods, to_remove, 'stable');
    end
end
