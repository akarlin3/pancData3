function [adc_thresh, info] = resolve_adc_thresh(source, preset_thresh, opt_results)
%RESOLVE_ADC_THRESH  Pick the ADC threshold that defines the sub-volume.
%
%   [adc_thresh, info] = resolve_adc_thresh(source, preset_thresh, opt_results)
%
%   Maps the config field `adc_thresh_source` to a concrete ADC threshold for
%   sub-volume creation, so the cut can be either pre-specified OR derived
%   from the threshold optimizer's tactics:
%
%     'preset'  -> the pre-specified config.adc_thresh (default; current
%                  behaviour, and the clean escape from selection bias).
%     'tactic1' -> reproducibility optimum (max Fx1-repeat Dice).
%     'tactic2' -> volume-inflection knee.
%     'tactic3' -> outcome-significance cut, BIAS-CORRECTED.  Tactic 3's raw
%                  in-sample threshold (opt_results.significance_thresh) is
%                  optimistically biased (selection-on-outcome; see
%                  docs/THRESHOLD_VALIDATION.md), so it MUST NOT feed
%                  sub-volume creation.  Instead the nested-CV recommended
%                  (out-of-fold) threshold is used -- repeated grouped K-fold
%                  if available, else leave-one-out.
%
%   ANY requested tactic that is unavailable or NaN (e.g. the optimizer did
%   not run, or the cohort lacked outcome labels) falls back to the
%   pre-specified preset with a logged note -- the function NEVER silently
%   substitutes the biased in-sample Tactic-3 value.
%
%   Inputs:
%       source       - char/string: 'preset'|'tactic1'|'tactic2'|'tactic3'
%       preset_thresh- scalar pre-specified threshold (config.adc_thresh)
%       opt_results  - struct from optimize_adc_threshold (optional; required
%                      only for tactic sources)
%
%   Outputs:
%       adc_thresh   - the resolved threshold (mm^2/s)
%       info         - struct: requested_source, source_used, value,
%                      fell_back (logical), note (human-readable)

    if nargin < 3, opt_results = struct(); end
    if isempty(source), source = 'preset'; end
    source = lower(char(source));

    info = struct();
    info.requested_source = source;
    info.source_used      = 'preset';
    info.value            = preset_thresh;
    info.fell_back        = false;
    info.note             = '';

    switch source
        case 'preset'
            adc_thresh = preset_thresh;
            info.note = sprintf('pre-specified threshold %.4g mm^2/s', preset_thresh);
            return;

        case 'tactic1'
            val = local_field(opt_results, 'optimal_thresh');
            label = 'Tactic 1 (reproducibility / max Fx1 Dice)';

        case 'tactic2'
            val = local_field(opt_results, 'inflection_thresh');
            label = 'Tactic 2 (volume inflection)';

        case 'tactic3'
            % Bias-corrected ONLY: nested-CV out-of-fold recommendation.
            % Never the raw in-sample significance_thresh.
            val = NaN; cv_label = '';
            if isfield(opt_results, 'nested_cv_repeated')
                val = local_field(opt_results.nested_cv_repeated, 'recommended_thresh');
                cv_label = 'repeated grouped K-fold';
            end
            if isnan(val) && isfield(opt_results, 'nested_cv')
                val = local_field(opt_results.nested_cv, 'recommended_thresh');
                cv_label = 'leave-one-out';
            end
            label = sprintf('Tactic 3 (outcome significance, nested-CV %s)', cv_label);

        otherwise
            adc_thresh = preset_thresh;
            info.fell_back = true;
            info.note = sprintf(['unrecognized adc_thresh_source "%s"; using ' ...
                'pre-specified %.4g mm^2/s'], source, preset_thresh);
            return;
    end

    if isnan(val)
        adc_thresh = preset_thresh;
        info.fell_back = true;
        info.note = sprintf(['%s threshold unavailable (optimizer not run or ' ...
            'inestimable); using pre-specified %.4g mm^2/s'], label, preset_thresh);
        return;
    end

    adc_thresh = val;
    info.source_used = source;
    info.value = val;
    info.note = sprintf('%s threshold %.4g mm^2/s', label, val);
end


function v = local_field(s, name)
% Return s.(name) as a finite scalar, or NaN if absent/empty/non-finite.
    v = NaN;
    if isstruct(s) && isfield(s, name)
        x = s.(name);
        if isscalar(x) && isnumeric(x) && isfinite(x)
            v = x;
        end
    end
end
