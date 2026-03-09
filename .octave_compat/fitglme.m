% FITGLME  Octave-compatible stub for MATLAB's fitglme (Statistics Toolbox).
%
%   mdl = fitglme(tbl, formula) fits a Generalized Linear Mixed-Effects
%   model in MATLAB (R2014b+). Octave does not have an equivalent function.
%
%   This stub returns a struct mimicking the key fields of a
%   GeneralizedLinearMixedModel object so that downstream code can access
%   mdl.Coefficients.Estimate, mdl.Coefficients.pValue, and mdl.Formula
%   without erroring.
%
%   Behavioral differences from MATLAB's fitglme:
%   - No actual model fitting is performed; estimates are zeros and p-values
%     are ones (maximally non-significant).
%   - Only two coefficients are returned (intercept + one predictor); the
%     real model may have more depending on the formula.
%   - Name-value pair arguments (Distribution, Link, etc.) are accepted
%     via varargin but silently ignored.
%   - This is a placeholder for test environments; real analyses require MATLAB.
function mdl = fitglme(tbl, formula, varargin)
    % Return a struct with placeholder coefficient estimates and p-values.
    mdl = struct();
    mdl.Coefficients = struct();
    mdl.Coefficients.Estimate = [0; 0];  % Intercept + one predictor, both zero.
    mdl.Coefficients.pValue = [1; 1];    % Non-significant placeholder p-values.
    mdl.Formula = formula;
end
