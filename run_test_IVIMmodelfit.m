% run_test_IVIMmodelfit.m
% 
% Summary:
%   A helper script to quickly run tests specifically targeting the 
%   IVIMmodelfit function and display the results in a tabular format.
%
% Usage:
%   Simply run this script to execute the tests and disp results.

results = runtests('tests/test_IVIMmodelfit.m');
disp(table(results));
