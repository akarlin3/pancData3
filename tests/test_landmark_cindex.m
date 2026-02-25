% test_landmark_cindex.m
clear; clc;
try
    % Let's load the essential variables for metrics.m
    % load_dwi_data_code pipeline usually outputs these
    load('dl_validation_manifest.mat', 'dl_provenance');
    
    % We will just run the script up to the new piece by running metrics.m
    % but maybe metrics.m takes too long to run entirely.
    % We will use evaluate_matlab_code to source metrics.m
catch
    disp('Error setting up the test.');
end
