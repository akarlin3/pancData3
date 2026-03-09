classdef test_vis_script < matlab.unittest.TestCase
% TEST_VIS_SCRIPT  Diagnostic smoke test for visualize_results.
%   Constructs a minimal mock dataset (1 patient, 1 DWI type) with
%   synthetic summary metrics and verifies that visualize_results runs
%   to completion without error. Does not validate figure content —
%   just confirms no crashes occur with the minimal input set.

    methods (Test)
        function testVisualizeResultsRuns(testCase)
        %TESTVISUALIZERESULTSRUNS Run visualize_results with mock data
        %   and verify it completes without throwing an error.

            patID = 'P01';

            % Create a synthetic 4D DWI image (10x10x10 volume, 4 b-values)
            % with decreasing signal to mimic diffusion attenuation
            dwi_img = zeros(10, 10, 10, 4, 'double');
            dwi_img(:,:,:,1) = 1000;  % b=0
            dwi_img(:,:,:,2) = 800;   % low b
            dwi_img(:,:,:,3) = 500;   % medium b
            dwi_img(:,:,:,4) = 100;   % high b

            % Create a temporary output directory (auto-cleaned on scope exit)
            outDir = tempname;
            mkdir(outDir);
            cleanup = onCleanup(@() rmdir(outDir, 's'));

            % Minimal config for visualization
            ConfigStruct.dataloc = pwd;
            ConfigStruct.output_folder = outDir;
            ConfigStruct.dwi_types_to_run = 1;

            % Minimal SummaryMetrics struct with one patient's worth of data
            SummaryMetrics.id_list = {patID};
            SummaryMetrics.mrn_list = {'MRN01'};
            SummaryMetrics.lf = [0];
            SummaryMetrics.adc_mean = 1.0e-3 * ones(1,1,1);
            SummaryMetrics.d_mean = 1.0e-3 * ones(1,1,1);
            SummaryMetrics.f_mean = 0.1 * ones(1,1,1);
            SummaryMetrics.dstar_mean = 0.05 * ones(1,1,1);
            SummaryMetrics.d95_gtvp = 40 * ones(1,1);
            SummaryMetrics.dmean_gtvp = 50 * ones(1,1);

            % Empty structs for fields not exercised by this smoke test
            CalculatedResults = struct();
            DataVectors = struct('adc_vector', {ones(10,1)});

            % Call the function under test
            visualize_results(DataVectors, SummaryMetrics, CalculatedResults, ConfigStruct);
            testCase.verifyTrue(true, 'visualize_results completed without error');
        end
    end
end
