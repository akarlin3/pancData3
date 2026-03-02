classdef test_vis_script < matlab.unittest.TestCase
    % Smoke test for visualize_results

    methods (Test)
        function testVisualizeResultsRuns(testCase)
            patID = 'P01';
            dwi_img = zeros(10, 10, 10, 4, 'double');
            dwi_img(:,:,:,1) = 1000;
            dwi_img(:,:,:,2) = 800;
            dwi_img(:,:,:,3) = 500;
            dwi_img(:,:,:,4) = 100;

            outDir = tempname;
            mkdir(outDir);
            cleanup = onCleanup(@() rmdir(outDir, 's'));

            ConfigStruct.dataloc = pwd;
            ConfigStruct.output_folder = outDir;
            ConfigStruct.dwi_types_to_run = 1;

            SummaryMetrics.id_list = {patID};
            SummaryMetrics.mrn_list = {'MRN01'};
            SummaryMetrics.lf = [0];
            SummaryMetrics.adc_mean = 1.0e-3 * ones(1,1,1);
            SummaryMetrics.d_mean = 1.0e-3 * ones(1,1,1);
            SummaryMetrics.f_mean = 0.1 * ones(1,1,1);
            SummaryMetrics.dstar_mean = 0.05 * ones(1,1,1);
            SummaryMetrics.d95_gtvp = 40 * ones(1,1);
            SummaryMetrics.dmean_gtvp = 50 * ones(1,1);

            CalculatedResults = struct();
            DataVectors = struct('adc_vector', {ones(10,1)});

            visualize_results(DataVectors, SummaryMetrics, CalculatedResults, ConfigStruct);
            testCase.verifyTrue(true, 'visualize_results completed without error');
        end
    end
end
