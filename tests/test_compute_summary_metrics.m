classdef test_compute_summary_metrics < matlab.unittest.TestCase
    % TEST_COMPUTE_SUMMARY_METRICS Unit tests for the compute_summary_metrics function.
    % Ensures that summary metrics are correctly calculated from raw vectors,
    % properly handling missing data and computing standard statistics.

    properties
        ConfigStruct
        DataVectors
        IDList
        MRNList
        LF
        Immuno
        GTVLoc
        DWILoc
        Dmean
        D95
        V50Gy
    end

    methods(TestMethodSetup)
        function createMockInputs(testCase)
            testCase.ConfigStruct = struct('dataloc', pwd, 'adc_thresh', 1.15e-3, ...
                'high_adc_thresh', 1e-3, 'd_thresh', 1e-3, 'f_thresh', 0.1, ...
                'dstar_thresh', 0.01, 'use_checkpoints', false, ...
                'min_vox_hist', 1, 'adc_max', 3e-3, 'dwi_types_to_run', 1);

            % Create mock data vectors for 2 patients, 3 timepoints.
            % Fields must match the names used in compute_summary_metrics.
            empty_entry = struct( ...
                'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
                'adc_vector_dncnn', [], 'd_vector_dncnn', [], 'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
                'd_vector_ivimnet', [], 'f_vector_ivimnet', [], 'dstar_vector_ivimnet', [], ...
                'vox_vol', 1);
            testCase.DataVectors = repmat(empty_entry, 2, 3);

            % Patient 1, Timepoint 1: Good data
            testCase.DataVectors(1,1).adc_vector = [0.001, 0.002, 0.0015]';
            testCase.DataVectors(1,1).d_vector   = [0.0008, 0.0012, 0.001]';
            testCase.DataVectors(1,1).f_vector   = [0.1, 0.2, 0.15]';
            testCase.DataVectors(1,1).dstar_vector = [0.01, 0.02, 0.015]';

            % Patient 1, Timepoint 2: Missing data (empty) — leave as []

            % Patient 1, Timepoint 3: Some NaN
            testCase.DataVectors(1,3).adc_vector = [0.001, NaN, 0.002]';
            testCase.DataVectors(1,3).d_vector   = [0.0008, NaN, 0.001]';
            testCase.DataVectors(1,3).f_vector   = [0.1, NaN, 0.2]';
            testCase.DataVectors(1,3).dstar_vector = [0.01, NaN, 0.02]';

            % Patient 2, Timepoint 1: f values below threshold for sub-volume test
            testCase.DataVectors(2,1).adc_vector = [0.001, 0.001, 0.001]';
            testCase.DataVectors(2,1).d_vector   = [0.001, 0.001, 0.001]';
            testCase.DataVectors(2,1).f_vector   = [0.05, 0.05, 0.05]';  % all < f_thresh(0.1)
            testCase.DataVectors(2,1).dstar_vector = [0.02, 0.02, 0.02]';

            % Other inputs
            testCase.IDList = {'PT1', 'PT2'};
            testCase.MRNList = {'1234', '5678'};
            testCase.LF = [0; 1];
            testCase.Immuno = [0; 0];
            testCase.GTVLoc = {'/path/gtv1', '/path/gtv2'};
            testCase.DWILoc = {'/path/dwi1', '/path/dwi2'};
            testCase.Dmean = [40, 50];
            testCase.D95 = [35, 45];
            testCase.V50Gy = [10, 20];
        end
    end

    methods(Test)
        function testBasicCalculation(testCase)
            % Ensure basic calculations run without error and output is correct structure
            summary = compute_summary_metrics(testCase.ConfigStruct, testCase.DataVectors, ...
                testCase.IDList, testCase.MRNList, testCase.LF, testCase.Immuno, ...
                testCase.GTVLoc, testCase.DWILoc, testCase.Dmean, testCase.D95, testCase.V50Gy);

            testCase.verifyTrue(isstruct(summary));
            testCase.verifyEqual(summary.id_list, testCase.IDList);
            testCase.verifyEqual(summary.mrn_list, testCase.MRNList);

            % Verify matrix sizes: 2 patients, 3 timepoints, 3 DWI types (always allocated as 3)
            sz = size(summary.adc_mean);
            testCase.verifyEqual(sz, [2, 3, 3]);

            % Verify exact mean values for PT1, TP1 (Standard pipeline = index 1)
            testCase.verifyEqual(summary.adc_mean(1,1,1), mean([0.001, 0.002, 0.0015]), 'AbsTol', 1e-10);
            testCase.verifyEqual(summary.f_mean(2,1,1), 0.05, 'AbsTol', 1e-10);

            % Verify NaN handling for PT1, TP3 (mean of non-NaN: [0.001, 0.002])
            testCase.verifyEqual(summary.adc_mean(1,3,1), mean([0.001, 0.002]), 'AbsTol', 1e-10);

            % Verify empty handling for PT1, TP2
            testCase.verifyTrue(isnan(summary.adc_mean(1,2,1)));
        end

        function testSubvolumeCalculation(testCase)
            summary = compute_summary_metrics(testCase.ConfigStruct, testCase.DataVectors, ...
                testCase.IDList, testCase.MRNList, testCase.LF, testCase.Immuno, ...
                testCase.GTVLoc, testCase.DWILoc, testCase.Dmean, testCase.D95, testCase.V50Gy);

            % PT1, TP1 ADC: [0.001, 0.002, 0.0015], Thresh: 0.00115.
            % Sub-voxels: only [0.001] qualifies (1 voxel). vox_vol=1 => adc_sub_vol=1.
            % adc_sub_vol_pc = 1/3.
            testCase.verifyEqual(summary.adc_sub_vol_pc(1,1,1), 1/3, 'AbsTol', 1e-5);

            % PT2, TP1 f: [0.05, 0.05, 0.05], Thresh: 0.1.
            % All 3 voxels are below threshold. vox_vol=1 => f_sub_vol=3.
            testCase.verifyEqual(summary.f_sub_vol(2,1,1), 3.0, 'AbsTol', 1e-5);
        end
    end
end
