classdef run_manual_test < matlab.unittest.TestCase
% RUN_MANUAL_TEST  Diagnostic smoke test for compute_summary_metrics.
%   Constructs a minimal inline mock dataset (2 patients, 3 timepoints,
%   1 DWI type) with realistic IVIM parameter values and verifies that
%   compute_summary_metrics runs to completion and returns a struct.
%
%   This test is useful for quick manual verification after changes to
%   compute_summary_metrics or its sub-functions (compute_adc_metrics,
%   compute_ivim_metrics).

    methods (Test)
        function testComputeSummaryMetricsRuns(testCase)
        %TESTCOMPUTESUMMARYMETRICSRUNS Build mock data and verify that
        %   compute_summary_metrics produces a valid struct output.

            % Minimal config with threshold values matching typical clinical ranges
            ConfigStruct = struct('dataloc', pwd, 'adc_thresh', 1e-3, ...
                'high_adc_thresh', 1.15e-3, 'd_thresh', 1e-3, 'f_thresh', 0.1, ...
                'dstar_thresh', 0.01, 'use_checkpoints', false, 'dwi_types_to_run', [1], ...
                'min_vox_hist', 1, 'adc_max', 3e-3, 'core_method', 'adc_threshold');

            % Initialize DataVectors: 2 patients x 3 timepoints x 1 DWI type
            % All fields start empty; specific entries are populated below.
            DataVectors = repmat(struct('adc_vector', [], 'd_vector', [], ...
                'f_vector', [], 'dstar_vector', [], 'vox_vol', 1), 2, 3, 1);

            % Patient 1, Timepoint 1: Good data (3 voxels, no NaN)
            DataVectors(1,1,1).adc_vector = [0.001, 0.002, 0.0015];
            DataVectors(1,1,1).d_vector = [0.0008, 0.0012, 0.001];
            DataVectors(1,1,1).f_vector = [0.1, 0.2, 0.15];
            DataVectors(1,1,1).dstar_vector = [0.01, 0.02, 0.015];

            % Patient 1, Timepoint 3: Partially missing data (NaN in 2nd voxel)
            DataVectors(1,3,1).adc_vector = [0.001, NaN, 0.002];
            DataVectors(1,3,1).d_vector = [0.0008, NaN, 0.001];
            DataVectors(1,3,1).f_vector = [0.1, NaN, 0.2];
            DataVectors(1,3,1).dstar_vector = [0.01, NaN, 0.02];

            % Patient 2, Timepoint 1: Uniform data (all voxels identical)
            DataVectors(2,1,1).adc_vector = [0.001, 0.001, 0.001];
            DataVectors(2,1,1).d_vector = [0.001, 0.001, 0.001];
            DataVectors(2,1,1).f_vector = [0.2, 0.2, 0.2];
            DataVectors(2,1,1).dstar_vector = [0.02, 0.02, 0.02];

            % Clinical metadata for the 2 patients
            IDList = {'PT1', 'PT2'};
            MRNList = {'1234', '5678'};
            LF = [0; 1];           % Local failure status (0=control, 1=failure)
            Immuno = [0; 0];       % Immunotherapy flag
            GTVLoc = {'/path/gtv1', '/path/gtv2'};
            DWILoc = {'/path/dwi1', '/path/dwi2'};
            Dmean = [40, 50];      % Mean dose (Gy)
            D95 = [35, 45];        % D95 dose (Gy)
            V50Gy = [10, 20];      % Volume receiving >= 50 Gy (cc)

            % Run the function under test and verify output type
            summary = compute_summary_metrics(ConfigStruct, DataVectors, ...
                IDList, MRNList, LF, Immuno, GTVLoc, DWILoc, Dmean, D95, V50Gy);
            testCase.verifyTrue(isstruct(summary), 'Expected struct output');
        end
    end
end
