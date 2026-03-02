pkg load statistics;
addpath('../../core');

% Emulate test_compute_summary_metrics
ConfigStruct = struct('dataloc', pwd, 'adc_thresh', 1.15e-3, ...
    'high_adc_thresh', 1e-3, 'd_thresh', 1e-3, 'f_thresh', 0.1, ...
    'dstar_thresh', 0.01, 'use_checkpoints', false, 'dwi_types_to_run', [1], 'min_vox_hist', 1, 'adc_max', 3e-3);

% Create mock data vectors for 2 patients, 3 timepoints
DataVectors = repmat(struct('adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], 'vox_vol', 1), 2, 3, 1);

% Patient 1, Timepoint 1: Good data
DataVectors(1,1,1).adc_vector = [0.001, 0.002, 0.0015];
DataVectors(1,1,1).d_vector = [0.0008, 0.0012, 0.001];
DataVectors(1,1,1).f_vector = [0.1, 0.2, 0.15];
DataVectors(1,1,1).dstar_vector = [0.01, 0.02, 0.015];

% Patient 1, Timepoint 2: Missing data (empty)
% Patient 1, Timepoint 3: Some NaN
DataVectors(1,3,1).adc_vector = [0.001, NaN, 0.002];
DataVectors(1,3,1).d_vector = [0.0008, NaN, 0.001];
DataVectors(1,3,1).f_vector = [0.1, NaN, 0.2];
DataVectors(1,3,1).dstar_vector = [0.01, NaN, 0.02];

% Patient 2, Timepoint 1: Good data
DataVectors(2,1,1).adc_vector = [0.001, 0.001, 0.001];
DataVectors(2,1,1).d_vector = [0.001, 0.001, 0.001];
DataVectors(2,1,1).f_vector = [0.2, 0.2, 0.2];
DataVectors(2,1,1).dstar_vector = [0.02, 0.02, 0.02];

IDList = {'PT1', 'PT2'};
MRNList = {'1234', '5678'};
LF = [0; 1];
Immuno = [0; 0];
GTVLoc = {'/path/gtv1', '/path/gtv2'};
DWILoc = {'/path/dwi1', '/path/dwi2'};
Dmean = [40, 50];
D95 = [35, 45];
V50Gy = [10, 20];

try
    summary = compute_summary_metrics(ConfigStruct, DataVectors, IDList, MRNList, LF, Immuno, GTVLoc, DWILoc, Dmean, D95, V50Gy);
    disp('Computed successfully');
catch e
    disp('Error:');
    disp(e.message);
    disp(e.stack(1));
end
