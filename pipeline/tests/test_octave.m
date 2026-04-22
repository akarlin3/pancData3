% TEST_OCTAVE Smoke test for the load_dwi_data reload path under Octave-like conditions.
%
% This script exercises the "skip_to_reload" code path in load_dwi_data.m,
% which reloads previously saved DWI vectors from a .mat checkpoint file
% instead of re-processing DICOM data from scratch. This is the fast-restart
% path used when checkpointed data already exists on disk.
%
% The test:
%   1. Constructs a minimal config_struct with skip_to_reload=true
%   2. Creates a dummy pipeline_voxels_Standard.mat with all required variables
%      (empty arrays for each field that load_dwi_data expects)
%   3. Calls load_dwi_data() and verifies it completes without error
%
% Note: Despite the filename "test_octave", this is not an Octave-specific
% test. It tests general reload functionality that must work in both MATLAB
% and Octave environments.

% Add core/ and utils/ to the MATLAB path
baseDir = fileparts(fileparts(which(mfilename)));
addpath(fullfile(baseDir, 'core'), fullfile(baseDir, 'utils'));

% Build a minimal config struct with skip_to_reload enabled.
% dataloc points to the system temp directory where we will place
% the checkpoint .mat file.
config_struct.dataloc = tempdir;
config_struct.dwi_type_name = 'Standard';
config_struct.skip_to_reload = true;       % Skip DICOM processing, load from checkpoint
config_struct.ivim_bthr = 100;             % IVIM b-value threshold (s/mm^2)
config_struct.adc_thresh = 0.001;          % ADC threshold for tumor core (mm^2/s)
config_struct.high_adc_thresh = 0.00115;   % High ADC threshold (mm^2/s)
config_struct.d_thresh = 0.001;            % D threshold for IVIM core (mm^2/s)
config_struct.f_thresh = 0.1;              % f threshold for IVIM core (fraction)
config_struct.adc_max = 0.003;             % Maximum plausible ADC (mm^2/s)
config_struct.min_vox_hist = 5;            % Minimum voxels for histogram analysis
config_struct.use_gpu = false;             % GPU acceleration (disabled for test)
config_struct.gpu_device = 1;              % GPU device index

% Single test patient
id_list = {'test1'};

% Initialize all variables that load_dwi_data expects to find in the
% checkpoint .mat file. All are empty since this is a minimal smoke test.
data_vectors_gtvn=[]; data_vectors_gtvp=[]; lf=[]; immuno=[]; mrn_list={}; fx_dates=[]; dwi_locations=[]; rtdose_locations=[]; gtv_locations=[]; gtvn_locations=[]; dmean_gtvp=[]; dmean_gtvn=[]; d95_gtvp=[]; d95_gtvn=[]; v50gy_gtvp=[]; v50gy_gtvn=[]; bad_dwi_locations=[]; bad_dwi_count=[];

% Save the checkpoint file that load_dwi_data will attempt to reload
save(fullfile(config_struct.dataloc, 'pipeline_voxels_Standard.mat'), 'data_vectors_gtvn', 'data_vectors_gtvp', 'lf', 'immuno', 'mrn_list', 'id_list', 'fx_dates', 'dwi_locations', 'rtdose_locations', 'gtv_locations', 'gtvn_locations', 'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', 'v50gy_gtvp', 'v50gy_gtvn', 'bad_dwi_locations', 'bad_dwi_count');

% Attempt the reload path; success means load_dwi_data can parse the
% checkpoint file and return without error.
try
    [dvp, dvn, summary_metrics] = load_dwi_data(config_struct);
    fprintf('SUCCESS\n');
catch e
    fprintf('ERROR: %s\n', e.message);
end
