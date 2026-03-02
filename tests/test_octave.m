baseDir = fileparts(fileparts(which(mfilename)));
addpath(fullfile(baseDir, 'core'), fullfile(baseDir, 'utils'));
config_struct.dataloc = tempdir;
config_struct.dwi_type_name = 'Standard';
config_struct.skip_to_reload = true;
config_struct.ivim_bthr = 100;
config_struct.adc_thresh = 0.00115;
config_struct.high_adc_thresh = 0.0015;
config_struct.d_thresh = 0.00115;
config_struct.f_thresh = 0.2;
config_struct.adc_max = 0.003;
config_struct.min_vox_hist = 5;

id_list = {'test1'};

data_vectors_gtvn=[]; data_vectors_gtvp=[]; lf=[]; immuno=[]; mrn_list={}; fx_dates=[]; dwi_locations=[]; rtdose_locations=[]; gtv_locations=[]; gtvn_locations=[]; dmean_gtvp=[]; dmean_gtvn=[]; d95_gtvp=[]; d95_gtvn=[]; v50gy_gtvp=[]; v50gy_gtvn=[]; bad_dwi_locations=[]; bad_dwi_count=[];

save(fullfile(config_struct.dataloc, 'dwi_vectors_Standard.mat'), 'data_vectors_gtvn', 'data_vectors_gtvp', 'lf', 'immuno', 'mrn_list', 'id_list', 'fx_dates', 'dwi_locations', 'rtdose_locations', 'gtv_locations', 'gtvn_locations', 'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', 'v50gy_gtvp', 'v50gy_gtvn', 'bad_dwi_locations', 'bad_dwi_count');

try
    [dvp, dvn, summary_metrics] = load_dwi_data(config_struct);
    fprintf('SUCCESS\n');
catch e
    fprintf('ERROR: %s\n', e.message);
end
