patID = 'P01';
dwi_img = zeros(10, 10, 10, 4, 'double');
dwi_img(:,:,:,1) = 1000; % b=0
dwi_img(:,:,:,2) = 800;  % b=30
dwi_img(:,:,:,3) = 500;  % b=150
dwi_img(:,:,:,4) = 100;  % b=550

ConfigStruct.dataloc = pwd;
ConfigStruct.output_folder = fullfile(pwd, 'saved_figures_test');
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

addpath(fullfile(pwd, 'core'));
addpath(fullfile(pwd, 'utils'));
addpath(fullfile(pwd, 'dependencies'));

visualize_results(DataVectors, SummaryMetrics, CalculatedResults, ConfigStruct);
