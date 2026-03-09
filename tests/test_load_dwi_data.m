classdef test_load_dwi_data < matlab.unittest.TestCase
    % TEST_LOAD_DWI_DATA Unit tests for data loading and reloading logic
    %
    % Validates the skip_to_reload branching, fallback mechanisms for missing
    % dwi_type_name variants, file-not-found error handling, and empty directory
    % discovery behavior in load_dwi_data.m.

    properties
        TempDataLoc      % Temporary directory acting as the pipeline's dataloc
        ConfigStruct     % Configuration struct with skip_to_reload=true by default
        OriginalPath     % Saved MATLAB path for restoration in teardown
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Creates a temporary directory and a minimal config struct
            % configured for reload-mode testing. The config uses
            % skip_to_reload=true so that load_dwi_data attempts to read
            % a pre-existing .mat file rather than running the full DICOM
            % conversion and model fitting pipeline.
            testCase.TempDataLoc = tempname;
            mkdir(testCase.TempDataLoc);

            % Clinical data sheet path — not actually written here; tests
            % that need it create their own via writetable.
            dummy_sheet = fullfile(testCase.TempDataLoc, 'mock.xlsx');

            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.dataloc = testCase.TempDataLoc;
            testCase.ConfigStruct.skip_to_reload = true;
            testCase.ConfigStruct.ivim_bthr = 100;
            testCase.ConfigStruct.dcm2nii_call = 'dummy';
            testCase.ConfigStruct.clinical_data_sheet = 'mock.xlsx';
            testCase.ConfigStruct.adc_thresh = 0.001;
            testCase.ConfigStruct.high_adc_thresh = 0.00115;
            testCase.ConfigStruct.d_thresh = 0.001;
            testCase.ConfigStruct.f_thresh = 0.1;
            testCase.ConfigStruct.adc_max = 0.003;
            testCase.ConfigStruct.min_vox_hist = 100;
            testCase.ConfigStruct.dwi_types_to_run = 1:3;
            testCase.ConfigStruct.use_checkpoints = false;
            testCase.ConfigStruct.core_method = 'adc_threshold';

            % Add core and utils to path so we can call load_dwi_data
            testCase.OriginalPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods(TestMethodTeardown)
        function teardownEnvironment(testCase)
            % Remove all temporary files and restore the original MATLAB
            % path to prevent test pollution between methods.
            if exist(testCase.TempDataLoc, 'dir')
                rmdir(testCase.TempDataLoc, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods
        function createDummySave(testCase, filename)
            % Creates a minimal but valid dwi_vectors .mat file containing
            % all variables that load_dwi_data expects when reloading:
            % data_vectors_gtvp/gtvn (struct arrays), patient metadata
            % (id_list, mrn_list, lf, immuno), DICOM dates, file locations,
            % and dose metrics. These mimic the output of Section 3 of
            % load_dwi_data.m (the full DICOM processing path).
            data_vectors_gtvp = repmat(struct('adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], 'adc_vector_dncnn', [], 'd_vector_dncnn', [], 'f_vector_dncnn', [], 'dstar_vector_dncnn', [], 'd_vector_ivimnet', [], 'f_vector_ivimnet', [], 'dstar_vector_ivimnet', [], 'ID', 'P01', 'MRN', 'MRN01', ...
                'LF', 0, 'Immuno', 0, 'Fraction', 1, 'Repeatability_index', 1, 'vox_vol', 1), [1, 6, 6]);
            data_vectors_gtvn = repmat(struct('adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], 'adc_vector_dncnn', [], 'd_vector_dncnn', [], 'f_vector_dncnn', [], 'dstar_vector_dncnn', [], 'd_vector_ivimnet', [], 'f_vector_ivimnet', [], 'dstar_vector_ivimnet', [], 'ID', 'P01', 'MRN', 'MRN01'), [1, 6, 6]);

            lf = 0;
            immuno = 0;
            mrn_list = {'MRN01'};
            id_list = {'P01'};
            fx_dates = {'20230101'};
            dwi_locations = cell(1,1,1);
            rtdose_locations = cell(1,1);
            gtv_locations = cell(1,1,1);
            gtvn_locations = cell(1,1,1);
            dmean_gtvp = 50;
            dmean_gtvn = nan;
            d95_gtvp = 45;
            d95_gtvn = nan;
            v50gy_gtvp = 90;
            v50gy_gtvn = nan;
            bad_dwi_locations = {};
            bad_dwi_count = 0;

            save(fullfile(testCase.TempDataLoc, filename), ...
                'data_vectors_gtvp', 'data_vectors_gtvn', 'lf', 'immuno', ...
                'mrn_list', 'id_list', 'fx_dates', 'dwi_locations', ...
                'rtdose_locations', 'gtv_locations', 'gtvn_locations', ...
                'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', ...
                'v50gy_gtvp', 'v50gy_gtvn', 'bad_dwi_locations', 'bad_dwi_count');
        end
    end

    methods(Test)

        function testSkipToReloadSuccess(testCase)
            % Verifies the happy path for skip_to_reload=true: when
            % dwi_vectors.mat exists in dataloc, load_dwi_data should
            % successfully load it, run compute_summary_metrics, and return
            % non-empty GTVp data vectors and a summary_metrics struct
            % with the expected patient ID.

            % 1. Create a dummy dwi_vectors.mat
            testCase.createDummySave('dwi_vectors.mat');

            % 2. Execute load_dwi_data — should succeed and return valid outputs
            try
                [gtvp, gtvn, summary] = load_dwi_data(testCase.ConfigStruct);

                % 3. Verify outputs
                testCase.verifyNotEmpty(gtvp, 'GTVp data vectors should not be empty');
                testCase.verifyNotEmpty(summary, 'Summary metrics should not be empty');
                testCase.verifyEqual(summary.id_list{1}, 'P01', 'ID list should match mock data');
            catch ME
                testCase.verifyFail(sprintf('Failed with error: %s', ME.message));
            end
        end

        function testSkipToReloadFallback(testCase)
            % Verifies that if a specific dwi_type_name is requested but missing,
            % the code falls back to loading the default 'dwi_vectors.mat'

            % 1. Request a specific DWI pipeline variant
            testCase.ConfigStruct.dwi_type_name = 'DnCNN';

            % 2. Create ONLY the default dwi_vectors.mat (simulate missing specific variant)
            testCase.createDummySave('dwi_vectors.mat');

            % 3. Execute - it should fallback and succeed
            try
                [gtvp, ~, summary] = load_dwi_data(testCase.ConfigStruct);
                testCase.verifyEqual(summary.id_list{1}, 'P01', 'Should fallback to default file');
            catch ME
                testCase.verifyFail(sprintf('Fallback failed with error: %s', ME.message));
            end
        end

        function testSkipToReloadMissingError(testCase)
            % Verifies that if no data file exists and skip_to_reload is true,
            % an appropriate error is thrown

            % 1. Ensure directory is empty
            % 2. Execute and expect error
            testCase.verifyError(@() load_dwi_data(testCase.ConfigStruct), ...
                ?MException, 'Should throw an error when dwi_vectors.mat is missing');
        end

        function testSkipToReloadMissingSpecificFallthrough(testCase)
            % Verifies that if a specific DWI type variant is requested with
            % skip_to_reload=true but the .mat file doesn't exist, the code
            % falls through to the full load path (skip_to_reload overridden
            % to false).  The full load will fail downstream for other
            % reasons (no clinical sheet, etc.), but NOT with the
            % load_dwi_data:fileNotFound error.

            testCase.ConfigStruct.dwi_type_name = 'IVIMnet';

            try
                load_dwi_data(testCase.ConfigStruct);
                testCase.verifyFail('Expected an error from the full load path');
            catch ME
                testCase.verifyNotEqual(ME.identifier, 'load_dwi_data:fileNotFound', ...
                    'Should NOT get fileNotFound error — should fall through to full load');
            end
        end

        function testDiscoverFilesEmpty(testCase)
            % Verifies behavior when skip_to_reload is false but the data
            % directory contains no valid patient DICOM folders. The function
            % should complete without error: Section 2 reads the clinical
            % sheet, Section 3's parfor loop finds no patients to process,
            % and the resulting empty data structures are saved and returned.

            testCase.ConfigStruct.skip_to_reload = false;

            % Create a minimal clinical Excel file so readtable does not crash
            T = table({'P01'}, 0, 0, 'VariableNames', {'Pat', 'LF', 'Immuno'});
            sheet_path = fullfile(testCase.TempDataLoc, testCase.ConfigStruct.clinical_data_sheet);
            writetable(T, sheet_path, 'Sheet', 'Clin List_MR');

            try
                [gtvp, gtvn, summary] = load_dwi_data(testCase.ConfigStruct);

                % Since there are no DICOM folders, output structs should be empty
                % or initialized to empty states.
                % Because skip_to_reload is false, Section 3 will save a new dwi_vectors.mat
                % and Section 4 will reload it, and compute_summary_metrics will process empty arrays.

                testCase.verifyTrue(isempty(gtvp) || numel(fieldnames(gtvp)) == 0, ...
                    'GTVp should be empty struct');
                testCase.verifyTrue(isempty(gtvn) || numel(fieldnames(gtvn)) == 0, ...
                    'GTVn should be empty struct');
                testCase.verifyEmpty(summary.id_list, 'Summary id_list should be empty');

                % Verify it saved the empty state
                testCase.verifyTrue(exist(fullfile(testCase.TempDataLoc, 'dwi_vectors.mat'), 'file') == 2, ...
                    'Should save empty state to dwi_vectors.mat');
            catch ME
                testCase.verifyFail(sprintf('Empty discovery failed: %s', ME.message));
            end
        end
    end
end
