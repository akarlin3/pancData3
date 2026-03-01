classdef test_load_dwi_data < matlab.unittest.TestCase
    % TEST_LOAD_DWI_DATA Unit tests for data loading and reloading logic
    %
    % Validates the skip_to_reload branching, fallback mechanisms for missing
    % dwi_type_name variants, file-not-found error handling, and empty directory
    % discovery behavior in load_dwi_data.m.

    properties
        TempDataLoc
        ConfigStruct
        OriginalPath
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Create a temporary directory to act as dataloc
            testCase.TempDataLoc = tempname;
            mkdir(testCase.TempDataLoc);

            % Create a dummy clinical data sheet as it's loaded in Section 2
            % readtable('mock.xlsx', 'Sheet', 'Clin List')
            % Use .xlsx because readtable/writetable use the 'Sheet' argument
            dummy_sheet = fullfile(testCase.TempDataLoc, 'mock.xlsx');

            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.dataloc = testCase.TempDataLoc;
            testCase.ConfigStruct.skip_to_reload = true;
            testCase.ConfigStruct.ivim_bthr = 100;
            testCase.ConfigStruct.dcm2nii_call = 'dummy';
            testCase.ConfigStruct.clinical_data_sheet = 'mock.xlsx';

            % Add core and utils to path so we can call load_dwi_data
            testCase.OriginalPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardownEnvironment(testCase)
            % Clean up temporary directory
            if exist(testCase.TempDataLoc, 'dir')
                rmdir(testCase.TempDataLoc, 's');
            end
            % Restore path
            path(testCase.OriginalPath);
        end
    end

    methods
        function createDummySave(testCase, filename)
            % Creates a valid dwi_vectors.mat file with all required variables
            % needed by compute_summary_metrics

            % We need to mock the variables saved at the end of Section 3
            data_vectors_gtvp = struct('adc_vector', [], 'ID', 'P01', 'MRN', 'MRN01', ...
                'LF', 0, 'Immuno', 0, 'Fraction', 1, 'Repeatability_index', 1, 'vox_vol', 1);
            data_vectors_gtvn = struct('adc_vector', [], 'ID', 'P01', 'MRN', 'MRN01');

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
            % Verifies that skip_to_reload = true successfully loads the default
            % dwi_vectors.mat and processes it through compute_summary_metrics

            % 1. Create a dummy dwi_vectors.mat
            testCase.createDummySave('dwi_vectors.mat');

            % 2. Execute load_dwi_data
            % Should not throw an error and should return summary_metrics
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

        function testSkipToReloadMissingSpecificError(testCase)
            % Verifies that if specific variant is requested, AND fallback is missing,
            % it errors out completely.

            testCase.ConfigStruct.dwi_type_name = 'IVIMnet';

            testCase.verifyError(@() load_dwi_data(testCase.ConfigStruct), ...
                ?MException, 'Should throw error when both specific and fallback are missing');
        end

        function testDiscoverFilesEmpty(testCase)
            % Verifies behavior when skip_to_reload is false but the data directory
            % contains no valid patient folders.

            testCase.ConfigStruct.skip_to_reload = false;

            % Create an empty clinical sheet to prevent readtable from crashing
            % We need a valid dummy Excel file or a mocked table if it errors
            % But if id_list is empty, the parfor loop won't run, but readtable
            % happens BEFORE the loop in Section 2.

            % Mocking the clinical sheet creation for readtable
            T = table({'P01'}, 0, 0, 'VariableNames', {'Pat', 'LF', 'Immuno'});
            sheet_path = fullfile(testCase.TempDataLoc, testCase.ConfigStruct.clinical_data_sheet);
            writetable(T, sheet_path, 'Sheet', 'Clin List');

            try
                [gtvp, gtvn, summary] = load_dwi_data(testCase.ConfigStruct);

                % Since there are no DICOM folders, output structs should be empty
                % or initialized to empty states.
                % Because skip_to_reload is false, Section 3 will save a new dwi_vectors.mat
                % and Section 4 will reload it, and compute_summary_metrics will process empty arrays.

                testCase.verifyEmpty(gtvp, 'GTVp should be empty struct');
                testCase.verifyEmpty(gtvn, 'GTVn should be empty struct');
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
