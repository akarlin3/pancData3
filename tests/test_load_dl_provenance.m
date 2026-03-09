classdef test_load_dl_provenance < matlab.unittest.TestCase
    % TEST_LOAD_DL_PROVENANCE Unit tests for load_dl_provenance.
    %
    % load_dl_provenance retrieves the Deep Learning training provenance
    % manifest (lists of patient IDs used to train DnCNN and IVIMnet
    % models) so the pipeline can detect and prevent data leakage between
    % DL training sets and the analysis cohort.
    %
    % Tests cover:
    %   - No arguments: returns default struct with empty ID lists
    %   - Non-existent manifest file: returns defaults gracefully
    %   - Valid .mat file: loads dncnn_train_ids and ivimnet_train_ids
    %   - Partial .mat file (missing one field): fills missing field with default
    %   - Explicit provenance struct passed as 2nd argument (caller pattern)
    %   - Explicit provenance struct with non-empty file path (base pattern)

    properties
        TempManifestFile   % Path to a temp .mat file used as mock manifest
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempManifestFile = [tempname '.mat'];

            current_dir = fileparts(mfilename('fullpath'));
            utils_dir = fullfile(current_dir, '..', 'utils');
            addpath(utils_dir);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            % Delete the temp manifest file if it was created
            if exist(testCase.TempManifestFile, 'file')
                delete(testCase.TempManifestFile);
            end

            % Clean up base workspace variable to avoid cross-test contamination
            evalin('base', 'clear dl_provenance_workspace;');
        end
    end

    methods(Test)
        function testNoArgumentsReturnsDefaults(testCase)
            % Calling with no arguments should return a struct with both
            % dncnn_train_ids and ivimnet_train_ids fields, both empty.
            dl_provenance = load_dl_provenance();

            testCase.verifyTrue(isfield(dl_provenance, 'dncnn_train_ids'));
            testCase.verifyTrue(isfield(dl_provenance, 'ivimnet_train_ids'));
            testCase.verifyEmpty(dl_provenance.dncnn_train_ids);
            testCase.verifyEmpty(dl_provenance.ivimnet_train_ids);
        end

        function testMissingFileReturnsDefaults(testCase)
            % A non-existent file path should not error; the function
            % should fall back to the default empty-ID struct.
            dl_provenance = load_dl_provenance('non_existent_file.mat');

            testCase.verifyTrue(isfield(dl_provenance, 'dncnn_train_ids'));
            testCase.verifyTrue(isfield(dl_provenance, 'ivimnet_train_ids'));
            testCase.verifyEmpty(dl_provenance.dncnn_train_ids);
            testCase.verifyEmpty(dl_provenance.ivimnet_train_ids);
        end

        function testLoadFromFile(testCase)
            % Verifies the primary loading path: a .mat file with both
            % dncnn_train_ids and ivimnet_train_ids should be loaded and
            % returned verbatim.
            dl_provenance.dncnn_train_ids = {'patient1', 'patient2'};
            dl_provenance.ivimnet_train_ids = {'patient3', 'patient4'};

            % Save to temp file
            save(testCase.TempManifestFile, 'dl_provenance');

            % Load using the function
            loaded_provenance = load_dl_provenance(testCase.TempManifestFile);

            % Verify
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'patient1', 'patient2'});
            testCase.verifyEqual(loaded_provenance.ivimnet_train_ids, {'patient3', 'patient4'});
        end

        function testLoadFromFileMissingFields(testCase)
            % When the .mat file contains dncnn_train_ids but is missing
            % ivimnet_train_ids, the function should fill the missing field
            % with an empty default rather than erroring.
            dl_provenance.dncnn_train_ids = {'patient1', 'patient2'};
            % ivimnet_train_ids is missing

            % Save to temp file
            save(testCase.TempManifestFile, 'dl_provenance');

            % Load using the function
            loaded_provenance = load_dl_provenance(testCase.TempManifestFile);

            % Verify populated defaults
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'patient1', 'patient2'});
            testCase.verifyTrue(isfield(loaded_provenance, 'ivimnet_train_ids'));
            testCase.verifyEmpty(loaded_provenance.ivimnet_train_ids);
        end

        function testLoadFromCallerWorkspace(testCase)
            % Tests the explicit parameter-passing interface: when a
            % provenance struct is supplied as the 2nd argument with an
            % empty file path, it should be used directly without disk I/O.
            dl_provenance_workspace.dncnn_train_ids = {'caller_pat1'};
            dl_provenance_workspace.ivimnet_train_ids = {'caller_pat2'};

            % Pass provenance struct as explicit second argument
            loaded_provenance = load_dl_provenance('', dl_provenance_workspace);

            % Verify
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'caller_pat1'});
            testCase.verifyEqual(loaded_provenance.ivimnet_train_ids, {'caller_pat2'});
        end

        function testLoadFromBaseWorkspace(testCase)
            % When both a file path and a provenance struct are provided,
            % the explicit struct takes precedence over the file.  This
            % verifies that the 2nd argument always wins.
            base_provenance.dncnn_train_ids = {'base_pat1'};
            base_provenance.ivimnet_train_ids = {'base_pat2'};

            % Pass provenance struct as explicit second argument
            loaded_provenance = load_dl_provenance('some_file.mat', base_provenance);

            % Verify
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'base_pat1'});
            testCase.verifyEqual(loaded_provenance.ivimnet_train_ids, {'base_pat2'});
        end
    end
end
