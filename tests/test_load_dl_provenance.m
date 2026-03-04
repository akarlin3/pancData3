classdef test_load_dl_provenance < matlab.unittest.TestCase
    % TEST_LOAD_DL_PROVENANCE Unit tests for loading DL provenance manifests

    properties
        TempManifestFile
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempManifestFile = [tempname '.mat'];

            % Add utils directory to path
            current_dir = fileparts(mfilename('fullpath'));
            utils_dir = fullfile(current_dir, '..', 'utils');
            addpath(utils_dir);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            if exist(testCase.TempManifestFile, 'file')
                delete(testCase.TempManifestFile);
            end

            % Clean up base workspace variable if it exists
            evalin('base', 'clear dl_provenance_workspace;');
        end
    end

    methods(Test)
        function testNoArgumentsReturnsDefaults(testCase)
            % Test when calling with no arguments
            dl_provenance = load_dl_provenance();

            testCase.verifyTrue(isfield(dl_provenance, 'dncnn_train_ids'));
            testCase.verifyTrue(isfield(dl_provenance, 'ivimnet_train_ids'));
            testCase.verifyEmpty(dl_provenance.dncnn_train_ids);
            testCase.verifyEmpty(dl_provenance.ivimnet_train_ids);
        end

        function testMissingFileReturnsDefaults(testCase)
            % Test passing a non-existent file path
            dl_provenance = load_dl_provenance('non_existent_file.mat');

            testCase.verifyTrue(isfield(dl_provenance, 'dncnn_train_ids'));
            testCase.verifyTrue(isfield(dl_provenance, 'ivimnet_train_ids'));
            testCase.verifyEmpty(dl_provenance.dncnn_train_ids);
            testCase.verifyEmpty(dl_provenance.ivimnet_train_ids);
        end

        function testLoadFromFile(testCase)
            % Create a mock .mat file containing a dl_provenance struct
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
            % Create a mock .mat file missing one of the ID fields
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
            % Test explicit parameter passing (replaces old evalin('caller') pattern)
            dl_provenance_workspace.dncnn_train_ids = {'caller_pat1'};
            dl_provenance_workspace.ivimnet_train_ids = {'caller_pat2'};

            % Pass provenance struct as explicit second argument
            loaded_provenance = load_dl_provenance('', dl_provenance_workspace);

            % Verify
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'caller_pat1'});
            testCase.verifyEqual(loaded_provenance.ivimnet_train_ids, {'caller_pat2'});
        end

        function testLoadFromBaseWorkspace(testCase)
            % Test explicit parameter passing (replaces old evalin('base') pattern)
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
