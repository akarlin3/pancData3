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
            % Define variable in the caller workspace (this function's scope)
            dl_provenance_workspace.dncnn_train_ids = {'caller_pat1'};
            dl_provenance_workspace.ivimnet_train_ids = {'caller_pat2'};

            % Assign to caller workspace is tricky in matlab tests directly,
            % but since the function checks evalin('caller', ...), and the test method
            % itself is the caller, having it as a local variable here works.

            % However, to ensure it works properly as intended by 'caller' evalin,
            % we create a local helper function
            loaded_provenance = testCase.callLoadHelper();

            % Verify
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'caller_pat1'});
            testCase.verifyEqual(loaded_provenance.ivimnet_train_ids, {'caller_pat2'});
        end

        function testLoadFromBaseWorkspace(testCase)
            % Define variable in the base workspace
            base_provenance.dncnn_train_ids = {'base_pat1'};
            base_provenance.ivimnet_train_ids = {'base_pat2'};
            assignin('base', 'dl_provenance_workspace', base_provenance);

            % Load using the function
            loaded_provenance = load_dl_provenance('some_file.mat');

            % Verify
            testCase.verifyEqual(loaded_provenance.dncnn_train_ids, {'base_pat1'});
            testCase.verifyEqual(loaded_provenance.ivimnet_train_ids, {'base_pat2'});
        end
    end

    methods(Access=private)
        function loaded_provenance = callLoadHelper(testCase)
            % Define the variable in the caller's workspace (this helper's scope)
            dl_provenance_workspace.dncnn_train_ids = {'caller_pat1'};
            dl_provenance_workspace.ivimnet_train_ids = {'caller_pat2'};

            % Call the target function, which will look at this workspace
            loaded_provenance = load_dl_provenance();
        end
    end
end
