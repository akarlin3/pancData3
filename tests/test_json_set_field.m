classdef test_json_set_field < matlab.unittest.TestCase
    % TEST_JSON_SET_FIELD Unit tests for json_set_field utility

    methods(TestMethodSetup)
        function setup(~)
            current_dir = fileparts(mfilename('fullpath'));
            addpath(fullfile(current_dir, '..', 'utils'));
        end
    end

    methods(Test)
        function testReplaceStringValue(testCase)
            s = '{"dwi_type": "Standard", "other": 1}';
            result = json_set_field(s, 'dwi_type', 'dnCNN');
            testCase.verifySubstring(result, '"dnCNN"');
            testCase.verifySubstring(result, '"other": 1');
        end

        function testReplaceBooleanTrue(testCase)
            s = '{"skip_to_reload": false}';
            result = json_set_field(s, 'skip_to_reload', true);
            testCase.verifySubstring(result, 'true');
            testCase.verifyTrue(~contains(result, 'false'));
        end

        function testReplaceBooleanFalse(testCase)
            s = '{"skip_to_reload": true}';
            result = json_set_field(s, 'skip_to_reload', false);
            testCase.verifySubstring(result, 'false');
        end

        function testPreservesFormatting(testCase)
            s = sprintf('{\n    "dwi_type": "Standard",\n    "dataloc": "/data/path",\n    "skip_to_reload": false\n}');
            result = json_set_field(s, 'dwi_type', 'IVIMnet');
            result = json_set_field(result, 'skip_to_reload', true);
            % Verify indentation and newlines are preserved
            testCase.verifySubstring(result, sprintf('    "dwi_type": "IVIMnet"'));
            testCase.verifySubstring(result, sprintf('    "dataloc": "/data/path"'));
            testCase.verifySubstring(result, sprintf('    "skip_to_reload": true'));
        end

        function testPreservesFieldOrder(testCase)
            s = '{"b": 2, "a": 1}';
            result = json_set_field(s, 'a', 99);
            % "b" should still appear before "a"
            b_pos = strfind(result, '"b"');
            a_pos = strfind(result, '"a"');
            testCase.verifyLessThan(b_pos(1), a_pos(1));
        end

        function testMissingFieldErrors(testCase)
            s = '{"dwi_type": "Standard"}';
            testCase.verifyError(@() json_set_field(s, 'nonexistent', 'val'), ...
                'json_set_field:fieldNotFound');
        end

        function testNumericValue(testCase)
            s = '{"ivim_bthr": 100}';
            result = json_set_field(s, 'ivim_bthr', 200);
            testCase.verifySubstring(result, '200');
        end

        function testPreservesEmptyStringFields(testCase)
            % Verifies that fields with "" values (like cause_of_death_column)
            % are not corrupted when other fields are modified.
            s = sprintf('{\n    "dwi_type": "Standard",\n    "cause_of_death_column": "",\n    "skip_to_reload": false\n}');
            result = json_set_field(s, 'dwi_type', 'dnCNN');
            result = json_set_field(result, 'skip_to_reload', true);
            testCase.verifySubstring(result, '"cause_of_death_column": ""');
        end

        function testRoundTripPreservesBytes(testCase)
            % Simulate the execute_all_workflows pattern: modify fields
            % then restore original. The original string must be unchanged.
            original = sprintf('{\n    "dwi_type": "Standard",\n    "skip_to_reload": false,\n    "dataloc": "/my/path"\n}');
            modified = json_set_field(original, 'dwi_type', 'dnCNN');
            modified = json_set_field(modified, 'skip_to_reload', true);
            restored = json_set_field(modified, 'dwi_type', 'Standard');
            restored = json_set_field(restored, 'skip_to_reload', false);
            testCase.verifyEqual(restored, original);
        end
    end
end
