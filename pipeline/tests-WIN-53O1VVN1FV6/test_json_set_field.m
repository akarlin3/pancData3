classdef test_json_set_field < matlab.unittest.TestCase
    % TEST_JSON_SET_FIELD Unit tests for the json_set_field utility.
    %
    % json_set_field performs targeted regex-based replacement of a single
    % field's value within a raw JSON string, without parsing the entire JSON
    % into a struct (which would lose formatting, field order, and comments).
    %
    % This utility is critical for execute_all_workflows.m, which modifies
    % config.json in-place between DWI type runs (e.g., changing dwi_type
    % from "Standard" to "dnCNN") while preserving the file's exact
    % formatting and field order.
    %
    % Tests cover:
    %   - String, boolean, and numeric value replacement
    %   - Formatting and indentation preservation
    %   - Field order preservation (no alphabetic re-sorting)
    %   - Error on missing field name
    %   - Empty string fields are not corrupted by adjacent edits
    %   - Round-trip: modify then restore produces identical output
    %
    % Run tests with:
    %   results = runtests('tests/test_json_set_field.m');

    methods(TestMethodSetup)
        function setup(~)
            current_dir = fileparts(mfilename('fullpath'));
            addpath(fullfile(current_dir, '..', 'utils'));
        end
    end

    methods(Test)
        function testReplaceStringValue(testCase)
            % Replace a string value while leaving other fields untouched.
            s = '{"dwi_type": "Standard", "other": 1}';
            result = json_set_field(s, 'dwi_type', 'dnCNN');
            testCase.verifySubstring(result, '"dnCNN"');
            testCase.verifySubstring(result, '"other": 1');
        end

        function testReplaceBooleanTrue(testCase)
            % Replace false with true. The word "false" must no longer
            % appear anywhere in the output.
            s = '{"skip_to_reload": false}';
            result = json_set_field(s, 'skip_to_reload', true);
            testCase.verifySubstring(result, 'true');
            testCase.verifyTrue(~contains(result, 'false'));
        end

        function testReplaceBooleanFalse(testCase)
            % Replace true with false (the reverse direction).
            s = '{"skip_to_reload": true}';
            result = json_set_field(s, 'skip_to_reload', false);
            testCase.verifySubstring(result, 'false');
        end

        function testPreservesFormatting(testCase)
            % Multi-line JSON with 4-space indentation. After modifying two
            % fields, the indentation, newlines, and unmodified fields must
            % remain byte-identical to the original formatting.
            s = sprintf('{\n    "dwi_type": "Standard",\n    "dataloc": "/data/path",\n    "skip_to_reload": false\n}');
            result = json_set_field(s, 'dwi_type', 'IVIMnet');
            result = json_set_field(result, 'skip_to_reload', true);
            testCase.verifySubstring(result, sprintf('    "dwi_type": "IVIMnet"'));
            testCase.verifySubstring(result, sprintf('    "dataloc": "/data/path"'));
            testCase.verifySubstring(result, sprintf('    "skip_to_reload": true'));
        end

        function testPreservesFieldOrder(testCase)
            % Unlike jsondecode/jsonencode, json_set_field must not re-sort
            % fields alphabetically. Field "b" appears before "a" in the
            % input and must remain before "a" in the output.
            s = '{"b": 2, "a": 1}';
            result = json_set_field(s, 'a', 99);
            b_pos = strfind(result, '"b"');
            a_pos = strfind(result, '"a"');
            testCase.verifyLessThan(b_pos(1), a_pos(1));
        end

        function testMissingFieldErrors(testCase)
            % Attempting to set a field that does not exist in the JSON
            % string must throw an error rather than silently failing.
            s = '{"dwi_type": "Standard"}';
            testCase.verifyError(@() json_set_field(s, 'nonexistent', 'val'), ...
                'json_set_field:fieldNotFound');
        end

        function testNumericValue(testCase)
            % Replace a numeric value (integer). The regex must match
            % both string and numeric JSON value patterns.
            s = '{"ivim_bthr": 100}';
            result = json_set_field(s, 'ivim_bthr', 200);
            testCase.verifySubstring(result, '200');
        end

        function testPreservesEmptyStringFields(testCase)
            % Edge case: fields with empty string values ("") can confuse
            % regex patterns. Modifying adjacent fields must not corrupt
            % the empty string value of cause_of_death_column.
            s = sprintf('{\n    "dwi_type": "Standard",\n    "cause_of_death_column": "",\n    "skip_to_reload": false\n}');
            result = json_set_field(s, 'dwi_type', 'dnCNN');
            result = json_set_field(result, 'skip_to_reload', true);
            testCase.verifySubstring(result, '"cause_of_death_column": ""');
        end

        function testRoundTripPreservesBytes(testCase)
            % Simulates the execute_all_workflows pattern: config.json is
            % modified for dnCNN/IVIMnet runs then restored to its original
            % state. After a full round-trip (modify -> restore), the output
            % string must be byte-identical to the original.
            original = sprintf('{\n    "dwi_type": "Standard",\n    "skip_to_reload": false,\n    "dataloc": "/my/path"\n}');
            modified = json_set_field(original, 'dwi_type', 'dnCNN');
            modified = json_set_field(modified, 'skip_to_reload', true);
            restored = json_set_field(modified, 'dwi_type', 'Standard');
            restored = json_set_field(restored, 'skip_to_reload', false);
            testCase.verifyEqual(restored, original);
        end
    end
end
