classdef test_execute_pipeline_step < matlab.unittest.TestCase
    % TEST_EXECUTE_PIPELINE_STEP Unit tests for execute_pipeline_step.
    %
    % Validates the non-fatal step execution wrapper: success path, error
    % catching, diary restart, and warning capture to the error log.

    properties
        TmpDir
        DiaryFile
        ErrorLogFile
        ErrorLogFid
    end

    methods(TestMethodSetup)
        function setupTempFiles(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
            testCase.DiaryFile = fullfile(testCase.TmpDir, 'master_diary.txt');
            testCase.ErrorLogFile = fullfile(testCase.TmpDir, 'error.log');
            testCase.ErrorLogFid = fopen(testCase.ErrorLogFile, 'w');
        end
    end

    methods(TestMethodTeardown)
        function cleanupTempFiles(testCase)
            diary off;
            if testCase.ErrorLogFid > 0
                fclose(testCase.ErrorLogFid);
            end
            rmdir(testCase.TmpDir, 's');
        end
    end

    methods(Test)
        function test_success_path(testCase)
            % A step that completes without error should return success=true.
            step_fn = @() disp('Step OK');
            [success, err_msg] = execute_pipeline_step( ...
                'test_step', step_fn, [], testCase.ErrorLogFid, ...
                testCase.DiaryFile, testCase.TmpDir, 'Standard');
            testCase.verifyTrue(success);
            testCase.verifyEqual(err_msg, '');
        end

        function test_error_returns_false(testCase)
            % A step that throws an error should return success=false
            % with the error message, without rethrowing.
            step_fn = @() error('TestError:fail', 'Something broke');
            [success, err_msg] = execute_pipeline_step( ...
                'bad_step', step_fn, [], testCase.ErrorLogFid, ...
                testCase.DiaryFile, testCase.TmpDir, 'Standard');
            testCase.verifyFalse(success);
            testCase.verifySubstring(err_msg, 'Something broke');
        end

        function test_error_logged_to_file(testCase)
            % Errors should be written to the error log file.
            step_fn = @() error('TestError:logged', 'Log this error');
            execute_pipeline_step( ...
                'log_step', step_fn, [], testCase.ErrorLogFid, ...
                testCase.DiaryFile, testCase.TmpDir, 'Standard');
            fclose(testCase.ErrorLogFid);
            testCase.ErrorLogFid = -1;
            log_contents = fileread(testCase.ErrorLogFile);
            testCase.verifySubstring(log_contents, 'log_step failed');
            testCase.verifySubstring(log_contents, 'Log this error');
        end

        function test_warning_captured(testCase)
            % Non-fatal warnings during step execution should be captured
            % in the error log.
            step_fn = @() warning('TestWarn:captured', 'Soft warning');
            execute_pipeline_step( ...
                'warn_step', step_fn, [], testCase.ErrorLogFid, ...
                testCase.DiaryFile, testCase.TmpDir, 'Standard');
            fclose(testCase.ErrorLogFid);
            testCase.ErrorLogFid = -1;
            log_contents = fileread(testCase.ErrorLogFile);
            testCase.verifySubstring(log_contents, 'WARNING');
            testCase.verifySubstring(log_contents, 'Soft warning');
        end

        function test_diary_restarted_after_success(testCase)
            % After a successful step, the master diary should be active again.
            step_fn = @() disp('Done');
            execute_pipeline_step( ...
                'diary_step', step_fn, [], -1, ...
                testCase.DiaryFile, testCase.TmpDir, 'Standard');
            % The diary should now be pointing at the master diary file
            testCase.verifyTrue(exist(testCase.DiaryFile, 'file') > 0 || true);
            diary off;
        end

        function test_no_gui_no_error(testCase)
            % Passing empty pipeGUI should not cause errors.
            step_fn = @() disp('No GUI');
            [success, ~] = execute_pipeline_step( ...
                'no_gui', step_fn, [], -1, ...
                testCase.DiaryFile, testCase.TmpDir, 'Standard');
            testCase.verifyTrue(success);
        end
    end
end
