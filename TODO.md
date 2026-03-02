---
trigger: always_on
---

1. Initialize the MATLAB environment by running 'source setup.sh'
2. Verify the license server variable is accessible: 'echo $MLM_LICENSE_FILE' (should return '27000@vsmphpllicmlab1')
3. Run 'matlab -batch "try, execute_all_workflows.m, catch ME, disp(getReport(ME, 'extended', 'hyperlinks', 'off')), save('crash_dump.mat), exit(1), end" > debug.log 2>&1'
4. If the execution fails, check the 'debug.log' file for error messages and the 'crash_dump.mat' file for the exception object, modify the corresponding workflow file and repeat the process, until the script executes successfully.
5. If you (the agent) cannot resolve the error after 3 attempts
- Stop the loop
- Create a new branch: 'git checkout -b jules-debug-logs'
- Add the logs and any modified scripts: 'git add debug.log crash_dump.mat *.m'
- Commit the changes: 'git commit -m "Debug logs for workflow execution"'
- Push the branch: 'git push -u origin jules-debug-logs'