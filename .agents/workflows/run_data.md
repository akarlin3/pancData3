---
description: Runs the DWI Workflow sequentially for all DWI types
---

1. Verify that `config.json` exists in the workspace root.
2. Run tests using the `mcp_matlab-mcp_run_matlab_test_file` tool on `tests/test_dwi_pipeline.m` to ensure the codebase is stable.
// turbo
3. Execute the full sequential pipeline orchestration using a background batch command in the terminal to avoid worker sync and timeout issues. Use the `run_command` tool with `CommandLine`: `matlab -batch "execute_all_workflows"`, and `Cwd`: `<repository_root>`. Set `SafeToAutoRun` to true, and `WaitMsBeforeAsync` to 5000.
4. Wait for the background command to finish using the `command_status` tool. Once done, verify the generated summary metric MAT files and output figures.
