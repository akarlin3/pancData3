---
description: Runs the DWI Workflow
---

1. Verify that `config.json` exists in the workspace root. 
2. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run the pipeline.
   - `project_path`: `c:\Users\karlina1\Desktop\pancData3`
   - `code`: `run_dwi_pipeline('config.json')`
3. Run tests using the `mcp_matlab-mcp_run_matlab_file` tool on `c:\Users\karlina1\Desktop\pancData3\run_all_tests.m`.