---
description: Runs the DWI Workflow sequentially for all DWI types
---

1. Verify that `config.json` exists in the workspace root.
2. Update `config.json` to set `"dwi_type": "Standard"`.
3. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run the pipeline:
   - `project_path`: `c:\Users\karlina1\Desktop\pancData3`
   - `code`: `run_dwi_pipeline('config.json')`
4. Update `config.json` to set `"dwi_type": "dnCNN"`.
5. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run the pipeline again.
6. Update `config.json` to set `"dwi_type": "IVIMnet"`.
7. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run the pipeline again.
8. Run tests using the `mcp_matlab-mcp_run_matlab_file` tool on `c:\Users\karlina1\Desktop\pancData3\run_all_tests.m`.