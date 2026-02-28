---
description: Runs the DWI Workflow sequentially for all DWI types
---

1. Verify that `config.json` exists in the workspace root.
2. Update `config.json` to set `"dwi_type": "Standard"`.
3. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run each step of the pipeline individually:
   - `project_path`: `c:\Users\karlina1\Desktop\pancData3`
   - `code`: `run_dwi_pipeline('config.json', {'load'})` % includes: discover_patient_files, load_dwi_data (extraction), compute_summary_metrics
   - `code`: `run_dwi_pipeline('config.json', {'sanity'})`
   - `code`: `run_dwi_pipeline('config.json', {'visualize'})`
   - `code`: `run_dwi_pipeline('config.json', {'metrics_baseline'})`
   - `code`: `run_dwi_pipeline('config.json', {'metrics_longitudinal'})`
   - `code`: `run_dwi_pipeline('config.json', {'metrics_dosimetry'})`
   - `code`: `run_dwi_pipeline('config.json', {'metrics_stats_comparisons'})`
   - `code`: `run_dwi_pipeline('config.json', {'metrics_stats_predictive'})`
   - `code`: `run_dwi_pipeline('config.json', {'metrics_survival'})`
4. Update `config.json` to set `"dwi_type": "dnCNN"`.
5. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run each step of the pipeline individually again (load, sanity, visualize, metrics_baseline, metrics_longitudinal, metrics_dosimetry, metrics_stats_comparisons, metrics_stats_predictive, metrics_survival).
6. Update `config.json` to set `"dwi_type": "IVIMnet"`.
7. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run each step of the pipeline individually again (load, sanity, visualize, metrics_baseline, metrics_longitudinal, metrics_dosimetry, metrics_stats_comparisons, metrics_stats_predictive, metrics_survival).
8. Run tests using the `mcp_matlab-mcp_run_matlab_file` tool on `c:\Users\karlina1\Desktop\pancData3\run_all_tests.m`.