---
description: Runs the DWI Workflow sequentially for all DWI types
---

1. Verify that `config.json` exists in the workspace root.
2. Run tests using the `mcp_matlab-mcp_run_matlab_file` tool on `c:\Users\karlina1\Desktop\pancData3\run_all_tests.m`.
3. Update `config.json` to set `"dwi_type": "Standard"` and `"skip_to_reload": false`.
4. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run each step of the pipeline individually:
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
5. Update `config.json` to set `"skip_to_reload": true` and `"dwi_type": "dnCNN"`.
6. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run each step of the pipeline individually again (load, sanity, visualize, metrics_baseline, metrics_longitudinal, metrics_dosimetry, metrics_stats_comparisons, metrics_stats_predictive, metrics_survival).
7. Update `config.json` to set `"dwi_type": "IVIMnet"`.
8. Use the `mcp_matlab-mcp_evaluate_matlab_code` tool to run each step of the pipeline individually again (load, sanity, visualize, metrics_baseline, metrics_longitudinal, metrics_dosimetry, metrics_stats_comparisons, metrics_stats_predictive, metrics_survival).
