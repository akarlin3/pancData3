## 2024-05-24 - Fix Insecure Deserialization in MATLAB `load()` Calls
**Vulnerability:** Transparent `load()` calls in multiple MATLAB scripts (`run_dwi_pipeline.m`, `core/load_dwi_data.m`, `core/compute_summary_metrics.m`).
**Learning:** Calling `load(filename)` without an output argument deserializes arbitrary variables from the `.mat` file directly into the local workspace, leading to potential variable clobbering, injection, or unexpected code execution.
**Prevention:** Use structural assignments (e.g., `tmp = load(filename); var = tmp.var;`) when loading variables from `.mat` files to restrict extraction to explicitly defined fields and prevent uncontrolled workspace pollution.
