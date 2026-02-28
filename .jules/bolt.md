## 2024-03-01 - Fast AUC calculation via Mann-Whitney U
**Learning:** `perfcurve` in MATLAB is extremely slow when only the AUC is needed because it computes the entire ROC curve across all thresholds.
**Action:** Replace `perfcurve` with a fast AUC computation using the Mann-Whitney U statistic via the `tiedrank` function when the full curve is not required (e.g., in loop-heavy operations like feature selection/filtering).
