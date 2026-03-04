
## 2024-11-20 - P-Value Formatting in Scientific Visualizations
**Learning:** In scientific and medical data visualization, formatting p-values directly with a fixed float precision (like `sprintf('p = %.3f', p)`) creates a poor user experience for highly significant results by displaying `p = 0.000`. This is mathematically misleading and violates APA/scientific reporting standards.
**Action:** Always wrap p-value display logic in a formatter function (e.g., `format_p_value`) that intercepts values smaller than the display precision limit (e.g., `< 0.001`) and returns a less-than string (`p < 0.001`) to maintain UI accuracy and trust.
