"""Tests for generate_report.py — HTML report generation.

This file has been split into focused modules:

- test_generate_report_formatting.py: helper/utility function tests
  (sig_tag, section, forest_plot, effect_size, copy helpers, series
  normalisation, figure/table numbering)

- test_generate_report_sections.py: section builder tests
  (data completeness, feature overlap, power analysis, patient flow,
  sensitivity analysis, feature stability, appendix)

- test_generate_report_manuscript.py: manuscript/publication tests
  (manuscript findings, reporting checklist, results draft, journal
  guide, effect sizes in findings)

- test_generate_report_figures.py: figure gallery tests
  (PNG embedding, DWI type grouping)

- test_generate_report_integration.py: full report integration tests
  (end-to-end generation, new section integration, BibTeX export,
  data quality, Cox PH direction, correlations context)

All tests are discovered automatically by pytest from their new
locations. This file is kept as a pointer for documentation purposes.
"""
