"""Report section builders — split into logical submodules.

All public section functions are re-exported here for backward
compatibility with existing ``from report_sections import …`` usage.
"""

from .metadata import (  # noqa: F401
    _section_publication_header,
    _section_data_availability,
    _section_table_index,
    _section_figure_index,
)
from .main_results import (  # noqa: F401
    _section_executive_summary,
    _section_hypothesis,
    _section_statistical_significance,
    _section_treatment_response,
    _section_predictive_performance,
    _section_manuscript_ready_findings,
    _section_results_draft,
)
from .data_sections import (  # noqa: F401
    _section_cohort_overview,
    _section_patient_flow,
    _section_data_completeness,
    _section_mat_data,
    _section_appendix,
    _section_figure_gallery,
)
from .analysis_sections import (  # noqa: F401
    _section_graph_overview,
    _section_graph_issues,
    _section_stats_by_graph_type,
    _section_cross_dwi_comparison,
    _section_correlations,
    _section_feature_overlap,
)
from .statistics import (  # noqa: F401
    _section_effect_sizes,
    _section_multiple_comparisons,
    _section_model_diagnostics,
    _section_sensitivity_analysis,
    _section_power_analysis,
)
from .discussion import (  # noqa: F401
    _section_methods,
    _section_limitations,
    _section_conclusions,
    _section_reporting_checklist,
    _section_journal_guide,
)
