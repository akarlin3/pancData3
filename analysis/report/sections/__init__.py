"""Report section builders — split into logical submodules.

All public section functions are re-exported here for backward
compatibility with existing ``from report_sections import …`` usage.
"""

from .metadata import (  # noqa: F401
    _part_break,
    _section_cover_page,
    _section_print_toc,
    _section_publication_header,
    _section_data_availability,
    _section_table_index,
    _section_figure_index,
)
from .main_results import (  # noqa: F401
    _section_executive_summary,
    _section_hypothesis,
    _section_treatment_response,
)
from .statistical_reporting import (  # noqa: F401
    _section_statistical_significance,
    _section_broad_statistical_overview,
)
from .manuscript import (  # noqa: F401
    _section_predictive_performance,
    _section_manuscript_ready_findings,
    _section_results_draft,
)
from .enrollment import (  # noqa: F401
    _section_cohort_overview,
    _section_patient_flow,
    _section_data_completeness,
)
from .supplemental import (  # noqa: F401
    _section_mat_data,
)
from .gallery import (  # noqa: F401
    _section_appendix,
    _section_figure_gallery,
)
from .graph_overview import (  # noqa: F401
    _section_graph_overview,
    _section_graph_issues,
    _section_stats_by_graph_type,
)
from .cross_dwi import (  # noqa: F401
    _section_cross_dwi_comparison,
    _section_feature_overlap,
)
from .correlations import (  # noqa: F401
    _section_correlations,
)
from .effect_sizes import (  # noqa: F401
    _section_effect_sizes,
    _section_multiple_comparisons,
)
from .model_diagnostics import (  # noqa: F401
    _section_model_diagnostics,
    _section_sensitivity_analysis,
)
from .power_analysis import (  # noqa: F401
    _section_power_analysis,
)
from .discussion import (  # noqa: F401
    _section_methods,
    _section_limitations,
    _section_conclusions,
)
from .publication import (  # noqa: F401
    _section_reporting_checklist,
    _section_journal_guide,
)
