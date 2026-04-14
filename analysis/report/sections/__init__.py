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
from .main_results_summary import _section_executive_summary  # noqa: F401
from .main_results_hypothesis import _section_hypothesis  # noqa: F401
from .main_results_trends import _section_treatment_response  # noqa: F401
from .statistical_reporting import (  # noqa: F401
    _section_statistical_significance,
    _section_broad_statistical_overview,
)
from .manuscript_performance import _section_predictive_performance  # noqa: F401
from .manuscript_findings import _section_manuscript_ready_findings  # noqa: F401
from .manuscript_results import _section_results_draft  # noqa: F401
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
    build_dca_section,
)
from .model_robustness import (  # noqa: F401
    _section_model_robustness,
)
from .statistics_robustness import (  # noqa: F401
    build_nri_idi_section,
)
from .analysis_features import (  # noqa: F401
    build_texture_section,
)
from .data_quality import (  # noqa: F401
    build_registration_quality_section,
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
from .cross_pipeline_dice import (  # noqa: F401
    _section_cross_pipeline_dice,
)
from .failure_rates import (  # noqa: F401
    _section_failure_rates,
)
from .pruning_results import (  # noqa: F401
    _section_pruning_results,
)
from .core_method_outcomes import (  # noqa: F401
    _section_core_method_outcomes,
)
from .forest_plot import (  # noqa: F401
    _section_forest_plot_figure,
    extract_hr_data,
    generate_forest_plot,
)
from .subvolume_stability import (  # noqa: F401
    _section_subvolume_stability,
)
from .per_method_cor import (  # noqa: F401
    _section_per_method_cor,
)
from .repeatability_dice import (  # noqa: F401
    _section_repeatability_dice,
)
from .dose_response_roc import (  # noqa: F401
    _section_dose_response_roc,
)
from .gtv_confounding import (  # noqa: F401
    _section_gtv_confounding,
)
from .risk_dose_concordance import (  # noqa: F401
    _section_risk_dose_concordance,
)
from .dose_context import (  # noqa: F401
    _section_dose_context,
)
from .subvolume_sizes import (  # noqa: F401
    _section_subvolume_sizes,
)
from .threshold_optimization import (  # noqa: F401
    _section_threshold_optimization,
)
from .baseline_vs_delta import (  # noqa: F401
    _section_baseline_vs_delta,
)
