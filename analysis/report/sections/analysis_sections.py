"""Report sections: analysis sections group.

This module re-exports functions that have been split into:
- ``graph_overview`` (graph overview, issues, stats by type)
- ``cross_dwi`` (cross-DWI comparison, feature overlap)
- ``correlations`` (notable correlations)

Import directly from the submodules for new code.
"""

from report.sections.graph_overview import (  # noqa: F401
    _section_graph_overview,
    _section_graph_issues,
    _section_stats_by_graph_type,
)
from report.sections.cross_dwi import (  # noqa: F401
    _section_cross_dwi_comparison,
    _section_feature_overlap,
)
from report.sections.correlations import (  # noqa: F401
    _section_correlations,
)
