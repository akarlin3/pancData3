"""Report sections: analysis group (graphs, cross-DWI, features).

Functions have been split into submodules; this file re-exports them
for backward compatibility.
"""

from .analysis_graphs import (  # noqa: F401
    _section_graph_overview,
    _section_graph_issues,
    _section_stats_by_graph_type,
)
from .analysis_cross_dwi import (  # noqa: F401
    _section_cross_dwi_comparison,
    _section_correlations,
)
from .analysis_features import (  # noqa: F401
    _section_feature_overlap,
)
