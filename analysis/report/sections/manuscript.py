"""Report sections: manuscript-ready outputs and predictive performance.

This module re-exports functions that have been split into dedicated
submodules for maintainability.  All original imports continue to work.
"""

from .manuscript_performance import _section_predictive_performance  # noqa: F401
from .manuscript_findings import _section_manuscript_ready_findings  # noqa: F401
from .manuscript_results import _section_results_draft  # noqa: F401
