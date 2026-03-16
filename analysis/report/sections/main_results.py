"""Report sections: main results group.

Functions have been split into submodules; this file re-exports them
for backward compatibility.
"""

from .main_results_summary import _section_executive_summary  # noqa: F401
from .main_results_hypothesis import _section_hypothesis  # noqa: F401
from .main_results_trends import _section_treatment_response  # noqa: F401
