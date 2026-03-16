"""Report sections: statistics group.

Functions have been split into submodules; this file re-exports them
for backward compatibility.
"""

from .statistics_effects import (  # noqa: F401
    _section_effect_sizes,
    _section_multiple_comparisons,
)
from .statistics_diagnostics import (  # noqa: F401
    _section_model_diagnostics,
)
from .statistics_robustness import (  # noqa: F401
    _section_sensitivity_analysis,
    _section_power_analysis,
)
