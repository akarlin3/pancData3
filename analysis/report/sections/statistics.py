"""Report sections: statistics group.

This module re-exports functions that have been split into:
- ``effect_sizes`` (effect size analysis, multiple comparisons)
- ``model_diagnostics`` (model diagnostics, sensitivity analysis)
- ``power_analysis`` (statistical power commentary)

Import directly from the submodules for new code.
"""

from report.sections.effect_sizes import (  # noqa: F401
    _section_effect_sizes,
    _section_multiple_comparisons,
)
from report.sections.model_diagnostics import (  # noqa: F401
    _section_model_diagnostics,
    _section_sensitivity_analysis,
)
from report.sections.power_analysis import (  # noqa: F401
    _section_power_analysis,
)
