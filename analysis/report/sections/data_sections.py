"""Report sections: data sections group (cohort, flow, completeness, MAT data).

This module re-exports all data section functions from their new submodules
for backward compatibility.
"""

from __future__ import annotations

from .data_overview import (  # noqa: F401
    _section_cohort_overview,
    _section_patient_flow,
)
from .data_quality import (  # noqa: F401
    _section_data_completeness,
)
from .data_supplemental import (  # noqa: F401
    _section_mat_data,
)
