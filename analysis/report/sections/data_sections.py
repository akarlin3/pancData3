"""Report sections: data sections group (cohort, flow, completeness, MAT data).

This module re-exports functions that have been split into:
- ``enrollment`` (cohort overview, patient flow, data completeness)
- ``supplemental`` (MAT file data: dosimetry, core method comparison)

Import directly from the submodules for new code.
"""

from report.sections.enrollment import (  # noqa: F401
    _section_cohort_overview,
    _section_patient_flow,
    _section_data_completeness,
)
from report.sections.supplemental import (  # noqa: F401
    _section_mat_data,
)
