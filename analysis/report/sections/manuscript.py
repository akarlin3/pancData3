"""Report sections: manuscript-ready outputs and predictive performance.

This module re-exports functions that have been split into dedicated
submodules for maintainability.  All original imports continue to work.
"""

from .manuscript_performance import _section_predictive_performance  # noqa: F401
from .manuscript_findings import _section_manuscript_ready_findings  # noqa: F401
from .manuscript_results import _section_results_draft  # noqa: F401

# Configuration for statistical significance thresholds
SIGNIFICANCE_THRESHOLDS = {
    'alpha': 0.05,
    'bonferroni_alpha': 0.01,
    'strict_alpha': 0.001,
    'trending': 0.10
}

# Effect size interpretation constants
EFFECT_SIZE_INTERPRETATIONS = {
    'cohens_d': {
        'small': 0.2,
        'medium': 0.5,
        'large': 0.8
    },
    'eta_squared': {
        'small': 0.01,
        'medium': 0.06,
        'large': 0.14
    },
    'r_squared': {
        'small': 0.02,
        'medium': 0.13,
        'large': 0.26
    },
    'correlation': {
        'small': 0.1,
        'medium': 0.3,
        'large': 0.5
    }
}

def get_significance_threshold(threshold_type='alpha'):
    """Get configurable significance threshold."""
    return SIGNIFICANCE_THRESHOLDS.get(threshold_type, 0.05)

def interpret_effect_size(effect_size_value, effect_type='cohens_d'):
    """Interpret effect size magnitude based on conventional thresholds."""
    thresholds = EFFECT_SIZE_INTERPRETATIONS.get(effect_type, EFFECT_SIZE_INTERPRETATIONS['cohens_d'])
    abs_value = abs(effect_size_value)
    
    if abs_value >= thresholds['large']:
        return 'large'
    elif abs_value >= thresholds['medium']:
        return 'medium'
    elif abs_value >= thresholds['small']:
        return 'small'
    else:
        return 'negligible'