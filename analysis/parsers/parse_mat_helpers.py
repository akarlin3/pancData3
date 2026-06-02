#!/usr/bin/env python3
"""Shared helpers and optional-dependency shims for MAT-metric parsing.

This module holds the small JSON-safe conversion helpers, the method
description lookup table, and the optional ``scipy`` / ``numpy`` imports
shared by :mod:`parsers.parse_mat_metrics` and its sibling parser modules.
It is a pure data / helper module with no side effects so it can be
imported safely from any of the MAT parsers.

If ``scipy`` and ``numpy`` are not installed, ``scipy_io`` and ``numpy_np``
are set to ``None`` and the callers degrade gracefully (returning empty
data) rather than crashing.
"""

import math
import typing

# Attempt to import scipy and numpy.  These are optional dependencies --
# if not installed, the callers degrade gracefully (returns empty data).
try:
    import scipy.io  # type: ignore
    import numpy as np  # type: ignore
    scipy_io: typing.Any = scipy.io
    numpy_np: typing.Any = np
except ImportError:
    print("Warning: scipy or numpy not found. Please pip install scipy numpy.")
    scipy_io = None
    numpy_np = None


# Known method name → human-readable description mapping.
METHOD_DESC: dict[str, str] = {
    "adc_threshold": "ADC threshold (< threshold mm²/s)",
    "d_threshold": "D threshold (< threshold mm²/s)",
    "df_intersection": "D–f intersection (D low AND f high)",
    "otsu": "Otsu automatic thresholding",
    "gmm": "Gaussian mixture model (2-component)",
    "kmeans": "K-means clustering (k=2)",
    "region_growing": "Region growing from seed voxels",
    "active_contours": "Active contour (snake) segmentation",
    "percentile": "Percentile-based threshold",
    "spectral": "Spectral clustering",
    "fdm": "Functional diffusion map",
}


def _safe_float(v: typing.Any) -> typing.Optional[float]:
    """Convert a value to a JSON-safe float, returning None for NaN/Inf.

    Parameters
    ----------
    v : any
        Input value to convert.

    Returns
    -------
    float or None
        Rounded float, or ``None`` if the value is NaN, Inf, or unconvertible.
    """
    if v is None:
        return None
    try:
        f = float(v)
        if math.isnan(f) or math.isinf(f):
            return None
        return float(f"{f:.4g}")
    except (TypeError, ValueError):
        return None


def _array_to_list(arr: typing.Any) -> typing.Any:
    """Recursively convert a numpy array to a nested list of JSON-safe floats.

    NaN and Inf values are replaced with ``None`` (JSON ``null``).

    Parameters
    ----------
    arr : numpy array or scalar
        The value to convert.

    Returns
    -------
    list or None
        Nested Python lists with ``None`` in place of NaN/Inf.
    """
    if not hasattr(arr, "tolist"):
        return _safe_float(arr)
    raw = arr.tolist()
    return _nested_safe(raw)


def _nested_safe(obj: typing.Any) -> typing.Any:
    """Recursively replace float NaN/Inf with None in nested lists/scalars."""
    if isinstance(obj, list):
        return [_nested_safe(x) for x in obj]
    if isinstance(obj, float):
        return _safe_float(obj)
    return obj
