"""Tests for cross-DWI agreement analysis (Bland-Altman, CCC, ICC)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cross_reference.cross_dwi_agreement import (
    bland_altman,
    icc_absolute,
    lins_ccc,
)


class TestBlandAltman:
    """Tests for Bland-Altman analysis."""

    def test_perfect_agreement(self):
        """Perfect agreement should have zero bias."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = bland_altman(x, y)
        assert result["bias"] == pytest.approx(0.0)
        assert result["sd"] == pytest.approx(0.0, abs=1e-10)
        assert result["n"] == 5

    def test_known_bias(self):
        """Systematic bias should appear in the bias statistic."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = x + 0.5  # systematic positive bias
        result = bland_altman(x, y)
        assert result["bias"] == pytest.approx(-0.5)

    def test_limits_of_agreement(self):
        """LoA should be bias ± 1.96 * SD."""
        rng = np.random.default_rng(42)
        x = rng.normal(10, 2, size=100)
        y = x + rng.normal(0, 0.5, size=100)
        result = bland_altman(x, y)
        expected_loa_width = 2 * 1.96 * result["sd"]
        actual_loa_width = result["loa_hi"] - result["loa_lo"]
        assert actual_loa_width == pytest.approx(expected_loa_width, rel=1e-10)

    def test_pct_outside(self):
        """Percentage outside LoA should be ~5% for normal differences."""
        rng = np.random.default_rng(42)
        n = 1000
        x = rng.normal(0, 1, size=n)
        y = x + rng.normal(0, 1, size=n)
        result = bland_altman(x, y)
        # Expect roughly 5% outside LoA (with some tolerance)
        assert 1.0 < result["pct_outside_loa"] < 10.0

    def test_missing_handling(self):
        """Empty arrays should return zero n."""
        result = bland_altman(np.array([]), np.array([]))
        assert result["n"] == 0


class TestLinsCCC:
    """Tests for Lin's concordance correlation coefficient."""

    def test_perfect_agreement(self):
        """Perfect agreement should give CCC = 1."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = lins_ccc(x, x)
        assert result["ccc"] == pytest.approx(1.0, abs=1e-10)

    def test_perfect_negative(self):
        """Perfect negative linear should give CCC = -1."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = -x + 6  # perfect negative linear (same scale)
        result = lins_ccc(x, y)
        assert result["ccc"] == pytest.approx(-1.0, abs=0.01)

    def test_systematic_bias_reduces_ccc(self):
        """Systematic bias (same correlation, different means) reduces CCC."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y_perfect = x.copy()
        y_biased = x + 10.0  # large systematic bias

        ccc_perfect = lins_ccc(x, y_perfect)["ccc"]
        ccc_biased = lins_ccc(x, y_biased)["ccc"]
        assert ccc_perfect > ccc_biased

    def test_ci_contains_point_estimate(self):
        """CI should contain the point estimate."""
        rng = np.random.default_rng(42)
        x = rng.normal(0, 1, size=50)
        y = x + rng.normal(0, 0.3, size=50)
        result = lins_ccc(x, y)
        assert result["ci_lo"] <= result["ccc"] <= result["ci_hi"]

    def test_few_points(self):
        """With <3 points, should return NaN."""
        result = lins_ccc(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
        assert math.isnan(result["ccc"])


class TestICC:
    """Tests for intraclass correlation coefficient."""

    def test_perfect_agreement(self):
        """Perfect agreement should give ICC ≈ 1."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = icc_absolute(x, x)
        assert result == pytest.approx(1.0, abs=0.01)

    def test_random_disagreement(self):
        """Random data should give ICC near 0."""
        rng = np.random.default_rng(42)
        x = rng.normal(0, 1, size=100)
        y = rng.normal(0, 1, size=100)
        result = icc_absolute(x, y)
        assert abs(result) < 0.3

    def test_few_points(self):
        """With <3 points, should return NaN."""
        result = icc_absolute(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
        assert math.isnan(result)
