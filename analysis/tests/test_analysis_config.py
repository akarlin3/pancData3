"""Tests for the centralised analysis configuration system.

Covers:
- _deep_merge: recursive dict merging
- load_analysis_config: layered loading (defaults → file → MATLAB config)
- get_config / reset_config_cache: caching behaviour
- CLI override integration via run_analysis arg parsing
- Backward compatibility: DWI_TYPES module constant still works
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from shared import (
    DWI_TYPES,
    _DEFAULTS,
    _deep_merge,
    get_config,
    load_analysis_config,
    reset_config_cache,
)


@pytest.fixture(autouse=True)
def _clean_config_cache():
    """Reset the config cache before and after every test."""
    reset_config_cache()
    yield
    reset_config_cache()


# ---------------------------------------------------------------------------
# _deep_merge
# ---------------------------------------------------------------------------

class TestDeepMerge:
    """Verify recursive dict merging behaviour."""

    def test_scalar_override(self):
        """Top-level scalar values in override replace base values."""
        base = {"a": 1, "b": 2}
        override = {"b": 99}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": 99}

    def test_nested_partial_override(self):
        """Nested dicts are merged recursively, not replaced wholesale."""
        base = {"vision": {"model": "old", "retries": 4}}
        override = {"vision": {"model": "new"}}
        result = _deep_merge(base, override)
        assert result["vision"]["model"] == "new"
        assert result["vision"]["retries"] == 4

    def test_new_keys_added(self):
        """Keys present in override but not in base are added."""
        base = {"a": 1}
        override = {"b": 2}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": 2}

    def test_does_not_mutate_base(self):
        """The original base dict is not modified."""
        base = {"a": {"x": 1}}
        _deep_merge(base, {"a": {"x": 99}})
        assert base["a"]["x"] == 1

    def test_list_replaced_not_merged(self):
        """Lists are replaced entirely (no element-wise merging)."""
        base = {"items": [1, 2, 3]}
        override = {"items": [4, 5]}
        result = _deep_merge(base, override)
        assert result["items"] == [4, 5]


# ---------------------------------------------------------------------------
# load_analysis_config
# ---------------------------------------------------------------------------

class TestLoadAnalysisConfig:
    """Verify layered config loading."""

    def test_defaults_returned_when_no_files(self, tmp_path: Path):
        """With no config files, built-in defaults are returned."""
        cfg = load_analysis_config(
            config_path=tmp_path / "nonexistent.json",
            matlab_config_path=tmp_path / "nonexistent.json",
        )
        assert cfg["dwi_types"] == ["Standard", "dnCNN", "IVIMnet"]
        assert cfg["vision"]["gemini_model"] == "gemini-3.5-pro"
        assert cfg["statistics"]["p_noteworthy"] == 0.05

    def test_analysis_config_overrides_defaults(self, tmp_path: Path):
        """Values in analysis_config.json override built-in defaults."""
        config = {"vision": {"gemini_model": "gemini-2.0-flash"}}
        config_path = tmp_path / "analysis_config.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=config_path,
            matlab_config_path=tmp_path / "nonexistent.json",
        )
        assert cfg["vision"]["gemini_model"] == "gemini-2.0-flash"
        # Other vision keys should still have their defaults.
        assert cfg["vision"]["max_retries"] == 4

    def test_partial_statistics_override(self, tmp_path: Path):
        """Overriding one statistics key preserves the others."""
        config = {"statistics": {"p_noteworthy": 0.10}}
        config_path = tmp_path / "analysis_config.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=config_path,
            matlab_config_path=tmp_path / "nonexistent.json",
        )
        assert cfg["statistics"]["p_noteworthy"] == 0.10
        assert cfg["statistics"]["p_significant"] == 0.01

    def test_priority_graphs_override(self, tmp_path: Path):
        """Priority graphs list can be fully replaced."""
        config = {"priority_graphs": ["Custom_Graph_A", "Custom_Graph_B"]}
        config_path = tmp_path / "analysis_config.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=config_path,
            matlab_config_path=tmp_path / "nonexistent.json",
        )
        assert cfg["priority_graphs"] == ["Custom_Graph_A", "Custom_Graph_B"]

    def test_bad_json_falls_back_to_defaults(self, tmp_path: Path):
        """Malformed JSON in the config file is silently ignored."""
        config_path = tmp_path / "analysis_config.json"
        config_path.write_text("{bad json", encoding="utf-8")

        cfg = load_analysis_config(
            config_path=config_path,
            matlab_config_path=tmp_path / "nonexistent.json",
        )
        assert cfg == load_analysis_config(
            config_path=tmp_path / "nonexistent.json",
            matlab_config_path=tmp_path / "nonexistent.json",
        )

    def test_matlab_config_unknown_dwi_type_added(self, tmp_path: Path):
        """A novel dwi_type from MATLAB config is appended to the list."""
        matlab_cfg = {"dwi_type": "CustomDWI"}
        matlab_path = tmp_path / "config.json"
        matlab_path.write_text(json.dumps(matlab_cfg), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=tmp_path / "nonexistent.json",
            matlab_config_path=matlab_path,
        )
        assert "CustomDWI" in cfg["dwi_types"]
        # Original types should still be present.
        assert "Standard" in cfg["dwi_types"]

    def test_matlab_config_known_dwi_type_not_duplicated(self, tmp_path: Path):
        """A known dwi_type from MATLAB config does not duplicate the list."""
        matlab_cfg = {"dwi_type": "Standard"}
        matlab_path = tmp_path / "config.json"
        matlab_path.write_text(json.dumps(matlab_cfg), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=tmp_path / "nonexistent.json",
            matlab_config_path=matlab_path,
        )
        assert cfg["dwi_types"].count("Standard") == 1

    def test_matlab_config_missing_dwi_type_ignored(self, tmp_path: Path):
        """MATLAB config without dwi_type field is silently ignored."""
        matlab_cfg = {"some_other_field": "value"}
        matlab_path = tmp_path / "config.json"
        matlab_path.write_text(json.dumps(matlab_cfg), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=tmp_path / "nonexistent.json",
            matlab_config_path=matlab_path,
        )
        assert cfg["dwi_types"] == ["Standard", "dnCNN", "IVIMnet"]

    def test_full_override_chain(self, tmp_path: Path):
        """Both analysis config and MATLAB config apply in order."""
        analysis_cfg = {
            "vision": {"gemini_model": "custom-model"},
            "statistics": {"p_noteworthy": 0.10},
        }
        analysis_path = tmp_path / "analysis_config.json"
        analysis_path.write_text(json.dumps(analysis_cfg), encoding="utf-8")

        matlab_cfg = {"dwi_type": "ExperimentalDWI"}
        matlab_path = tmp_path / "config.json"
        matlab_path.write_text(json.dumps(matlab_cfg), encoding="utf-8")

        cfg = load_analysis_config(
            config_path=analysis_path,
            matlab_config_path=matlab_path,
        )
        assert cfg["vision"]["gemini_model"] == "custom-model"
        assert cfg["statistics"]["p_noteworthy"] == 0.10
        assert "ExperimentalDWI" in cfg["dwi_types"]


# ---------------------------------------------------------------------------
# get_config / reset_config_cache
# ---------------------------------------------------------------------------

class TestGetConfig:
    """Verify config caching and reset behaviour."""

    def test_returns_dict_with_expected_keys(self):
        """get_config returns a dict with all top-level keys from defaults."""
        reset_config_cache()
        cfg = get_config()
        assert "dwi_types" in cfg
        assert "vision" in cfg
        assert "statistics" in cfg
        assert "priority_graphs" in cfg

    def test_caching(self):
        """Repeated calls return the same object (identity check)."""
        reset_config_cache()
        cfg1 = get_config()
        cfg2 = get_config()
        assert cfg1 is cfg2

    def test_reset_clears_cache(self):
        """After reset, get_config reloads (new object)."""
        reset_config_cache()
        cfg1 = get_config()
        reset_config_cache()
        cfg2 = get_config()
        # Not the same object, but equal in value.
        assert cfg1 is not cfg2
        assert cfg1 == cfg2


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """Ensure the module-level DWI_TYPES constant still works."""

    def test_dwi_types_constant(self):
        """The legacy DWI_TYPES module constant matches the default config."""
        assert DWI_TYPES == ["Standard", "dnCNN", "IVIMnet"]

    def test_defaults_have_all_expected_sections(self):
        """_DEFAULTS contains all sections that scripts rely on."""
        assert "dwi_types" in _DEFAULTS
        assert "vision" in _DEFAULTS
        assert "statistics" in _DEFAULTS
        assert "priority_graphs" in _DEFAULTS
        assert "gemini_model" in _DEFAULTS["vision"]
        assert "p_noteworthy" in _DEFAULTS["statistics"]
