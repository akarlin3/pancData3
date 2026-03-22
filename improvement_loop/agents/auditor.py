# auditor.py
"""Auditor agent — runs code audits and returns parsed findings."""

import json
import os
import re
from typing import List

from improvement_loop.agents._api import api_call_with_retry
from improvement_loop.evaluator import Finding
from improvement_loop.loop_config import get_config as _get_loop_config

# Repo root is two levels up from this file's directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

AUDIT_SYSTEM_PROMPT = """\
You are a senior code auditor for a MATLAB + Python medical-imaging research pipeline \
(pancreatic DWI analysis at Memorial Sloan Kettering).

## Repository overview
- MATLAB pipeline in `pipeline/` — IVIM and ADC diffusion model fitting, survival analysis, \
competing risks, treatment response prediction.
- Python analysis in `analysis/` — post-hoc parsing, cross-DWI comparison, HTML/PDF reports.
- Config-driven via `config.json`; backwards-compatible defaults in `pipeline/utils/parse_config.m`.
- `pipeline/dependencies/` is READ-ONLY — never suggest modifications there.
- Patient data must NEVER be exposed; no PHI in code, logs, or suggestions.

## Your task
Audit the codebase and return a JSON array of findings. Each finding must have:
{
  "dimension": one of ("performance", "correctness", "error_handling", "modularity", \
"memory", "code_quality", "test_coverage", "security", "cross_platform"),
  "file": "path/to/file.m or .py",
  "function_name": "optional — specific function if applicable",
  "description": "what the problem is — cite specific lines",
  "fix": "concise fix description (1-3 sentences, NO full code blocks — just describe the change)",
  "importance": integer 1-10 (8-10 = critical bugs/data integrity, 4-7 = moderate, 1-3 = minor),
  "branch_name": "improvement/<slug>" (slug: lowercase, hyphens, no spaces, max 50 chars after prefix)
}

## Rules
- Be specific: cite files, functions, and line numbers.
- Do NOT suggest changes to `pipeline/dependencies/`.
- Do NOT suggest changes that would introduce data leakage (removing patient-stratified CV, \
cross-timepoint imputation, etc.).
- Do NOT inflate importance scores — style nits are 1-2, not 7-8.
- Cover multiple dimensions if possible (performance, correctness, security, tests, etc.).
- Return ONLY the JSON array — no markdown fences, no commentary.
"""


def collect_source_files() -> str:
    """Collect key source files as context for the audit."""
    cfg = _get_loop_config()
    # Smaller focused set to stay under rate limits
    key_files = [
        "pipeline/run_dwi_pipeline.m",
        "pipeline/core/fit_models.m",
        "pipeline/core/metrics_survival.m",
        "pipeline/core/metrics_stats_predictive.m",
        "pipeline/utils/parse_config.m",
        "pipeline/utils/knn_impute_train_test.m",
        "pipeline/utils/safe_load_mask.m",
        "pipeline/utils/escape_shell_arg.m",
        "analysis/run_analysis.py",
        "analysis/shared.py",
    ]

    parts = []
    for rel_path in key_files:
        full_path = os.path.join(REPO_ROOT, rel_path)
        if not os.path.exists(full_path):
            continue
        try:
            with open(full_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            if len(content) > cfg.max_file_chars:
                content = content[:cfg.max_file_chars] + "\n... [truncated]"
            parts.append(f"=== {rel_path} ===\n{content}")
        except OSError:
            continue
    return "\n\n".join(parts)


def _call_audit_api(iteration: int, context: str) -> str:
    """Call the Claude API to run a code audit. Returns raw audit text."""
    cfg = _get_loop_config()

    if cfg.rag_enabled:
        try:
            from improvement_loop.rag.retriever import get_context_for_audit
            source_context = get_context_for_audit(context)
        except Exception as e:
            print(f"⚠️  RAG context failed, falling back to static files: {e}")
            source_context = collect_source_files()
    else:
        source_context = collect_source_files()

    user_message = (
        f"## Iteration context\n{context}\n\n"
        f"## Source files\n{source_context}\n\n"
        "Return your findings as a JSON array."
    )

    return api_call_with_retry({
        "model": cfg.audit_model,
        "max_tokens": cfg.audit_max_tokens,
        "system": AUDIT_SYSTEM_PROMPT,
        "messages": [{"role": "user", "content": user_message}],
    })


def parse_findings(audit_output: str) -> List[Finding]:
    """Parse findings from audit output JSON."""
    # Strip markdown fences if present (handles preamble before the fence)
    text = audit_output.strip()
    if "```" in text:
        # Extract content between first pair of ``` fences
        parts = text.split("```")
        if len(parts) >= 3:
            # parts[1] is the content between the first pair of fences
            inner = parts[1]
            if inner.startswith("json"):
                inner = inner[4:]
            text = inner.strip()
        elif text.startswith("```"):
            text = parts[1]
            if text.startswith("json"):
                text = text[4:]
            text = text.strip()

    # If text still doesn't look like JSON, try to find a JSON array directly
    if not text.startswith("["):
        bracket_pos = text.find("[")
        if bracket_pos >= 0:
            text = text[bracket_pos:]



    # Truncation guard: check if the JSON array appears complete
    if not text.rstrip().endswith("]"):
        print("⚠️  Audit response appears truncated — consider increasing max_tokens further")
        # Attempt repair: find the last complete finding object by looking
        # for "},\n  {" or "}\n  {" boundaries (the split between findings).
        # Then try progressively shorter substrings.
        # Find all positions where a finding object boundary occurs
        boundaries = [m.end() - 1 for m in re.finditer(r'\}\s*,\s*\{', text)]
        # Also try the simple last-} approach as fallback
        last_brace = text.rfind("}")
        candidates = sorted(set(boundaries + ([last_brace] if last_brace >= 0 else [])),
                            reverse=True)
        for pos in candidates:
            # Take everything up to and including the } at this boundary
            # For boundary matches, pos points to the comma; use the } before it
            chunk = text[:pos + 1].rstrip().rstrip(",")
            if not chunk.lstrip().startswith("["):
                chunk = "[" + chunk
            chunk = chunk + "]"
            try:
                raw_list = json.loads(chunk)
                if not isinstance(raw_list, list):
                    continue
                print(f"⚠️  Recovered {len(raw_list)} findings from truncated response")
                findings = []
                for i, raw in enumerate(raw_list):
                    try:
                        finding = Finding(**raw)
                        findings.append(finding)
                    except Exception as e:
                        print(f"WARNING: Skipping finding {i}: {e}")
                return findings
            except json.JSONDecodeError:
                continue
        print("WARNING: Could not recover any findings from truncated response")
        print(f"Raw output (last 200 chars): {text[-200:]}")
        return []

    try:
        raw_list = json.loads(text)
    except json.JSONDecodeError as e:
        print(f"WARNING: Could not parse audit findings as JSON: {e}")
        print(f"Raw output (first 500 chars): {audit_output[:500]}")
        return []

    if not isinstance(raw_list, list):
        print(f"WARNING: Expected JSON array, got {type(raw_list).__name__}")
        return []

    findings = []
    for i, raw in enumerate(raw_list):
        try:
            finding = Finding(**raw)
            findings.append(finding)
        except Exception as e:
            print(f"WARNING: Skipping finding {i}: {e}")
    return findings


def audit(iteration: int, context: str, dry_run: bool = False) -> List[Finding]:
    """Run a code audit and return parsed findings."""
    if dry_run:
        return []

    raw_output = _call_audit_api(iteration, context)
    return parse_findings(raw_output)
