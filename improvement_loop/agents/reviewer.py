# reviewer.py
"""Reviewer agent — LLM-based code review quality gate for proposed fixes."""

import difflib
import json
from dataclasses import dataclass, field
from typing import List, Literal, Optional

from improvement_loop.agents._api import api_call_with_retry
from improvement_loop.evaluator import Finding
from improvement_loop.loop_config import get_config as _get_loop_config

REVIEW_SYSTEM_PROMPT = """\
You are a senior code reviewer for a MATLAB + Python medical-imaging research pipeline \
(pancreatic DWI analysis at Memorial Sloan Kettering).

You are reviewing a proposed code change (shown as a unified diff) that was generated to \
address a specific finding from a code audit. Evaluate the change against these criteria:

## Review criteria

1. **Relevance** — Does the diff actually address the finding description? Reject changes \
that are unrelated to the stated problem.

2. **Correctness** — Does the change introduce obvious bugs, regressions, or logic errors? \
Check for off-by-one errors, broken control flow, missing edge cases, incorrect MATLAB or \
Python semantics.

3. **Forbidden modifications** — Does the change modify anything under `pipeline/dependencies/`? \
This directory is READ-ONLY. Any modification here is an automatic reject.

4. **Data leakage safety** — Could the change introduce data leakage? Examples: removing \
patient-stratified cross-validation (`make_grouped_folds`), allowing cross-timepoint \
imputation bleed, removing train/test separation in `knn_impute_train_test`.

5. **PHI safety** — Could the change expose patient health information? Check for logging \
of patient IDs, writing patient data to non-secure locations, or removing PHI guards.

6. **Scope** — Is the change minimal and focused? Reject changes that include unrelated \
refactoring, gratuitous style changes, or feature additions beyond the stated fix.

7. **Config compatibility** — Does the change preserve backward compatibility of \
`config.json` fields? New fields must have defaults in `parse_config.m`; removed fields \
must not break existing configs.

## Risk flags
Include any applicable flags from this list:
- "LEAKAGE_RISK" — change could introduce data leakage
- "PHI_RISK" — change could expose patient data
- "DEPS_MODIFIED" — change modifies pipeline/dependencies/
- "INCORRECT_FINDING" — the original finding is technically wrong and the fix is misguided

## Output format
Return ONLY a JSON object with exactly these keys:
{
  "verdict": "approve" | "request_changes" | "reject",
  "reasoning": "1-3 sentence explanation",
  "risk_flags": ["FLAG1", ...]
}

Use "approve" when the change is correct, focused, and safe.
Use "request_changes" when the change has minor issues that could be fixed.
Use "reject" when the change is fundamentally flawed, dangerous, or out of scope.
Return an empty list for risk_flags if none apply.
"""

# Flags that force rejection regardless of the LLM's verdict
_CRITICAL_FLAGS = frozenset({"LEAKAGE_RISK", "PHI_RISK"})


@dataclass
class ReviewVerdict:
    """Structured result from a code review."""
    verdict: Literal["approve", "request_changes", "reject"]
    reasoning: str
    risk_flags: List[str]


def _generate_diff(original_content: str, new_content: str, filename: str = "file") -> str:
    """Generate a unified diff between original and new content."""
    original_lines = original_content.splitlines(keepends=True)
    new_lines = new_content.splitlines(keepends=True)
    diff_lines = difflib.unified_diff(
        original_lines, new_lines,
        fromfile=f"a/{filename}", tofile=f"b/{filename}",
    )
    return "".join(diff_lines)


def _parse_review_verdict(raw_text: str) -> Optional[ReviewVerdict]:
    """Parse a JSON review verdict from raw API response text.

    Returns None if parsing fails.
    """
    cleaned = raw_text.strip()
    if "```" in cleaned:
        parts = cleaned.split("```")
        if len(parts) >= 3:
            inner = parts[1]
            if inner.startswith("json"):
                inner = inner[4:]
            cleaned = inner.strip()
        elif cleaned.startswith("```"):
            cleaned = parts[1]
            if cleaned.startswith("json"):
                cleaned = cleaned[4:]
            cleaned = cleaned.strip()

    # Try to find a JSON object if text doesn't start with {
    if not cleaned.startswith("{"):
        brace_pos = cleaned.find("{")
        if brace_pos >= 0:
            cleaned = cleaned[brace_pos:]

    try:
        data = json.loads(cleaned)
    except json.JSONDecodeError:
        return None

    if not isinstance(data, dict):
        return None

    verdict = data.get("verdict")
    reasoning = data.get("reasoning")
    risk_flags = data.get("risk_flags")

    if verdict not in ("approve", "request_changes", "reject"):
        return None
    if not isinstance(reasoning, str):
        return None
    if not isinstance(risk_flags, list):
        return None

    return ReviewVerdict(
        verdict=verdict,
        reasoning=reasoning,
        risk_flags=[str(f) for f in risk_flags],
    )


def review(finding: Finding, original_content: str, new_content: str,
           dry_run: bool = False) -> ReviewVerdict:
    """Review a proposed code change and return a structured verdict."""
    if dry_run:
        return ReviewVerdict(verdict="approve", reasoning="dry-run", risk_flags=[])

    diff_text = _generate_diff(original_content, new_content, filename=finding.file)

    cfg = _get_loop_config()
    user_message = (
        f"## Finding being addressed\n"
        f"**Dimension:** {finding.dimension}\n"
        f"**File:** {finding.file}\n"
        f"**Description:** {finding.description}\n"
        f"**Proposed fix:** {finding.fix}\n\n"
        f"## Unified diff\n```diff\n{diff_text}\n```\n\n"
        f"Review this change and return your verdict as JSON."
    )

    raw_response = api_call_with_retry({
        "model": cfg.review_model,
        "max_tokens": cfg.review_max_tokens,
        "system": REVIEW_SYSTEM_PROMPT,
        "messages": [{"role": "user", "content": user_message}],
    })

    verdict = _parse_review_verdict(raw_response)
    if verdict is None:
        return ReviewVerdict(
            verdict="request_changes",
            reasoning="Review parse failed",
            risk_flags=["EVALUATION_FAILED"],
        )

    # Enforce critical flag override: LEAKAGE_RISK or PHI_RISK → always reject
    if _CRITICAL_FLAGS & set(verdict.risk_flags):
        triggered = _CRITICAL_FLAGS & set(verdict.risk_flags)
        if verdict.verdict != "reject":
            verdict = ReviewVerdict(
                verdict="reject",
                reasoning=(
                    f"Auto-rejected due to critical risk flags: {sorted(triggered)}. "
                    f"Original reasoning: {verdict.reasoning}"
                ),
                risk_flags=verdict.risk_flags,
            )

    return verdict
