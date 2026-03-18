# evaluator.py
import anthropic  # type: ignore
import json
import re
import time
import os
from typing import List, Literal, Optional

from pydantic import BaseModel, field_validator


VALID_DIMENSIONS = (
    "performance", "correctness", "error_handling", "modularity",
    "memory", "code_quality", "test_coverage", "security", "cross_platform",
)

VALID_STATUSES = ("pending", "implemented", "merged")

# Characters forbidden in git branch names (see git-check-ref-format)
_GIT_BRANCH_INVALID_RE = re.compile(
    r"[\x00-\x1f\x7f ~^:?*\[\\]"  # control chars + special chars
    r"|\.\.+"                        # consecutive dots
    r"|@\{"                          # @{ sequence
    r"|\.$"                          # trailing dot
    r"|\.lock$"                      # .lock suffix
)


class Finding(BaseModel):
    """Typed schema for a single codebase improvement finding."""

    dimension: Literal[
        "performance", "correctness", "error_handling", "modularity",
        "memory", "code_quality", "test_coverage", "security", "cross_platform",
    ]
    file: str
    function_name: Optional[str] = None
    description: str
    fix: str
    importance: int
    branch_name: str
    status: Optional[Literal["pending", "implemented", "merged"]] = None

    @field_validator("importance")
    @classmethod
    def _importance_range(cls, v: int) -> int:
        if not (1 <= v <= 10):
            raise ValueError("importance must be between 1 and 10 inclusive")
        return v

    @field_validator("branch_name")
    @classmethod
    def _valid_branch_name(cls, v: str) -> str:
        if not v.startswith("improvement/"):
            raise ValueError("branch_name must start with 'improvement/'")
        slug = v[len("improvement/"):]
        if not slug:
            raise ValueError("branch_name slug must not be empty")
        if len(slug) > 50:
            raise ValueError("branch_name slug must be at most 50 characters")
        if " " in v:
            raise ValueError("branch_name must not contain spaces")
        # No extra slashes beyond the leading "improvement/"
        if "/" in slug:
            raise ValueError(
                "branch_name must not contain slashes after 'improvement/'"
            )
        if _GIT_BRANCH_INVALID_RE.search(v):
            raise ValueError(
                "branch_name contains characters invalid in git branch names"
            )
        return v

    def to_log_dict(self) -> dict:
        """Serialize to the dict format used by loop_tracker's JSON log."""
        d: dict = {
            "dimension": self.dimension,
            "file": self.file,
            "description": self.description,
            "fix": self.fix,
            "importance": self.importance,
            "branch_name": self.branch_name,
        }
        if self.function_name is not None:
            d["function_name"] = self.function_name
        if self.status is not None:
            d["status"] = self.status
        return d

client = anthropic.Anthropic()

JUDGE_SYSTEM_PROMPT = """\
You are an expert code-audit judge for a MATLAB + Python medical-imaging research pipeline \
(pancreatic DWI analysis at Memorial Sloan Kettering). Your job is to evaluate the quality \
of a codebase improvement audit and return a structured JSON score.

## Domain context
- The pipeline fits IVIM and ADC diffusion models to MRI data, performs survival analysis, \
competing risks modeling, and treatment response prediction.
- It uses MATLAB (R2021a+) with parfor parallelism, checkpointing, and cross-platform \
(Windows/macOS/Linux) support. Post-hoc analysis is in Python 3.12+.
- Key concerns: data leakage prevention, patient-stratified CV, temporal leakage in \
imputation, file deletion safety (sentinel-based provenance), shell injection prevention, \
and PHI protection.

## Scoring dimensions (each 0-10)

**specificity** — Are findings concrete and actionable? Do they cite specific files, \
functions, and line numbers? A vague "improve error handling" scores low; \
"add try-catch around the dcm2niix system() call in convert_dicom.m:47 to handle \
missing executables" scores high.

**accuracy** — Are the findings technically correct? Would implementing them actually \
improve the code? Penalise findings that misunderstand MATLAB semantics, confuse Octave \
shims with bugs, suggest changes to pipeline/dependencies/ (read-only), or recommend \
patterns that would introduce data leakage.

**coverage** — Does the audit examine all relevant dimensions: performance, correctness, \
error handling, modularity, memory usage, code quality, test coverage, security, and \
cross-platform compatibility? An audit that only looks at formatting scores low.

**prioritization** — Are importance scores well-calibrated? Critical bugs and data-integrity \
issues should score 8-10; minor style nits should score 1-2. Penalise audits that inflate \
trivial findings or underrate serious ones.

**domain_appropriateness** — Do the findings respect the medical-physics domain? \
Suggestions must not compromise: patient data safety, statistical rigour (leakage guards, \
FDR correction), clinical interpretability, or backwards compatibility of config.json. \
Bonus for findings that improve reproducibility or audit trails.

**overall** — Holistic quality of the audit as a guide for the next improvement iteration. \
Weight accuracy and domain_appropriateness most heavily.

## Flags
Return a list of string flags for any of these conditions:
- "LEAKAGE_RISK" — a suggested change could introduce data leakage
- "PHI_RISK" — a suggestion could expose patient data
- "DEPS_MODIFIED" — audit suggests modifying pipeline/dependencies/
- "INFLATED_SCORES" — importance scores are systematically too high
- "DEFLATED_SCORES" — importance scores are systematically too low
- "MISSING_DIMENSION" — an entire audit dimension was skipped
- "INCORRECT_FINDING" — one or more findings are technically wrong

Return an empty list if no flags apply.

## Output format
Return ONLY a JSON object with exactly these keys:
{
  "specificity": <float 0-10>,
  "accuracy": <float 0-10>,
  "coverage": <float 0-10>,
  "prioritization": <float 0-10>,
  "domain_appropriateness": <float 0-10>,
  "overall": <float 0-10>,
  "flags": [<string>, ...],
  "reasoning": "<1-3 sentence justification for the overall score>"
}
"""

CALIBRATION_EXAMPLES = """\
## Calibration examples

### Example 1 — High-quality audit (expected overall ~8.5)
Audit excerpt:
"Finding 1 (importance 8/10): fit_models.m:132 — the IVIM bi-exponential fit uses lsqcurvefit \
with default TolFun=1e-6, but for low-SNR voxels this causes non-convergence. Setting \
TolFun=1e-4 and adding a convergence check reduces failed fits by ~12% on synthetic data. \
Finding 2 (importance 3/10): text_progress_bar.m:28 — the bar string is rebuilt every \
iteration via string concatenation; pre-allocating a char array would save ~0.2s for 500 \
patients but is cosmetic."
Expected scores: specificity=9, accuracy=8, coverage=7, prioritization=9, \
domain_appropriateness=8, overall=8.5, flags=[], \
reasoning="Specific file:line citations, correct MATLAB semantics, well-calibrated importance."

### Example 2 — Medium-quality audit (expected overall ~5)
Audit excerpt:
"Finding 1 (importance 7/10): The pipeline should use more vectorization for speed. \
Finding 2 (importance 6/10): Error messages could be more descriptive. \
Finding 3 (importance 5/10): Consider adding type hints to all functions."
Expected scores: specificity=2, accuracy=5, coverage=4, prioritization=4, \
domain_appropriateness=5, overall=5, flags=["MISSING_DIMENSION"], \
reasoning="Findings are vague with no file references; MATLAB does not have type hints; \
security and leakage dimensions not examined."

### Example 3 — Poor audit with risky suggestions (expected overall ~2)
Audit excerpt:
"Finding 1 (importance 9/10): Replace make_grouped_folds.m with simple random k-fold CV \
for cleaner code. Finding 2 (importance 8/10): Cache patient features in a global variable \
to avoid recomputation across pipeline steps. Finding 3 (importance 7/10): Store intermediate \
results in pipeline/dependencies/cache/ for faster reruns."
Expected scores: specificity=3, accuracy=1, coverage=2, prioritization=1, \
domain_appropriateness=0, overall=2, flags=["LEAKAGE_RISK", "DEPS_MODIFIED", "INFLATED_SCORES"], \
reasoning="Suggests removing patient-stratified CV (leakage), global state (violates architecture), \
and writing to dependencies/ (read-only). Importance scores are wildly inflated for harmful changes."
"""

FULL_SYSTEM_PROMPT = JUDGE_SYSTEM_PROMPT + "\n\n" + CALIBRATION_EXAMPLES

REQUIRED_KEYS = {
    "specificity", "accuracy", "coverage",
    "prioritization", "domain_appropriateness",
    "overall", "flags", "reasoning"
}

SCORE_KEYS = {
    "specificity", "accuracy", "coverage",
    "prioritization", "domain_appropriateness", "overall"
}

MAX_RETRIES = 3
RETRY_DELAY = 2.0


def parse_and_validate(raw_text: str) -> Optional[dict]:
    cleaned = raw_text.strip()
    if cleaned.startswith("```"):
        cleaned = cleaned.split("```")[1]
        if cleaned.startswith("json"):
            cleaned = cleaned[4:]  # type: ignore
    cleaned = cleaned.strip()
    
    try:
        data = json.loads(cleaned)
    except json.JSONDecodeError as e:
        print(f"JSON parse failed: {e}")
        return None
    
    missing = REQUIRED_KEYS - set(data.keys())
    if missing:
        print(f"Missing keys: {missing}")
        return None
    
    for key in SCORE_KEYS:
        val = data[key]
        if not isinstance(val, (int, float)) or not (0 <= val <= 10):
            print(f"Invalid score for '{key}': {val}")
            return None
    
    return data


def score_audit(audit_output: str) -> dict:
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=1000,
                system=FULL_SYSTEM_PROMPT,
                messages=[
                    {
                        "role": "user",
                        "content": f"Score this audit:\n\n{audit_output}"
                    }
                ]
            )
            result = parse_and_validate(response.content[0].text)
            if result is not None:
                return result
        except anthropic.APIError as e:
            print(f"Attempt {attempt} API error: {e}")
        
        time.sleep(RETRY_DELAY * attempt)
    
    return {
        "specificity": 5.0, "accuracy": 5.0, "coverage": 5.0,
        "prioritization": 5.0, "domain_appropriateness": 5.0,
        "overall": 5.0,
        "flags": ["EVALUATION_FAILED"],
        "reasoning": "Fallback scores — evaluator failed after max retries."
    }


def should_continue_loop(scores: dict, findings: List[Finding]) -> bool:
    """
    Returns True if the loop should continue, False if it should exit.
    Combines test results with judge scores.
    """
    # Any finding above threshold — keep going
    high_priority = [f for f in findings if f.importance >= 2]
    if high_priority:
        print(f"Continuing — {len(high_priority)} findings at importance >= 2")
        return True
    
    # Judge thinks coverage is poor — keep going
    if scores["coverage"] < 6:
        print(f"Continuing — audit coverage score {scores['coverage']}/10 is too low")
        return True
    
    # Flag critical issues
    if scores["flags"] and "EVALUATION_FAILED" not in scores["flags"]:
        print(f"Continuing — judge flagged issues: {scores['flags']}")
        return True
    
    print("Exit condition met — no findings above threshold and audit quality sufficient")
    return False