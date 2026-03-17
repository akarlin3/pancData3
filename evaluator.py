# evaluator.py
import anthropic  # type: ignore
import json
import time
import os
from typing import Optional

client = anthropic.Anthropic()

JUDGE_SYSTEM_PROMPT = """..."""  # Full prompt from section 1
CALIBRATION_EXAMPLES = """..."""  # Examples from section 4

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


def should_continue_loop(scores: dict, findings: list) -> bool:
    """
    Returns True if the loop should continue, False if it should exit.
    Combines test results with judge scores.
    """
    # Any finding above threshold — keep going
    high_priority = [f for f in findings if f.get("importance", 0) >= 2]
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