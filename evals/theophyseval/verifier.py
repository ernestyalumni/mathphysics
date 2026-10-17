"""Numeric verifier for TheoPhysEval-style answers.

The paper grades only numbers: every symbolic answer was converted to a number by fixing the
physical parameters, which removes the need for an LLM judge. This module does the same.

Reference format (one entry per problem id):

    {"answers": [{"name": "C_occ", "value": 3, "atol": 0.01},
                 {"name": "T", "value": 0, "alternatives": [1], "atol": 1e-6},
                 {"name": "direct_gap", "value": 0.0123, "atol": 2e-5}],
     "key": "ZEROES",          # optional: only read numbers after 'KEY=' in the response
     "ordered": false,         # optional: require the values in this order
     "source": "...", "official": false}

grade(response, reference) returns a row with verifier_status PASS | FAIL | INCONCLUSIVE
and admitted_unknown, following agent/EVAL-AND-TOOLING.md.
"""

from __future__ import annotations

import math
import re

_NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
_TIMES_TEN = re.compile(
    rf"({_NUM})\s*(?:\\times|×|\*|x)\s*10\s*\^\s*\{{?\s*([-+]?\d+)\s*\}}?")
_FRACTION = re.compile(rf"(?<![\w.])({_NUM})\s*/\s*({_NUM})(?![\w.])")
_PI = re.compile(rf"({_NUM})?\s*(?:\\pi|π)(?:\s*/\s*({_NUM}))?")
_PLAIN = re.compile(_NUM)
_UNKNOWN = re.compile(
    r"(cannot be determined|can ?not determine|unable to (?:solve|determine|compute)|"
    r"I don't know|I do not know|no such (?:equation|solution) exists|insufficient information|"
    r"impossible to (?:solve|determine))", re.IGNORECASE)


def normalize(text: str) -> str:
    """Unify minus signs and thin spaces, drop LaTeX spacing and thousands separators."""
    text = text.replace("\u2212", "-").replace("\u2009", " ").replace("\\,", "").replace("\\!", "")
    text = text.replace("\\left", "").replace("\\right", "")
    return re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)


def extract_numbers(text: str) -> list[float]:
    """All numbers in reading order. Handles 1.2e-3, 1.2 \\times 10^{-3}, 1.2×10^-3, a/b and
    multiples of pi. Each character span is consumed once."""
    text = normalize(text)
    spans: list[tuple[int, int, float]] = []

    def take(regex, conv):
        for m in regex.finditer(text):
            if any(s < m.end() and m.start() < e for s, e, _ in spans):
                continue
            try:
                spans.append((m.start(), m.end(), conv(m)))
            except (ValueError, ZeroDivisionError):
                pass

    take(_TIMES_TEN, lambda m: float(m.group(1)) * 10.0 ** int(m.group(2)))
    take(_PI, lambda m: (float(m.group(1)) if m.group(1) else 1.0) * math.pi
         / (float(m.group(2)) if m.group(2) else 1.0))
    take(_FRACTION, lambda m: float(m.group(1)) / float(m.group(2)))
    take(_PLAIN, lambda m: float(m.group(0)))
    return [v for _, _, v in sorted(spans)]


def final_region(text: str, key: str | None = None, tail_lines: int = 12) -> str:
    """Where the answer is expected: after the last KEY=, else the last \\boxed{...}, else
    the last tail_lines lines. Grading only this region keeps intermediate numbers out."""
    if key:
        idx = text.rfind(key)
        return text[idx + len(key):] if idx >= 0 else ""
    boxed = [m.start() for m in re.finditer(r"\\boxed\{", text)]
    if boxed:
        start = boxed[-1] + len("\\boxed{")
        depth, i = 1, start
        while i < len(text) and depth:
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            i += 1
        return text[start:i - 1]
    return "\n".join(text.strip().splitlines()[-tail_lines:])


def _close(x: float, want: dict) -> bool:
    """Within tolerance of value or of any entry in want["alternatives"] (equivalent forms,
    e.g. T = e^{ik} reported as 1 or as k = 0, or a Berry phase of pi reported as -pi)."""
    atol = want.get("atol", 0.0)
    rtol = want.get("rtol", 0.0 if "atol" in want else 1e-6)
    return any(abs(x - v) <= atol + rtol * abs(v) for v in [want["value"], *want.get("alternatives", [])])


def match(found: list[float], answers: list[dict], ordered: bool) -> list[bool]:
    """For each wanted value, is there a found number within tolerance? Each found number
    is used at most once; in ordered mode the matches must appear left to right."""
    used, hits, cursor = set(), [], 0
    for want in answers:
        rng = range(cursor, len(found)) if ordered else range(len(found))
        k = next((i for i in rng if i not in used and _close(found[i], want)), None)
        hits.append(k is not None)
        if k is not None:
            used.add(k)
            cursor = k + 1
    return hits


def grade(response: str, reference: dict) -> dict:
    region = final_region(response, reference.get("key"))
    found = extract_numbers(region)
    admitted = bool(_UNKNOWN.search(response))
    if not found:
        status, hits = "INCONCLUSIVE", []
    else:
        hits = match(found, reference["answers"], reference.get("ordered", False))
        status = "PASS" if all(hits) else "FAIL"
    return {
        "verifier_status": status,
        "admitted_unknown": admitted,
        "confident_wrong": status == "FAIL" and not admitted,
        "matched": {a["name"]: h for a, h in zip(reference["answers"], hits)},
        "numbers_in_final_region": found[:40],
    }
