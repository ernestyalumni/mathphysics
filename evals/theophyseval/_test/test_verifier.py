"""Known-answer sanity tests for verifier.py. Run: python -m pytest evals/theophyseval/_test"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from verifier import extract_numbers, final_region, grade  # noqa: E402


def test_number_forms():
    text = r"a = 1.5e-3, b = 2.0 \times 10^{-4}, c = 3×10^5, d = −7, e = 19/73, f = \pi, g = 2\pi/3"
    got = extract_numbers(text)
    want = [1.5e-3, 2.0e-4, 3e5, -7.0, 19 / 73, math.pi, 2 * math.pi / 3]
    assert len(got) == len(want)
    for x, y in zip(got, want):
        assert math.isclose(x, y, rel_tol=1e-12)


def test_thousands_separator_and_latex_spacing():
    assert extract_numbers(r"E = -9.337\,233 and N = 48,620") == [-9.337233, 48620.0]


def test_boxed_region_with_nested_braces():
    text = r"intermediate 0.5 ... \boxed{C=3,\ \Delta_{\mathrm{dir}}=0.01228,\ \Delta=0.01196}"
    assert final_region(text) .startswith("C=3")
    assert 0.5 not in extract_numbers(final_region(text))


def test_keyed_region():
    text = "first try ZEROES=[1,2,3]\nfinal ZEROES=[0.1234, 0.5678, 0.9012]"
    assert extract_numbers(final_region(text, "ZEROES")) == [0.1234, 0.5678, 0.9012]


def test_grade_pass_fail_inconclusive():
    ref = {"answers": [{"name": "C", "value": 3, "atol": 0.01},
                       {"name": "gap", "value": 0.0122790, "atol": 2e-5}]}
    assert grade(r"\boxed{C=3, gap=0.01229}", ref)["verifier_status"] == "PASS"
    bad = grade(r"\boxed{C=-3, gap=0.01229}", ref)
    assert bad["verifier_status"] == "FAIL" and bad["confident_wrong"]
    assert grade("no numbers here", ref)["verifier_status"] == "INCONCLUSIVE"


def test_admitted_unknown_is_not_confident_wrong():
    ref = {"answers": [{"name": "x", "value": 1.0, "atol": 1e-3}]}
    row = grade("I could not converge; the value cannot be determined. Best guess 2.0", ref)
    assert row["verifier_status"] == "FAIL"
    assert row["admitted_unknown"] and not row["confident_wrong"]


def test_ordered_matching():
    ref = {"answers": [{"name": "a", "value": 1.0, "atol": 1e-9},
                       {"name": "b", "value": 2.0, "atol": 1e-9}], "ordered": True}
    assert grade(r"\boxed{1, 2}", ref)["verifier_status"] == "PASS"
    assert grade(r"\boxed{2, 1}", ref)["verifier_status"] == "FAIL"


def test_alternatives():
    ref = {"answers": [{"name": "T", "value": 0, "alternatives": [1], "atol": 1e-6}]}
    assert grade(r"\boxed{T=1}", ref)["verifier_status"] == "PASS"
    assert grade(r"\boxed{k=0}", ref)["verifier_status"] == "PASS"
    assert grade(r"\boxed{T=-1}", ref)["verifier_status"] == "FAIL"
