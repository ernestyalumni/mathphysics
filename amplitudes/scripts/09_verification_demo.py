"""
09_verification_demo -- Two-minute proposer / independent-verifier demo.

Orchestrates the two existing, independently implemented code paths in
mhvamplitudes (closed-form SMGA formula vs. Berends-Giele recursion) to show
the core product behavior: candidate -> independent verifier -> PASS / FAIL
ledger, including one deliberately wrong candidate so the demo shows the
verifier actually rejecting something, not just displaying green checks.

Does not reimplement any physics: it imports and composes the existing
mhvamplitudes.amplitudes / mhvamplitudes.kinematics modules. The one
"broken" candidate below reuses the real sg_ij / sg_i_set primitives and
introduces a single deliberate bug in how they're composed (a wrong
normalization exponent) -- exactly the class of subtle, plausible-looking
error this product exists to catch.

Reference: arXiv:2602.12176 (Guevara, Lupsasca, Skinner, Strominger, Weil).

Truth label: this reproduces a known result from a published paper through
n = 10. It is a verification/reproducibility milestone, not a new physics
result and not evidence beyond the paper's own all-n formula and proof.

Usage (from repo root):
    uv run --project amplitudes/MHVamplitudes python amplitudes/scripts/09_verification_demo.py
"""

import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "MHVamplitudes" / "src"))

from mhvamplitudes.amplitudes.smga import sg_ij, sg_i_set, smga_formula_r1
from mhvamplitudes.amplitudes.berends_giele import compute_stripped_amplitude
from mhvamplitudes.kinematics.phase_space import make_r1_config

TOL = 1e-10
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = REPO_ROOT / "amplitudes" / "results" / "demo-latest.md"


def broken_smga_formula_r1(n, omegas, tilde_zs):
    """Deliberately wrong candidate: normalizes by 2**(n-1) instead of the
    correct 2**(n-2). Same sg_ij / sg_i_set primitives as the real formula."""
    result = 1.0
    for m in range(2, n):
        i1, i2 = m - 1, m
        sg_mm1 = sg_ij(omegas[i1], tilde_zs[i1], omegas[i2], tilde_zs[i2])
        sg_1_S = sg_i_set(omegas[0], tilde_zs[0], omegas[1:m], tilde_zs[1:m])
        result *= (sg_mm1 + sg_1_S)
    return result / (2 ** (n - 1))  # BUG: should be 2 ** (n - 2)


def git_info():
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True
    ).stdout.strip()
    branch = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True
    ).stdout.strip()
    return sha, branch


def run_comparison(n, seed, num_trials, formula_fn):
    rng = np.random.default_rng(seed)
    max_residual = 0.0
    status = "PASS"
    for _ in range(num_trials):
        omegas, tilde_zs = make_r1_config(n, rng)
        A_candidate = formula_fn(n, omegas, tilde_zs)
        A_bg = compute_stripped_amplitude(n, omegas, tilde_zs, minus_particle=0)
        residual = abs(A_candidate - A_bg)
        max_residual = max(max_residual, residual)
        if residual >= TOL:
            status = "FAIL"
    return status, max_residual


def main():
    sha, branch = git_info()
    now = datetime.now(timezone.utc).isoformat()

    print("=" * 70)
    print("Verification demo -- proposer -> independent verifier -> ledger")
    print("=" * 70)
    print()
    print("Candidate formula: SMGA closed form (arXiv:2602.12176, eq. 16)")
    print("  A_{1...n}|_R1 = (1/2^(n-2)) * prod_{m=2}^{n-1} (sg_{m,m+1} + sg_{1,2...m})")
    print("Independent verifier: Berends-Giele recursion (same paper, eq. 8 / eq. A17),")
    print("  a separate code path: mhvamplitudes.amplitudes.berends_giele")
    print()

    rows = []

    print("--- Known-correct candidate: closed form vs. recursion, n = 3..10 ---")
    for n in range(3, 11):
        num_trials = 10 if n <= 7 else 3
        seed = 1000 + n
        status, max_residual = run_comparison(n, seed, num_trials, smga_formula_r1)
        print(f"  n={n:2d}  trials={num_trials:2d}  seed={seed}  max_residual={max_residual:.2e}  {status}")
        rows.append({"candidate": "correct closed form", "n": n, "trials": num_trials,
                     "seed": seed, "max_residual": max_residual, "status": status})

    print()
    print("--- Deliberately wrong candidate: off-by-one normalization exponent ---")
    # Restricted to small n: as n grows, the product has more sign factors and the
    # true amplitude is increasingly likely to be exactly 0 in a given random trial,
    # which would spuriously "pass" the broken candidate too (0/2 == 0). Small n keeps
    # the demo's false-negative rate negligible without touching any physics code.
    for n in (4, 5, 6):
        num_trials = 10
        seed = 2000 + n
        status, max_residual = run_comparison(n, seed, num_trials, broken_smga_formula_r1)
        print(f"  n={n:2d}  trials={num_trials:2d}  seed={seed}  max_residual={max_residual:.2e}  "
              f"{status}  (expected FAIL)")
        rows.append({"candidate": "broken closed form (2^(n-1) bug)", "n": n, "trials": num_trials,
                     "seed": seed, "max_residual": max_residual, "status": status})

    correct_rows = [r for r in rows if r["candidate"] == "correct closed form"]
    broken_rows = [r for r in rows if r["candidate"].startswith("broken")]
    overall_pass = all(r["status"] == "PASS" for r in correct_rows)
    caught_the_bug = all(r["status"] == "FAIL" for r in broken_rows)

    print()
    print("=" * 70)
    print(f"Correct candidate: {'PASS on all n' if overall_pass else 'FAIL -- investigate'}")
    print(f"Broken candidate:  "
          f"{'correctly caught (FAIL) on all cases' if caught_the_bug else 'NOT caught -- verifier bug'}")
    print("=" * 70)
    print()
    print("LABEL: verified reproduction of a known result (arXiv:2602.12176); no novelty claim.")

    write_report(sha, branch, now, rows, overall_pass, caught_the_bug)
    print(f"\nReport written to {RESULTS_PATH.relative_to(REPO_ROOT)}")

    sys.exit(0 if (overall_pass and caught_the_bug) else 1)


def write_report(sha, branch, now, rows, overall_pass, caught_the_bug):
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Verification demo -- latest run",
        "",
        f"- Timestamp (UTC): {now}",
        f"- Git SHA: {sha}",
        f"- Branch: {branch}",
        f"- Backend: python {sys.version.split()[0]}, numpy {np.__version__}",
        "- Reference: arXiv:2602.12176 (Guevara, Lupsasca, Skinner, Strominger, Weil)",
        "",
        "## Pass/fail criteria",
        "",
        f"- Residual tolerance: {TOL:.0e}",
        "- A candidate PASSES a case if the max residual across trials is below tolerance.",
        "- A candidate FAILS a case if any trial's residual is at or above tolerance.",
        "",
        "## Results",
        "",
        "| Candidate | n | trials | seed | max residual | status |",
        "|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['candidate']} | {r['n']} | {r['trials']} | {r['seed']} | "
            f"{r['max_residual']:.2e} | {r['status']} |"
        )
    lines += [
        "",
        "## Summary",
        "",
        f"- Correct candidate: {'PASS on all n' if overall_pass else 'FAIL -- investigate'}",
        f"- Broken candidate (deliberate off-by-one exponent bug): "
        f"{'correctly caught (FAIL) on all cases' if caught_the_bug else 'NOT caught -- verifier bug'}",
        "",
        "**LABEL: verified reproduction of a known result; no novelty claim.**",
        "",
    ]
    RESULTS_PATH.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
