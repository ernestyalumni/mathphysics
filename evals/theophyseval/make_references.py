"""Write references.local.json from the solvers in solutions/.

These are self-computed reference answers, cross-checked inside each solver, NOT the official
keys (those are in the gated Hugging Face bundle). The file is deliberately untracked: the
benchmark authors keep answers off the public web to avoid training-data contamination, and
this repo has a public remote.

Run: python make_references.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "solutions"))

import maple_leaf as ml  # noqa: E402
import theta_logistic as tl  # noqa: E402


def dirac_reference():
    nodes, _ = ml.find_nodes(0.0)
    coords = [c for x, _ in nodes for c in x]
    energy = float(np.mean([e for _, e in nodes]))
    return {
        "answers": [{"name": f"node_coord_{i}", "value": float(c), "atol": 5e-6} for i, c in enumerate(coords)]
        + [{"name": "nodal_energy", "value": energy, "atol": 1e-8},
           {"name": "n_nodes", "value": len(nodes), "atol": 0.01},
           {"name": "symmetry_enforced", "value": 0, "atol": 0.01},
           {"name": "berry_phase_mod_2pi", "value": math.pi, "alternatives": [-math.pi], "atol": 1e-3}],
        "source": "solutions/maple_leaf.py dirac", "official": False,
    }


def chern_reference():
    ex = ml.band_extrema(0.15)
    return {
        "answers": [{"name": "C_occ", "value": round(ml.chern_occupied(0.15, 96)), "atol": 0.01},
                    {"name": "direct_gap", "value": float(ex["direct"]), "atol": 2e-5},
                    {"name": "indirect_gap", "value": float(ex["e4_min"] - ex["e3_max"]), "atol": 2e-5}],
        "source": "solutions/maple_leaf.py chern (link-variable sign convention of the prompt)",
        "official": False,
    }


def ed_reference():
    w0 = np.sort(eigsh(ml.heisenberg_sector(9)[0], k=3, which="SA", tol=1e-14)[0])
    w1 = np.sort(eigsh(ml.heisenberg_sector(10)[0], k=3, which="SA", tol=1e-14)[0])
    E0 = float(w0[0])
    return {
        "answers": [{"name": "E0", "value": E0, "atol": 1e-8},
                    {"name": "E0_per_site", "value": E0 / 18, "atol": 1e-8},
                    {"name": "S", "value": 0, "atol": 0.01},
                    {"name": "k_or_T", "value": 0, "alternatives": [1], "atol": 1e-6},  # k = 0 or T = 1
                    {"name": "spin_gap", "value": float(w1[0] - E0), "atol": 1e-8}],
        "source": "solutions/maple_leaf.py ed", "official": False,
    }


def theta_reference():
    ratio = tl.detection_volume(2.0) / tl.detection_volume(1.0)
    return {"answers": [{"name": "V2_over_V1", "value": round(ratio, 1), "atol": 0.051}],
            "source": "solutions/theta_logistic.py", "official": False}


def main():
    refs = {
        "maple-leaf-dirac-nodes": dirac_reference(),
        "maple-leaf-chern-number": chern_reference(),
        "maple-leaf-exact-diagonalization": ed_reference(),
        "theta-logistics": theta_reference(),
    }
    out = HERE / "references.local.json"
    out.write_text(json.dumps(refs, indent=1))
    print(f"wrote {len(refs)} references to {out}")


if __name__ == "__main__":
    main()
