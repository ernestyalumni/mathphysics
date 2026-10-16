"""Differential parity with the pre-existing NumPy implementation, not a proof.
Usage: python python_parity.py /absolute/path/to/qft-symbolics
Uses existing checkout/package, never installs dependencies or modifies old code.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT / "amplitudes/MHVamplitudes/src"
sys.path.insert(0, str(PACKAGE))
from mhvamplitudes.amplitudes.smga import smga_formula_r1
from mhvamplitudes.amplitudes.berends_giele import compute_stripped_amplitude

binary = Path(sys.argv[1]).resolve()
assert binary.is_file()
total = 0
zero = nonzero = 0
for n in range(3, 11):
    for seed in range(25):
        row = json.loads(subprocess.check_output(
            [str(binary), "smga-fixture", str(n), str(seed)], text=True
        ))
        # These fixtures have integral w and w*z. Keep conversion explicit.
        vectors = [[float(x) for x in p] for p in row["vectors"]]
        w = [v[0] for v in vectors]
        z = [v[1] / v[0] for v in vectors]
        old_formula = smga_formula_r1(n, w, z)
        old_recursion = compute_stripped_amplitude(n, w, z)
        new_formula = float(row["formula"])
        new_recursion = float(row["recursion"])
        assert row["status"] == "PASS", row
        assert abs(old_formula - new_formula) <= 1e-10, (n, seed, old_formula, row)
        assert abs(old_recursion - new_recursion) <= 1e-10, (n, seed, old_recursion, row)
        zero += new_formula == 0
        nonzero += new_formula != 0
        total += 1
source_hashes = {
    str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in [PACKAGE / "mhvamplitudes/amplitudes/smga.py",
              PACKAGE / "mhvamplitudes/amplitudes/berends_giele.py",
              PACKAGE / "mhvamplitudes/kinematics/spinor_helicity.py",
              PACKAGE / "mhvamplitudes/kinematics/phase_space.py"]
}
print(json.dumps({
    "status": "PASS", "cases": total, "zero": zero, "nonzero": nonzero,
    "multiplicities": [3, 10], "seeds": [0, 24], "tolerance": 1e-10,
    "sources": source_hashes,
    "scope": "Rust exact vs pre-existing NumPy paths at finite points only",
}, indent=2))
