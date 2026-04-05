#!/usr/bin/env python3
"""
export_dashboard_json.py
========================
Parse LaTeX output from Cadabra2/Srednicki export scripts and produce
a single chapters.json that qft-dashboard can consume at runtime.

Usage:
    python3 export_dashboard_json.py [--output /path/to/chapters.json]

Default output: ../../../JavaScript/qft-dashboard/public/chapters.json
"""

import re
import json
import sys
import os
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Optional

SR_DIR = Path(__file__).parent

# ── Chapter metadata (mirrors what was in chapters.ts) ───────────────────────
CHAPTER_META = {
    "ch34_left_right_spinors": {
        "id": "ch34", "number": "34",
        "title": "Left/Right Weyl Spinors",
        "status": "done",
        "scriptFile": "ch34_left_right_spinors.py",
        "description": (
            "Introduces two-component Weyl spinor formalism. Defines undotted "
            "(left-handed) and dotted (right-handed) spinor indices, the epsilon "
            "tensor, van der Waerden notation. Mostly-plus metric g=diag(-1,+1,+1,+1)."
        ),
    },
    "ch35_sigma_algebra": {
        "id": "ch35", "number": "35",
        "title": "Sigma Matrix Algebra",
        "status": "done",
        "scriptFile": "ch35_sigma_algebra.py",
        "description": (
            "Sigma and sigma-bar matrix identities, trace relations, completeness. "
            "Combines Weyl spinors into Dirac spinors and derives gamma matrices."
        ),
    },
    "ch36_weyl_lagrangian": {
        "id": "ch36", "number": "36",
        "title": "Weyl Lagrangian & Symmetries",
        "status": "done",
        "scriptFile": "ch36_weyl_lagrangian.py",
        "description": (
            "Symmetries of the Weyl Lagrangian: U(1) fermion number, C, P, T. "
            "Noether current for fermion number. Mass terms."
        ),
    },
    "ch37_canonical_quantization": {
        "id": "ch37", "number": "37",
        "title": "Canonical Quantization of Weyl Fields",
        "status": "done",
        "scriptFile": "ch37_canonical_quantization.py",
        "description": (
            "Promotes Weyl field to quantum operator via canonical anticommutation. "
            "Creation/annihilation operators, Fock space, fermion propagator, "
            "spin-sum completeness, Feynman slash."
        ),
    },
    "ch38_spinor_technology": {
        "id": "ch38", "number": "38",
        "title": "Spinor Technology & Feynman Rules",
        "status": "done",
        "scriptFile": "ch38_spinor_technology.py",
        "description": (
            "Spinor algebra toolkit: sigma-completeness, trace identities, "
            "Fierz rearrangement. Feynman rules for spinor QED."
        ),
    },
}

# Chapters without .tex files (future / not-started)
STUB_CHAPTERS = [
    {
        "id": "ch39", "number": "39",
        "title": "Feynman Rules for Dirac Fields",
        "status": "not-started",
        "scriptFile": "ch39_feynman_rules_dirac.py",
        "description": "Dirac propagator, QED Feynman rules, spin sums, trace theorems.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch40", "number": "40",
        "title": "Spin Sums",
        "status": "not-started",
        "scriptFile": "ch40_spin_sums.py",
        "description": "Spin-averaging, trace theorems, helicity vs trace approach.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch41", "number": "41",
        "title": "Gamma Matrix Technology",
        "status": "not-started",
        "scriptFile": "ch41_gamma_technology.py",
        "description": "Weyl representation, Clifford algebra, chiral projectors, Fierz identities.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch42", "number": "42",
        "title": "Spinor-Helicity (Core of MHV)",
        "status": "not-started",
        "scriptFile": "ch42_spinor_helicity.py",
        "description": "THE foundational chapter: massless momentum as spinor outer product, angle/square brackets, little group scaling.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch43", "number": "43",
        "title": "Feynman Rules for Majorana Fields",
        "status": "not-started",
        "scriptFile": "ch43_feynman_rules_majorana.py",
        "description": "Majorana condition, propagator with 1/2 factor, four-fermion contact terms.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch48", "number": "48",
        "title": "Massless Particles & Spinor-Helicity",
        "status": "not-started",
        "scriptFile": "ch48_massless_spinor_helicity.py",
        "description": "Angle and square spinors for massless momenta, polarization vectors, Schouten identity.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch60", "number": "60",
        "title": "MHV Amplitudes (Parke-Taylor)",
        "status": "done",
        "scriptFile": "ch60_spinor_helicity.py",
        "description": "Parke-Taylor formula for MHV gluon scattering amplitudes.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "bcfw", "number": "BCFW",
        "title": "BCFW Recursion Relations",
        "status": "not-started",
        "scriptFile": "bcfw_recursion.py",
        "description": "On-shell recursion for tree-level amplitudes via complex momentum shifts.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "adscft", "number": "AdS/CFT",
        "title": "AdS/CFT Correspondence",
        "status": "not-started",
        "scriptFile": "adscft_intro.py",
        "description": "Maldacena duality, holographic dictionary, GKPW relation.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch76", "number": "76",
        "title": "Nonabelian Gauge Theory (SU(N))",
        "status": "not-started",
        "scriptFile": "ch76_nonabelian_gauge.py",
        "description": "SU(N) gauge fields, covariant derivative, field strength, QCD Lagrangian.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch77", "number": "77",
        "title": "Group Representations (SU(N))",
        "status": "not-started",
        "scriptFile": "ch77_group_representations.py",
        "description": "SU(N) reps, Casimir operators, color factor algebra, Fierz identities.",
        "equations": [],
        "sections": [],
    },
    {
        "id": "ch80", "number": "80",
        "title": "Feynman Rules for Nonabelian Gauge Theory",
        "status": "not-started",
        "scriptFile": "ch80_feynman_rules_nonabelian.py",
        "description": "Gluon propagator, ghost-gluon vertex, three/four-gluon vertices, color-ordered rules.",
        "equations": [],
        "sections": [],
    },
]


def parse_tex_file(tex_path: Path) -> dict:
    """Parse a single .tex file into structured chapter data."""
    with open(tex_path, "r") as f:
        content = f.read()

    stem = tex_path.stem  # e.g. "ch37_canonical_quantization"
    meta = CHAPTER_META.get(stem, {})
    if not meta:
        return None

    # ── Extract document structure ───────────────────────────────────────
    sections = []
    current_section = {"title": "", "subsections": [], "equations": []}

    lines = content.split("\n")
    in_equation = False
    eq_buffer = []

    for i, line in enumerate(lines):
        stripped = line.strip()

        # Section headers
        sec_match = re.match(r"\\section\{(.+?)\}", stripped)
        if sec_match:
            if current_section["title"]:
                sections.append(current_section)
            current_section = {
                "title": sec_match.group(1),
                "subsections": [],
                "equations": [],
            }
            continue

        subsec_match = re.match(r"\\subsection\{(.+?)\}", stripped)
        if subsec_match:
            current_section["subsections"].append(subsec_match.group(1))
            continue

    # Don't forget the last section
    if current_section["title"]:
        sections.append(current_section)

    # ── Extract tagged equations ─────────────────────────────────────────
    equations = []

    # Strategy: find \tag{...}, walk backward to find the equation environment
    tag_positions = list(re.finditer(r"\\tag\{([^}]+)\}", content))

    for tag_match in tag_positions:
        tag = tag_match.group(1)
        tag_pos = tag_match.start()

        # Walk backward to find the start of the equation environment
        # Look for \begin{equation}, \begin{align}, or \[
        search_start = max(0, tag_pos - 600)
        before = content[search_start:tag_pos]

        # Find the innermost equation environment start
        env_starts = []
        for env_match in re.finditer(
            r"\\begin\{(equation|align|gather)\}", before
        ):
            env_starts.append((env_match.start() + search_start, env_match.group(1)))

        if env_starts:
            env_abs_start, env_type = env_starts[-1]  # innermost
            # Find the matching \end{...}
            end_pattern = rf"\\end\{{{env_type}\}}"
            end_match = re.search(end_pattern, content[env_abs_start:])
            if end_match:
                env_abs_end = env_abs_start + end_match.end()
                raw_latex = content[env_abs_start:env_abs_end]

                # Strip the environment wrapper, keep inner content
                inner = re.sub(
                    rf"\\begin\{{{env_type}\}}", "", raw_latex
                )
                inner = re.sub(rf"\\end\{{{env_type}\}}", "", inner)
                inner = re.sub(r"\\tag\{[^}]+\}", "", inner)
                inner = inner.strip()

                # Find which section this belongs to
                section_title = ""
                for sec in sections:
                    # Check if section title appears before this equation
                    sec_pos = content.find(f"\\section{{{sec['title']}}}")
                    if sec_pos >= 0 and sec_pos < tag_pos:
                        section_title = sec["title"]

                equations.append({
                    "label": tag,
                    "latex": inner,
                    "section": section_title,
                })

    # ── Extract boxed equations (key results) ────────────────────────────
    for boxed_match in re.finditer(
        r"\\boxed\{((?:[^{}]|\{[^{}]*\})*)\}", content
    ):
        boxed_latex = boxed_match.group(1).strip()
        # Check if this boxed eq already has a tag nearby
        nearby = content[boxed_match.start():boxed_match.end() + 100]
        tag_in_nearby = re.search(r"\\tag\{([^}]+)\}", nearby)
        if not tag_in_nearby:
            equations.append({
                "label": "boxed",
                "latex": boxed_latex,
                "section": "",
            })

    # ── Extract verified results ─────────────────────────────────────────
    verifications = []
    for ver_match in re.finditer(r"\\verified\{([^}]+)\}", content):
        verifications.append(ver_match.group(1))

    # ── Build chapter dict ───────────────────────────────────────────────
    docker_base = (
        "docker run --rm -v "
        "/home/propdev/.openclaw/workspace/repos/Monoclaw:/Monoclaw "
        "cadabra2-ubuntu:24.04 python3 /Monoclaw/Python/Cadabra2/Srednicki"
    )

    chapter = {
        "id": meta["id"],
        "number": meta["number"],
        "title": meta["title"],
        "status": meta["status"],
        "scriptFile": meta["scriptFile"],
        "description": meta["description"],
        "equations": equations,
        "sections": [
            {"title": s["title"], "subsections": s["subsections"]}
            for s in sections
        ],
        "verifications": verifications,
        "texFile": tex_path.name,
        "dockerCmd": f"{docker_base}/{meta['scriptFile']}",
    }
    return chapter


def main():
    output_path = SR_DIR / ".." / ".." / ".." / "JavaScript" / "qft-dashboard" / "public" / "chapters.json"

    # Check for --output flag
    if "--output" in sys.argv:
        idx = sys.argv.index("--output")
        if idx + 1 < len(sys.argv):
            output_path = Path(sys.argv[idx + 1])

    # Parse all .tex files
    chapters = []
    for tex_file in sorted(SR_DIR.glob("ch*_*.tex")):
        chapter = parse_tex_file(tex_file)
        if chapter:
            chapters.append(chapter)
            print(f"  ✅ {chapter['id']}: {len(chapter['equations'])} equations, "
                  f"{len(chapter['sections'])} sections, "
                  f"{len(chapter['verifications'])} verifications")

    # Add stub chapters (no .tex yet)
    for stub in STUB_CHAPTERS:
        # Don't duplicate if we already parsed a .tex for this id
        if not any(c["id"] == stub["id"] for c in chapters):
            chapters.append(stub)
            print(f"  ⏳ {stub['id']}: stub (no .tex)")

    # Sort by chapter number
    def sort_key(ch):
        try:
            return int(ch["number"])
        except ValueError:
            return 9999

    chapters.sort(key=sort_key)

    # Build output
    result = {
        "_generated": "export_dashboard_json.py",
        "_source": str(SR_DIR),
        "roadmap": [ch["id"] for ch in chapters],
        "chapters": chapters,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\n📄 Wrote {output_path} ({len(chapters)} chapters)")


if __name__ == "__main__":
    main()
