"""Build problems.jsonl from the TheoPhysEval arXiv source (arXiv:2610.16733).

The appendix carries a lightly edited copy of every prompt; the exact prompts, data files and
answer keys are in the gated Hugging Face bundle (RameshChandra1/TheoPhysEval). Each row:

  id          problem slug from problem_names.tex (e.g. "maple-leaf-chern-number")
  title       short name used in the paper's tables
  label       LaTeX label in appendix.tex
  theme       appendix subsection (one of 11 themes)
  statement   LaTeX source of the appendix prompt
  paper_grades {agent: {"R": 0|1, "I": 0|1}} from the paper's per-problem table
               (R = realistic track with web, I = isolated track)

Run: python extract_problems.py [--source DIR] [--out problems.jsonl]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

DEFAULT_SOURCE = Path("/media/propdev/Expansion/openclaw/.openclaw/workspace/Data/Public/papers/AI/"
                      "arXiv-2610.16733v1")
AGENTS = ["Astra", "Opus", "Gemini", "DeepSeek", "GLM"]


def problem_names(source: Path) -> dict[str, tuple[str, str]]:
    """slug -> (title, label) from \\DeclareProblemName{slug}{title}{label}."""
    text = (source / "problem_names.tex").read_text()
    out = {}
    for slug, title, label in re.findall(r"\\DeclareProblemName\{([^}]*)\}\{([^}]*)\}\{([^}]*)\}", text):
        out[slug] = (title, label.strip())
    return out


def paper_grades(source: Path) -> dict[str, dict[str, dict[str, int]]]:
    text = (source / "detailed_scores.tex").read_text()
    out = {}
    for slug, cells in re.findall(r"\\problemname\{([^}]*)\}((?:\s*&\s*[01]){10})", text):
        bits = [int(b) for b in re.findall(r"[01]", cells)]
        out[slug] = {a: {"R": bits[2 * k], "I": bits[2 * k + 1]} for k, a in enumerate(AGENTS)}
    return out


def split_appendix(source: Path) -> list[dict]:
    text = (source / "appendix.tex").read_text()
    # Cut at the electronic-structure / later themes too; every problem is a \subsubsection.
    pieces = re.split(r"(\\subsection\{[^}]*\}|\\subsubsection\{[^\n]*\})", text)
    theme, rows, current = None, [], None
    for piece in pieces:
        m_sub = re.match(r"\\subsection\{([^}]*)\}", piece)
        m_subsub = re.match(r"\\subsubsection\{(.*)\}$", piece)
        if m_sub:
            theme = m_sub.group(1)
            current = None
        elif m_subsub:
            current = {"heading": m_subsub.group(1), "theme": theme, "body": ""}
            rows.append(current)
        elif current is not None:
            current["body"] += piece
    for r in rows:
        m = re.search(r"\\label\{([^}]*)\}", r["body"])
        r["label"] = m.group(1).strip() if m else None
        body = re.sub(r"\\label\{[^}]*\}\s*", "", r["body"], count=1)
        body = body.split("\\endgroup")[0]
        r["statement"] = body.strip()
    return rows


def build(source: Path) -> list[dict]:
    names = problem_names(source)
    grades = paper_grades(source)
    by_label = {label: slug for slug, (_, label) in names.items()}
    rows = []
    for r in split_appendix(source):
        slug = by_label.get(r["label"])
        if slug is None:
            continue                                   # intro text, not a problem
        title, _ = names[slug]
        rows.append({
            "id": slug,
            "title": title,
            "label": r["label"],
            "theme": r["theme"],
            "heading": r["heading"],
            "statement": r["statement"],
            "paper_grades": grades.get(slug),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("problems.jsonl"))
    args = ap.parse_args()
    rows = build(args.source)
    with args.out.open("w") as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    missing = set(problem_names(args.source)) - {r["id"] for r in rows}
    print(f"wrote {len(rows)} problems to {args.out}; missing: {sorted(missing) or 'none'}")


if __name__ == "__main__":
    main()
