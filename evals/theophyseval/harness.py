"""Run TheoPhysEval problems through an agent and grade the answers numerically.

Backends
  claude   Claude Code headless (`claude -p`), with tools, in a fresh work directory per
           problem. Two tracks as in the paper: --track realistic (web allowed) or
           --track isolated (WebSearch/WebFetch disallowed). The paper's isolated track also
           cut the network inside a Docker sandbox; here the agent's Bash still has network.
  openai   any OpenAI-compatible chat endpoint (e.g. the local llama-cpp server) in a single
           closed-book turn with no tools: a lower baseline than either paper track.

Safety: the claude backend gives the agent Bash on this machine and costs real money, so
nothing runs unless --execute is passed; the default prints the plan. --budget-usd is passed
through as `claude --max-budget-usd` for every problem.

Prompts: by default the appendix text from problems.jsonl (lightly edited by the authors).
With --task-dir pointing at the downloaded Hugging Face bundle, the official
tasks-<track>/<id>-<track>/instruction.md is used instead.

Results: results/YYYY-MM-DD-<backend>-<model>-<track>.json, one row per problem in the format
of agent/EVAL-AND-TOOLING.md, graded against references.local.json when an entry exists.

Examples
  python harness.py --list
  python harness.py --problems maple-leaf-chern-number                      # dry run
  python harness.py --problems maple-leaf-chern-number --execute --model sonnet --budget-usd 3
  python harness.py --backend openai --base-url http://localhost:8080/v1 --model local \\
                    --problems theta-logistics --execute
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

from verifier import grade

HERE = Path(__file__).resolve().parent

PREAMBLE = (
    "You are solving a research-level theoretical physics problem autonomously. Work in the "
    "current directory; you may write and run scripts. Check your result independently "
    "(convergence, a second method, limiting cases) before reporting it. If you cannot obtain "
    "a reliable answer, say so explicitly instead of guessing. Use exactly the parameter values "
    "and conventions in the prompt, even if a reference you know uses different ones.\n\n"
    "End your response with a single line\n"
    "FINAL ANSWER: \\boxed{...}\n"
    "containing every requested number.\n\n"
    "Problem:\n\n"
)


def load_problems() -> dict[str, dict]:
    path = HERE / "problems.jsonl"
    if not path.exists():
        sys.exit("problems.jsonl missing: run `python extract_problems.py` first")
    return {r["id"]: r for r in map(json.loads, path.read_text().splitlines())}


def load_references() -> dict[str, dict]:
    path = HERE / "references.local.json"
    return json.loads(path.read_text()) if path.exists() else {}


def prompt_for(problem: dict, track: str, task_dir: Path | None) -> str:
    if task_dir:
        inst = task_dir / f"tasks-{track}" / f"{problem['id']}-{track}" / "instruction.md"
        if inst.exists():
            return inst.read_text()
        print(f"  [warn] {inst} not found; falling back to the appendix text", file=sys.stderr)
    return PREAMBLE + problem["statement"]


def claude_command(prompt: str, model: str, track: str, budget: float, effort: str | None) -> list[str]:
    cmd = ["claude", "-p", prompt, "--output-format", "json", "--model", model,
           "--max-budget-usd", f"{budget:g}", "--permission-mode", "acceptEdits",
           "--allowedTools", "Bash", "Read", "Write", "Edit", "Glob", "Grep"]
    if track == "realistic":
        cmd += ["WebSearch", "WebFetch"]
    else:
        cmd += ["--disallowedTools", "WebSearch", "WebFetch"]
    if effort:
        cmd += ["--effort", effort]
    return cmd


def run_claude(prompt, model, track, budget, effort, workdir: Path, timeout_s: int) -> dict:
    cmd = claude_command(prompt, model, track, budget, effort)
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=timeout_s)
    latency = int(1000 * (time.time() - t0))
    (workdir / "agent_stdout.json").write_text(proc.stdout)
    (workdir / "agent_stderr.txt").write_text(proc.stderr)
    try:
        out = json.loads(proc.stdout)
    except json.JSONDecodeError:
        out = {"result": proc.stdout}
    usage = out.get("usage") or {}
    return {
        "model_answer": out.get("result") or "",
        "latency_ms": latency,
        "tokens_in": usage.get("input_tokens"),
        "tokens_cache_read": usage.get("cache_read_input_tokens"),
        "tokens_out": usage.get("output_tokens"),
        "cost_usd": out.get("total_cost_usd"),
        "num_turns": out.get("num_turns"),
        "exit_code": proc.returncode,
    }


def run_openai(prompt, model, base_url, workdir: Path, timeout_s: int) -> dict:
    body = json.dumps({"model": model, "messages": [{"role": "user", "content": prompt}],
                       "temperature": 0}).encode()
    req = urllib.request.Request(base_url.rstrip("/") + "/chat/completions", data=body,
                                 headers={"Content-Type": "application/json"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        out = json.loads(resp.read())
    (workdir / "agent_response.json").write_text(json.dumps(out, indent=1))
    usage = out.get("usage") or {}
    return {"model_answer": out["choices"][0]["message"]["content"],
            "latency_ms": int(1000 * (time.time() - t0)),
            "tokens_in": usage.get("prompt_tokens"), "tokens_out": usage.get("completion_tokens"),
            "cost_usd": 0.0}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="list problems with paper grades and exit")
    ap.add_argument("--problems", nargs="*", default=[], help="problem ids (default: those with references)")
    ap.add_argument("--backend", choices=["claude", "openai"], default="claude")
    ap.add_argument("--model", default="sonnet")
    ap.add_argument("--track", choices=["realistic", "isolated"], default="isolated")
    ap.add_argument("--effort", default=None, help="claude --effort level (the paper used max)")
    ap.add_argument("--budget-usd", type=float, default=2.0, help="per-problem cap (claude backend)")
    ap.add_argument("--timeout-min", type=int, default=120)
    ap.add_argument("--base-url", default="http://localhost:8080/v1")
    ap.add_argument("--task-dir", type=Path, default=None, help="downloaded HF bundle root")
    ap.add_argument("--execute", action="store_true", help="actually run the agent")
    args = ap.parse_args()

    problems, refs = load_problems(), load_references()
    if args.list:
        for pid, p in problems.items():
            g = p["paper_grades"]
            row = " ".join(f"{a}:{g[a]['R']}{g[a]['I']}" for a in g)
            print(f"{pid:40s} {'REF' if pid in refs else '   '} {p['theme'][:28]:28s} {row}")
        return

    ids = args.problems or sorted(refs)
    unknown = [i for i in ids if i not in problems]
    if unknown:
        sys.exit(f"unknown problem ids: {unknown}")

    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    run_root = HERE / "runs" / stamp
    rows = []
    for pid in ids:
        prompt = prompt_for(problems[pid], args.track, args.task_dir)
        print(f"== {pid} ({len(prompt)} chars, reference: {'yes' if pid in refs else 'no'})")
        if not args.execute:
            if args.backend == "claude":
                cmd = claude_command("<prompt>", args.model, args.track, args.budget_usd, args.effort)
                print("   dry run:", " ".join(cmd), f"(cwd: runs/{stamp}/{pid})")
            else:
                print(f"   dry run: POST {args.base_url}/chat/completions model={args.model}")
            continue
        if args.backend == "claude" and shutil.which("claude") is None:
            sys.exit("claude CLI not on PATH")
        workdir = run_root / pid
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "prompt.md").write_text(prompt)
        try:
            if args.backend == "claude":
                res = run_claude(prompt, args.model, args.track, args.budget_usd, args.effort,
                                 workdir, 60 * args.timeout_min)
            else:
                res = run_openai(prompt, args.model, args.base_url, workdir, 60 * args.timeout_min)
        except subprocess.TimeoutExpired:
            res = {"model_answer": "", "latency_ms": 60_000 * args.timeout_min, "timeout": True}
        row = {"problem_id": pid, **res}
        row.update(grade(res["model_answer"], refs[pid]) if pid in refs
                   else {"verifier_status": "UNGRADED"})
        rows.append(row)
        print(f"   {row['verifier_status']}  cost={row.get('cost_usd')}  {row.get('latency_ms')} ms")

    if args.execute and rows:
        out = HERE / "results" / f"{dt.date.today()}-{args.backend}-{args.model}-{args.track}.json"
        out.parent.mkdir(exist_ok=True)
        graded = [r for r in rows if r["verifier_status"] != "UNGRADED"]
        summary = {
            "harness": "evals/theophyseval/harness.py", "backend": args.backend, "model": args.model,
            "track": args.track, "effort": args.effort, "references": "self-computed (not official)",
            "accuracy": sum(r["verifier_status"] == "PASS" for r in graded) / max(len(graded), 1),
            "confident_wrong": sum(r.get("confident_wrong", False) for r in graded),
            "honest_unknown": sum(r.get("admitted_unknown", False) for r in graded),
            "total_cost_usd": sum((r.get("cost_usd") or 0) for r in rows),
            "rows": rows,
        }
        out.write_text(json.dumps(summary, indent=1))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
