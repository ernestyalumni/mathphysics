#!/usr/bin/env python3
"""YM-F01 step 3: can a local Lean prover fill in the YMFormal proofs? Lean is the only judge.

For each theorem of the frozen statement commit, the model sees the module up to that theorem
(earlier declarations with their proofs, imported modules with proofs omitted) and the
theorem with `sorry`. From its answer only the proof of that theorem, and any new helper
`theorem`/`lemma` it writes, are kept. They are placed after our own prefix and under our own
statement text, so the model cannot change what is proved. A candidate passes when it compiles
in Docker (no network) and `#print axioms` lists only propext, Classical.choice and
Quot.sound: `sorry`, new axioms and `native_decide` all fail.

Subcommands:
  tasks    --statements-rev REV --out tasks.json
  generate --tasks tasks.json --out DIR --url URL --attempts K --max-tokens N
  judge    --tasks tasks.json --out DIR [--jobs J]
  report   --tasks tasks.json --out DIR
Standard library only. Resumable: finished attempts and verdicts are skipped.
"""
import argparse, concurrent.futures, hashlib, json, os, pathlib, re, subprocess, sys, time, urllib.request

HERE = pathlib.Path(__file__).resolve().parent
MODULES = ["YMFormal.Purity", "YMFormal.Defect", "YMFormal.ConditionalGap", "YMFormal.Casimir"]
IMAGE = os.environ.get("YM_FORMAL_IMAGE", "claw-work/lean-agent:lean4.35.0-mathlib")
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
SEED = 20261120
# Prompt of the Goedel-Prover-V2 model card, verbatim around the code.
PROMPT = ("Complete the following Lean 4 code:\n\n```lean4\n{code}```\n\n"
          "Before producing the Lean 4 code to formally prove the given theorem, provide a detailed proof plan "
          "outlining the main proof steps and strategies.\nThe plan should highlight key ideas, intermediate "
          "lemmas, and proof structures that will guide the construction of the final formal proof.")
DECLARATION = re.compile(r"^(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+)?(?:noncomputable\s+)?"
                         r"(theorem|lemma|def|abbrev|structure|instance|example|axiom|opaque)\b\s*([^\s:({\[]*)", re.M)
TOP_LEVEL = re.compile(r"^(?:theorem|lemma|def|abbrev|structure|instance|example|axiom|opaque|noncomputable|"
                       r"namespace|section|end\b|open|variable|#|/--|/-!|@\[|import|private|protected)", re.M)


def module_path(module):
    return HERE / (module.replace(".", "/") + ".lean")


def git_show(rev, path):
    relative = path.relative_to(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=HERE, text=True).strip())
    return subprocess.check_output(["git", "show", "%s:%s" % (rev, relative)], cwd=HERE, text=True)


def theorem_spans(text):
    """(name, start of its docstring or declaration, start of `theorem`) for each theorem."""
    spans = []
    for match in re.finditer(r"^theorem\s+(\S+)", text, re.M):
        start = match.start()
        doc = text.rfind("/--", 0, start)
        if doc != -1 and text[doc:start].rstrip().endswith("-/") and "\ntheorem" not in text[doc:start]:
            start = doc
        spans.append((match.group(1), start, match.start()))
    return spans


def statement_only(text):
    """The module with every proof replaced by `sorry` (for imported context)."""
    return re.sub(r":= by\n(?:  .*\n|\n)*?(?=\n?(?:/--|theorem|end |noncomputable|def|structure))", ":= by\n  sorry\n", text)


def make_tasks(rev):
    tasks = []
    for module in MODULES:
        frozen = git_show(rev, module_path(module))
        proved = module_path(module).read_text()
        imported = [m for m in MODULES if ("import " + m) in proved]
        context = "".join("-- Contents of %s (imported; proofs omitted)\n%s\n" % (m, git_show(rev, module_path(m)).split("\n", 1)[1])
                          for m in imported)
        proved_spans = {name: start for name, start, _ in theorem_spans(proved)}
        for name, start, decl in theorem_spans(frozen):
            end = frozen.index(":= by\n  sorry", decl)
            statement = frozen[start:end].rstrip()
            prefix = proved[:proved_spans[name]]
            tasks.append({"id": "%s.%s" % (module, name), "module": module, "name": name, "statement": statement,
                          "prefix": prefix,
                          "prompt": PROMPT.format(code=context + prefix + statement + " := by\n  sorry\n\nend YMFormal\n")})
    return tasks


def chat(url, prompt, max_tokens, seed):
    body = json.dumps({"messages": [{"role": "user", "content": prompt}], "max_tokens": max_tokens,
                       "temperature": 0.6, "top_p": 0.95, "top_k": 20, "min_p": 0.0, "seed": seed}).encode()
    request = urllib.request.Request(url + "/v1/chat/completions", body, {"Content-Type": "application/json"})
    began = time.time()
    with urllib.request.urlopen(request, timeout=7200) as response:
        answer = json.load(response)
    choice = answer["choices"][0]
    return {"content": choice["message"].get("content") or "", "reasoning": choice["message"].get("reasoning_content") or "",
            "finish_reason": choice.get("finish_reason"), "usage": answer.get("usage"), "seconds": time.time() - began}


def last_code_block(text):
    blocks = re.findall(r"```(?:lean4|lean)?\s*\n(.*?)```", text, re.S)
    return blocks[-1] if blocks else None


def split_top_level(code):
    """Top-level blocks (start offset, text) of a Lean source."""
    starts = [m.start() for m in TOP_LEVEL.finditer(code)] + [len(code)]
    return [(a, code[a:b]) for a, b in zip(starts, starts[1:])]


def assemble(task, answer):
    """Candidate file from our prefix and statement and the model's proof, or (None, reason)."""
    code = last_code_block(answer)
    if code is None:
        return None, "no code block"
    declared = {m.group(2).split(".")[-1] for m in DECLARATION.finditer(task["prefix"])}
    helpers, proof = [], None
    for _, block in split_top_level(code):
        match = DECLARATION.match(block)
        if not match:
            continue
        kind, name = match.group(1), match.group(2).split(".")[-1]
        if name == task["name"] and kind in ("theorem", "lemma"):
            at = block.find(":=")
            if at == -1:
                return None, "target has no proof"
            proof = block[at + 2:].rstrip()
        elif kind in ("theorem", "lemma") and name and name not in declared:
            helpers.append(block.rstrip())
    if proof is None:
        return None, "target theorem not found in the answer"
    head = task["statement"]
    if head.endswith(":= by"):
        head = head[:-len(":= by")].rstrip()
    text = (task["prefix"] + "".join(h + "\n\n" for h in helpers) + head + " :=" + proof + "\n\n"
            + "#print axioms YMFormal.%s\n\nend YMFormal\n" % task["name"])
    return text, None


def judge_file(path):
    command = ["docker", "run", "--rm", "--network", "none", "-u", "%d:%d" % (os.getuid(), os.getgid()),
               "-v", "%s:/work:ro" % HERE, "-v", "%s:/candidate:ro" % path.parent, "-w", "/work", "--entrypoint", "bash",
               IMAGE, "-c", "LEAN_PATH=/work/build:$LEAN_PATH timeout 600 lean /candidate/%s" % path.name]
    run = subprocess.run(command, capture_output=True, text=True)
    output = run.stdout + run.stderr
    errors = [line for line in output.splitlines() if ": error" in line]
    axioms = re.search(r"depends on axioms: \[([^\]]*)\]", output)
    used = set(a.strip() for a in axioms.group(1).split(",")) if axioms else set()
    no_axioms = "does not depend on any axioms" in output
    passed = run.returncode == 0 and not errors and (no_axioms or (axioms is not None and used <= ALLOWED_AXIOMS))
    return {"passed": passed, "returncode": run.returncode, "errors": errors[:5], "axioms": sorted(used),
            "output_tail": output[-1500:]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["tasks", "generate", "judge", "report"])
    parser.add_argument("--statements-rev")
    parser.add_argument("--tasks", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path)
    parser.add_argument("--url", default="http://127.0.0.1:8098")
    parser.add_argument("--attempts", type=int, default=4)
    parser.add_argument("--max-tokens", type=int, default=16384)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()

    if args.command == "tasks":
        tasks = make_tasks(args.statements_rev)
        args.out.write_text(json.dumps(tasks, indent=1, ensure_ascii=False))
        print("%d tasks, sha256 %s" % (len(tasks), hashlib.sha256(args.out.read_bytes()).hexdigest()))
        return
    tasks = json.loads(args.tasks.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    answers = args.out / "answers.jsonl"
    done = {}
    if answers.exists():
        for line in answers.read_text().splitlines():
            row = json.loads(line)
            done[(row["task"], row["attempt"])] = row

    if args.command == "generate":
        total = len(tasks) * args.attempts
        count = 0
        for attempt in range(args.attempts):          # attempt-major: every theorem gets attempt 0 first
            for index, task in enumerate(tasks):
                count += 1
                if (task["id"], attempt) in done:
                    continue
                seed = SEED + 1000 * attempt + index
                try:
                    result = chat(args.url, task["prompt"], args.max_tokens, seed)
                    error = None
                except Exception as exc:  # recorded, not retried: the run is resumable
                    result, error = {}, repr(exc)
                row = dict(task=task["id"], attempt=attempt, seed=seed, error=error, **result)
                with answers.open("a") as handle:
                    handle.write(json.dumps(row, ensure_ascii=False) + "\n")
                usage = (result.get("usage") or {}).get("completion_tokens")
                print("%s %d/%d %s#%d tokens=%s finish=%s seconds=%.0f error=%s" % (
                    time.strftime("%H:%M:%S"), count, total, task["id"], attempt, usage, result.get("finish_reason"),
                    result.get("seconds", 0), error), flush=True)
        return

    verdicts_path = args.out / "verdicts.jsonl"
    verdicts = {}
    if verdicts_path.exists():
        for line in verdicts_path.read_text().splitlines():
            row = json.loads(line)
            verdicts[(row["task"], row["attempt"])] = row
    by_id = {task["id"]: task for task in tasks}

    if args.command == "judge":
        candidates = args.out / "candidates"
        candidates.mkdir(exist_ok=True)
        pending = []
        for key, row in sorted(done.items()):
            if key in verdicts:
                continue
            text, reason = (None, "generation error") if row.get("error") else assemble(by_id[row["task"]], row.get("content", ""))
            if text is None:
                verdict = dict(task=key[0], attempt=key[1], passed=False, reason=reason)
                with verdicts_path.open("a") as handle:
                    handle.write(json.dumps(verdict) + "\n")
                continue
            path = candidates / ("%s__%d.lean" % (row["task"].replace(".", "_"), key[1]))
            path.write_text(text)
            pending.append((key, path))
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
            for (key, path), result in zip(pending, pool.map(lambda item: judge_file(item[1]), pending)):
                verdict = dict(task=key[0], attempt=key[1], reason="lean", **result)
                with verdicts_path.open("a") as handle:
                    handle.write(json.dumps(verdict) + "\n")
                print("%s#%d %s %s" % (key[0], key[1], "PASS" if result["passed"] else "fail",
                                       (result["errors"] or [""])[0][:160]), flush=True)
        return

    if args.command == "report":
        rows = []
        for task in tasks:
            attempts = [verdicts.get((task["id"], k)) for k in sorted({a for (t, a) in done if t == task["id"]})]
            passed = sum(1 for v in attempts if v and v["passed"])
            reasons = sorted({v["reason"] if v["reason"] != "lean" else "lean error" for v in attempts if v and not v["passed"]})
            rows.append({"task": task["id"], "attempts": len(attempts), "passed": passed, "failures": reasons})
        solved = sum(1 for r in rows if r["passed"] > 0)
        summary = {"theorems": len(rows), "solved_any_attempt": solved, "rows": rows}
        (args.out / "report.json").write_text(json.dumps(summary, indent=1))
        for r in rows:
            print("%-48s %d/%d %s" % (r["task"], r["passed"], r["attempts"], ", ".join(r["failures"])))
        print("solved (pass@%d): %d of %d" % (max((r["attempts"] for r in rows), default=0), solved, len(rows)))


if __name__ == "__main__":
    sys.exit(main())
