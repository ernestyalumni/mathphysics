#!/usr/bin/env python3
"""YM-F01 experiment 2: can Goedel-Formalizer-V2-8B state our propositions in Lean?

Inputs are the natural-language statements in `informal-statements.json`. Each corresponds to
Lean statements already written and checked by hand (`reference`). Prompt and sampling follow the
model card. An answer is scored in two ways:
- compiles: its last ```lean4 block, after `import Mathlib`, compiles in Docker against our Mathlib
  (`sorry` in the proof is expected).
- faithful: graded by hand against the informal statement and the reference, recorded in
  `grades.json` as equivalent / stronger hypotheses / weaker hypotheses / wrong.
  The sheet written by `sheet` is the grading input.

Subcommands: generate, judge, sheet, report. Standard library only; resumable.
"""
import argparse, json, pathlib, re, sys, time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prover  # noqa: E402  (chat, judge_file)

SEED = 20261121
PROMPT = ("Please autoformalize the following natural language problem statement in Lean 4. "
          "Use the following theorem name: {name}\n"
          "The natural language statement is: \n"
          "{statement}"
          "Think before you provide the lean statement.")


def chat(url, prompt, max_tokens, seed):
    import urllib.request
    body = json.dumps({"messages": [{"role": "user", "content": prompt}], "max_tokens": max_tokens,
                       "temperature": 0.9, "top_p": 0.95, "top_k": 20, "min_p": 0.0, "seed": seed}).encode()
    request = urllib.request.Request(url + "/v1/chat/completions", body, {"Content-Type": "application/json"})
    began = time.time()
    with urllib.request.urlopen(request, timeout=7200) as response:
        answer = json.load(response)
    choice = answer["choices"][0]
    return {"content": choice["message"].get("content") or "", "reasoning": choice["message"].get("reasoning_content") or "",
            "finish_reason": choice.get("finish_reason"), "usage": answer.get("usage"), "seconds": time.time() - began}


def extract(text):
    blocks = re.findall(r"```lean4\n(.*?)\n```", text, re.S)
    return blocks[-1].strip() if blocks else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["generate", "judge", "sheet", "report"])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--url", default="http://127.0.0.1:8084")
    parser.add_argument("--attempts", type=int, default=4)
    parser.add_argument("--max-tokens", type=int, default=14336)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    statements = json.loads((HERE / "informal-statements.json").read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    answers_path = args.out / "answers.jsonl"
    answers = {}
    if answers_path.exists():
        for line in answers_path.read_text().splitlines():
            row = json.loads(line)
            answers[(row["name"], row["attempt"])] = row

    if args.command == "generate":
        total, count = len(statements) * args.attempts, 0
        for attempt in range(args.attempts):
            for index, item in enumerate(statements):
                count += 1
                if (item["name"], attempt) in answers:
                    continue
                seed = SEED + 1000 * attempt + index
                try:
                    result, error = chat(args.url, PROMPT.format(name=item["name"], statement=item["statement"]),
                                         args.max_tokens, seed), None
                except Exception as exc:
                    result, error = {}, repr(exc)
                row = dict(name=item["name"], attempt=attempt, seed=seed, error=error, **result)
                with answers_path.open("a") as handle:
                    handle.write(json.dumps(row, ensure_ascii=False) + "\n")
                print("%s %d/%d %s#%d tokens=%s finish=%s seconds=%.0f error=%s" % (
                    time.strftime("%H:%M:%S"), count, total, item["name"], attempt,
                    (result.get("usage") or {}).get("completion_tokens"), result.get("finish_reason"),
                    result.get("seconds", 0), error), flush=True)
        return

    verdicts_path = args.out / "compile.jsonl"
    verdicts = {}
    if verdicts_path.exists():
        for line in verdicts_path.read_text().splitlines():
            row = json.loads(line)
            verdicts[(row["name"], row["attempt"])] = row

    if args.command == "judge":
        candidates = args.out / "candidates"
        candidates.mkdir(exist_ok=True)
        for key, row in sorted(answers.items()):
            if key in verdicts:
                continue
            code = None if row.get("error") else extract(row.get("content", ""))
            verdict = {"name": key[0], "attempt": key[1], "code": code}
            if code is None:
                verdict.update(compiles=False, reason="no lean4 block")
            else:
                body = "\n".join(line for line in code.splitlines() if not line.startswith("import "))
                path = candidates / ("%s__%d.lean" % key)
                path.write_text("import Mathlib\n\n" + body + "\n")
                result = prover.judge_file(path)
                verdict.update(compiles=result["returncode"] == 0 and not result["errors"], reason="lean",
                               errors=result["errors"][:3], declares_name=bool(re.search(r"\b%s\b" % key[0], code)))
            with verdicts_path.open("a") as handle:
                handle.write(json.dumps(verdict, ensure_ascii=False) + "\n")
            print("%s#%d %s %s" % (key[0], key[1], "compiles" if verdict["compiles"] else "fails",
                                   (verdict.get("errors") or [verdict["reason"]])[0][:150]), flush=True)
        return

    if args.command == "sheet":
        lines = ["# Formalizer grading sheet\n"]
        for item in statements:
            lines.append("## %s (%s)\n\nInformal: %s\n\nReference: %s\n" % (
                item["name"], item["source"], item["statement"], ", ".join(item["reference"])))
            for (name, attempt), verdict in sorted(verdicts.items()):
                if name == item["name"]:
                    lines.append("### attempt %d: %s\n\n```lean\n%s\n```\n" % (
                        attempt, "compiles" if verdict["compiles"] else "does not compile (%s)" % (
                            (verdict.get("errors") or [verdict["reason"]])[0][:200]), verdict.get("code") or ""))
        (args.out / "sheet.md").write_text("\n".join(lines))
        print(args.out / "sheet.md")
        return

    if args.command == "report":
        grades = json.loads((args.out / "grades.json").read_text()) if (args.out / "grades.json").exists() else {}
        rows = []
        for item in statements:
            mine = [v for (n, _), v in sorted(verdicts.items()) if n == item["name"]]
            graded = [grades.get("%s#%d" % (item["name"], v["attempt"])) for v in mine]
            rows.append({"name": item["name"], "attempts": len(mine), "compiles": sum(v["compiles"] for v in mine),
                         "faithful": sum(1 for g in graded if g and g["grade"] == "equivalent"),
                         "grades": graded})
        (args.out / "report.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False))
        for r in rows:
            print("%-18s compiles %d/%d  faithful %d/%d" % (r["name"], r["compiles"], r["attempts"], r["faithful"], r["attempts"]))


if __name__ == "__main__":
    sys.exit(main())
