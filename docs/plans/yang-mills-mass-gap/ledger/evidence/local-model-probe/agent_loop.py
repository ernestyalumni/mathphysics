#!/usr/bin/env python3
"""Tool-use and loop tests for a local OpenAI-compatible model.

Usage: agent_loop.py --lean <LEAN checkout> --out <directory outside git> [--endpoint http://127.0.0.1:8080/v1]
                     [--runs 3] [--only name,name]
Standard library only. The model can call only the fixed tools below; nothing it writes is executed
as code, and the file tools read two whitelisted files. Endpoint must be loopback or private-range.
"""
import argparse, hashlib, ipaddress, json, pathlib, re, time, urllib.parse, urllib.request
from fractions import Fraction

def tool_schema(name, description, properties, required):
    return {"type": "function", "function": {"name": name, "description": description,
            "parameters": {"type": "object", "properties": properties, "required": required}}}

INTEGER = {"type": "integer"}
TOOLS = [
    tool_schema("add", "Add two integers exactly.", {"a": INTEGER, "b": INTEGER}, ["a", "b"]),
    tool_schema("multiply", "Multiply two integers exactly.", {"a": INTEGER, "b": INTEGER}, ["a", "b"]),
    tool_schema("collatz_step", "One Collatz step: n/2 if n is even, 3n+1 if n is odd.", {"n": INTEGER}, ["n"]),
    tool_schema("read_lines", "Read lines start..end (1-based, inclusive, at most 40 lines) of an allowed file.",
                {"file": {"type": "string", "enum": ["ChallengeSpec.lean", "Solution.lean"]}, "start": INTEGER, "end": INTEGER},
                ["file", "start", "end"]),
    tool_schema("find_lines", "Return the 1-based line numbers of an allowed file that contain the exact text.",
                {"file": {"type": "string"}, "text": {"type": "string"}}, ["file", "text"]),
    tool_schema("sphere_monomial_average", "Exact average of n_1^e_1 ... n_d^e_d over the unit sphere of R^d, d = len(exponents). Returns a fraction.",
                {"exponents": {"type": "array", "items": INTEGER}}, ["exponents"]),
]
ALLOWED_FILES = ("ChallengeSpec.lean", "Solution.lean")

def chain_token(length, index):
    """Token number `index` of a hidden chain of the given length; cannot be guessed without calling."""
    return f"L{length}-{index}-" + hashlib.sha256(f"ym-chain:{length}:{index}".encode()).hexdigest()[:8]

TOOLS.append(tool_schema("next_token", "Return the token that follows the given token in a hidden sequence, or DONE when the sequence has ended. "
                         "Tokens cannot be computed or guessed; only this tool knows them.", {"token": {"type": "string"}}, ["token"]))

def run_tool(name, arguments, lean_root):
    try:
        if name == "add":
            return str(int(arguments["a"]) + int(arguments["b"]))
        if name == "multiply":
            return str(int(arguments["a"]) * int(arguments["b"]))
        if name == "collatz_step":
            n = int(arguments["n"])
            if n < 1:
                return "ERROR: n must be positive"
            return str(n // 2 if n % 2 == 0 else 3 * n + 1)
        if name in ("read_lines", "find_lines"):
            if arguments.get("file") not in ALLOWED_FILES:
                return "ERROR: file not allowed. Allowed files: " + ", ".join(ALLOWED_FILES)
            lines = pathlib.Path(lean_root, arguments["file"]).read_text().splitlines()
            if name == "find_lines":
                hits = [str(i + 1) for i, line in enumerate(lines) if str(arguments["text"]) in line]
                return ", ".join(hits[:50]) if hits else "no line contains that text"
            start, end = int(arguments["start"]), int(arguments["end"])
            if start < 1 or end < start or end - start >= 40:
                return "ERROR: need 1 <= start <= end and at most 40 lines"
            return "\n".join(f"{i}: {lines[i - 1]}" for i in range(start, min(end, len(lines)) + 1))
        if name == "sphere_monomial_average":
            exponents = [int(e) for e in arguments["exponents"]]
            if not exponents or any(e < 0 for e in exponents):
                return "ERROR: need a nonempty list of nonnegative integers"
            if any(e % 2 for e in exponents):
                return "0"
            numerator, half_degree = 1, 0
            for e in exponents:
                half_degree += e // 2
                for odd in range(1, e, 2):
                    numerator *= odd
            denominator = 1
            for step in range(half_degree):
                denominator *= len(exponents) + 2 * step
            return str(Fraction(numerator, denominator))
        if name == "next_token":
            match = re.fullmatch(r"L(\d+)-(\d+)-[0-9a-f]{8}", str(arguments["token"]))
            if not match or chain_token(int(match.group(1)), int(match.group(2))) != arguments["token"]:
                return "ERROR: unknown token"
            length, index = int(match.group(1)), int(match.group(2))
            return "DONE" if index + 1 >= length else chain_token(length, index + 1)
        return f"ERROR: unknown tool {name}"
    except (KeyError, TypeError, ValueError) as error:
        return f"ERROR: bad arguments ({error})"

TESTS = {
    "t1_single_tool": dict(prompt="What is 123456789 * 987654321? Use the multiply tool, then state the product.",
        check=lambda text, calls: "121932631112635269" in text.replace(",", ""), expect="121932631112635269, 1 call"),
    "t2_add_chain": dict(prompt="Using only the add tool (two numbers at a time), compute 17 + 29 + 46 + 58 + 73. Then give the total on a line 'TOTAL: <number>'.",
        check=lambda text, calls: bool(re.search(r"TOTAL:\s*223\b", text)), expect="223, 4 add calls"),
    "t3_collatz_loop": dict(prompt="Start from n = 7. Call collatz_step repeatedly, feeding each result back in, until the result is 1. Do not compute steps yourself. Then answer on a line 'STEPS: <number of calls you made>'.",
        check=lambda text, calls: bool(re.search(r"STEPS:\s*16\b", text)) and sum(1 for c in calls if c[0] == "collatz_step") == 16, expect="16 calls, STEPS: 16"),
    "t4_read_file": dict(prompt="Use read_lines on ChallengeSpec.lean, lines 186 to 188. Then answer in two lines: 'NAME: <name of the definition>' and 'FORALL: <how many times the symbol ∀ occurs in those lines>'.",
        check=lambda text, calls: "JointExistence" in text and bool(re.search(r"FORALL:\s*2\b", text)), expect="NAME: JointExistence, FORALL: 2"),
    "t5_search_then_read": dict(prompt="In Solution.lean, find the line that contains the text 'theorem' with find_lines, then read that single line with read_lines. Answer: 'LINE: <number>' and 'THEOREM: <full name of the theorem>'.",
        check=lambda text, calls: bool(re.search(r"LINE:\s*4\b", text)) and "YangMills.jointExistence" in text, expect="LINE: 4, THEOREM: YangMills.jointExistence"),
    "t6_tool_for_arithmetic": dict(prompt="Get the average of n1^2 n2^2 n3^2 over the unit sphere of R^6 with the sphere_monomial_average tool (exponents [2,2,2,0,0,0]). Then write one Rust line of the form assert_eq!(sphere_monomial_average(&[2, 2, 2, 0, 0, 0]).unwrap(), q(<numerator>, <denominator>));",
        check=lambda text, calls: bool(re.search(r"q\(\s*1\s*,\s*480\s*\)", text)), expect="q(1, 480)"),
    "t7_refused_file": dict(prompt="Read lines 1 to 3 of the file Secrets.lean with read_lines and quote them. If you cannot, say exactly 'CANNOT READ' and explain in one sentence.",
        check=lambda text, calls: "CANNOT READ" in text, expect="CANNOT READ (the tool refuses; no invented content)"),
    "t9_hidden_chain_12": dict(prompt=f"A hidden sequence starts with the token {chain_token(12, 0)}. Call next_token repeatedly, each time with the token you just received, until the tool returns DONE. Then answer 'LAST: <the last token you received before DONE>'.",
        check=lambda text, calls: chain_token(12, 11) in text and sum(1 for c in calls if c[0] == "next_token" and not c[2].startswith("ERROR")) == 12, expect="12 valid calls; LAST = token 11"),
    "t10_hidden_chain_25": dict(prompt=f"A hidden sequence starts with the token {chain_token(25, 0)}. Call next_token repeatedly, each time with the token you just received, until the tool returns DONE. Then answer 'LAST: <the last token you received before DONE>'.",
        check=lambda text, calls: chain_token(25, 24) in text and sum(1 for c in calls if c[0] == "next_token" and not c[2].startswith("ERROR")) == 25, expect="25 valid calls; LAST = token 24"),
    "t8_two_tools_plan": dict(prompt="Compute (12 + 30) * (5 + 9) using the add and multiply tools only. Finish with 'RESULT: <number>'.",
        check=lambda text, calls: bool(re.search(r"RESULT:\s*588\b", text)), expect="588; add, add, multiply"),
}

def post(endpoint, payload, timeout=600):
    request = urllib.request.Request(endpoint + "/chat/completions", data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(request, timeout=timeout))

def run_test(endpoint, model, test, lean_root, max_turns):
    messages = [{"role": "user", "content": test["prompt"]}]
    calls, bad_arguments, completion_tokens, started, final = [], 0, 0, time.time(), ""
    stop_reason = "max_turns"
    for _ in range(max_turns):
        reply = post(endpoint, {"model": model, "messages": messages, "tools": TOOLS, "max_tokens": 1024})
        completion_tokens += reply.get("usage", {}).get("completion_tokens", 0)
        message = reply["choices"][0]["message"]
        tool_calls = message.get("tool_calls") or []
        messages.append({"role": "assistant", "content": message.get("content") or "", **({"tool_calls": tool_calls} if tool_calls else {})})
        if not tool_calls:
            final, stop_reason = message.get("content") or "", reply["choices"][0].get("finish_reason")
            break
        for call in tool_calls:
            name = call["function"]["name"]
            try:
                arguments = json.loads(call["function"]["arguments"] or "{}")
                if not isinstance(arguments, dict):
                    raise ValueError("arguments are not an object")
                result = run_tool(name, arguments, lean_root)
            except ValueError as error:
                arguments, result, bad_arguments = call["function"]["arguments"], f"ERROR: arguments are not valid JSON ({error})", bad_arguments + 1
            calls.append((name, arguments, result[:120]))
            messages.append({"role": "tool", "tool_call_id": call.get("id", ""), "content": result})
    return dict(final=final, calls=calls, bad_arguments=bad_arguments, completion_tokens=completion_tokens,
                seconds=round(time.time() - started, 1), stop_reason=stop_reason, passed=bool(test["check"](final, calls)))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lean", required=True); parser.add_argument("--out", required=True)
    parser.add_argument("--endpoint", default="http://127.0.0.1:8080/v1")
    parser.add_argument("--runs", type=int, default=3); parser.add_argument("--only", default="")
    parser.add_argument("--max-turns", type=int, default=30)
    args = parser.parse_args()
    host = urllib.parse.urlparse(args.endpoint).hostname or ""
    if host != "localhost" and not (ipaddress.ip_address(host).is_loopback or ipaddress.ip_address(host).is_private):
        raise SystemExit("refusing a non-local endpoint")
    model = json.load(urllib.request.urlopen(args.endpoint + "/models", timeout=10))["data"][0]["id"]
    out = pathlib.Path(args.out); out.mkdir(parents=True, exist_ok=True)
    names = [n for n in args.only.split(",") if n] or list(TESTS)
    for name in names:
        results = [run_test(args.endpoint, model, TESTS[name], args.lean, args.max_turns) for _ in range(args.runs)]
        (out / f"{model}__{name}.json").write_text(json.dumps({"expect": TESTS[name]["expect"], "runs": results}, indent=1, default=str))
        passed = sum(r["passed"] for r in results)
        print(f"{name:24s} pass {passed}/{len(results)}  calls {[len(r['calls']) for r in results]}  bad-json {sum(r['bad_arguments'] for r in results)}"
              f"  seconds {[r['seconds'] for r in results]}  stop {[r['stop_reason'] for r in results]}", flush=True)

if __name__ == "__main__":
    main()
