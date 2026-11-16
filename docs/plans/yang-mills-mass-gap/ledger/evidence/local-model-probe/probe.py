#!/usr/bin/env python3
"""Probe a local OpenAI-compatible model with tasks whose answers the program already knows.

Usage: probe.py --n1 N1.txt --n2 N2.txt --lean <LEAN checkout> --out <directory outside git>
               [--endpoint http://127.0.0.1:8080/v1] [--only a,b,...]
Standard library only. Sends prompts to a loopback or private-range endpoint only.
Nothing a model writes is executed. Grading here is by keyword; read the saved answers too.
"""
import argparse, ipaddress, json, pathlib, re, time, urllib.parse, urllib.request

PURITY_SETUP = (
    "Let T be a self-adjoint trace-class operator with eigenvalues l_1 >= l_2 >= ... >= 0, l_1 > 0. "
    "Fix n >= 1. Put r_j = (l_j/l_1)^n for j >= 2, s = sum_{j>=2} r_j, q = sum_{j>=2} r_j^2, "
    "Z(n) = sum_j l_j^n = l_1^n (1+s), Z(2n) = l_1^{2n} (1+q), and P = Z(2n)/Z(n)^2 = (1+q)/(1+s)^2.\n"
)
PROOF_CORRECT = PURITY_SETUP + (
    "Claim: P <= 1.\nProof. Step 1: 0 <= r_j <= 1 for every j >= 2. Step 2: hence r_j^2 <= r_j, so q <= s. "
    "Step 3: therefore P = (1+q)/(1+s)^2 <= (1+s)/(1+s)^2 = 1/(1+s). Step 4: since s >= 0, 1/(1+s) <= 1. QED."
)
PROOF_WRONG = PURITY_SETUP + (
    "Claim: P >= 1.\nProof. Step 1: 0 <= r_j <= 1 for every j >= 2. Step 2: hence r_j^2 >= r_j, so q >= s. "
    "Step 3: therefore P = (1+q)/(1+s)^2 >= (1+s)/(1+s)^2 = 1/(1+s). Step 4: since s >= 0, 1/(1+s) >= 1. QED."
)

TELESCOPE_SETUP = (
    "Let Z(a,b,c,d) > 0 be symmetric in its four arguments, and define\n"
    "P1 = Z(2n,n,n,n)/Z(n,n,n,n)^2, P2 = Z(2n,2n,n,n)/Z(2n,n,n,n)^2, "
    "P3 = Z(2n,2n,2n,n)/Z(2n,2n,n,n)^2, P4 = Z(2n,2n,2n,2n)/Z(2n,2n,2n,n)^2.\n"
)
TELESCOPE_WRONG = TELESCOPE_SETUP + (
    "Claim: P1^4 * P2^4 * P3^2 * P4 = Z(2n,2n,2n,2n) / Z(n,n,n,n)^16.\n"
    "Proof. Step 1: P1^4 = Z(2n,n,n,n)^4 / Z(n,n,n,n)^8. Step 2: P2^4 = Z(2n,2n,n,n)^4 / Z(2n,n,n,n)^8. "
    "Step 3: P3^2 = Z(2n,2n,2n,n)^2 / Z(2n,2n,n,n)^4. Step 4: multiplying Steps 1-3 with P4, every intermediate "
    "factor cancels and the product is Z(2n,2n,2n,2n) / Z(n,n,n,n)^16. QED."
)
EXPONENTIAL_WRONG = (
    "Claim: if 0 <= d <= 1/4 then log(1/(e^d - 1)) >= log(1/d).\n"
    "Proof. Step 1: for 0 <= x <= 1 the function e^x is convex. Step 2: hence e^x - 1 <= x for 0 <= x <= 1. "
    "Step 3: so e^d - 1 <= d, and 1/(e^d - 1) >= 1/d. Step 4: log is increasing, which gives the claim. QED."
)
TRACE_WRONG = (
    "Let A, B be N x N complex matrices and consider the linear map L(Z) = A Z B on the space of N x N matrices, "
    "with basis the matrix units E_ij.\nClaim: tr(L) = tr(AB).\n"
    "Proof. Step 1: tr(L) = sum_{i,j} (A E_ij B)_{ij}. Step 2: (A E_ij B)_{ij} = A_ii B_jj. "
    "Step 3: so tr(L) = sum_{i,j} A_ii B_jj = sum_i (AB)_ii = tr(AB). QED."
)
LAYER_CORRECT = (
    "Let nu be a finite measure on (0, infinity), N(u) = nu((0,u]), and s > 0.\n"
    "Claim: integral of e^{-sE} d nu(E) = integral_0^infinity s e^{-su} N(u) du.\n"
    "Proof. Step 1: e^{-sE} = integral_E^infinity s e^{-su} du. Step 2: substitute and write the inner integral as "
    "integral_0^infinity 1{E <= u} s e^{-su} du. Step 3: all integrands are nonnegative, so by Tonelli the order of "
    "integration may be exchanged. Step 4: integral 1{E <= u} d nu(E) = nu((0,u]) = N(u). QED."
)

REVIEW_QUESTION = (
    "\n\nYou are a referee. Check each step. Answer in exactly this form:\n"
    "VERDICT: CORRECT or INCORRECT\nFIRST WRONG STEP: <number or none>\nREASON: <one sentence>"
)

def build_tasks(args):
    n1 = pathlib.Path(args.n1).read_text(errors="replace")
    n2_lines = pathlib.Path(args.n2).read_text(errors="replace").splitlines()
    n2_excerpt = "\n".join(n2_lines[150:245])
    lean = pathlib.Path(args.lean, "ChallengeSpec.lean").read_text().splitlines()
    joint = "\n".join(lean[185:188])
    return {
        "a_extract_short": dict(
            prompt="Below is an excerpt from the introduction of a mathematics paper (plain text extracted from a PDF).\n\n"
            + n2_excerpt + "\n\nAnswer from the excerpt only, in three numbered lines:\n"
            "1. The definition of the quantity Delta_beta(n) in equation (1.6).\n"
            "2. The recurrence inequality for Delta_beta(2n) that Section 2 proves.\n"
            "3. The threshold eta_* and the lower bound on m_lat(beta) that follows.",
            keys=[r"4\s*log", r"64", r"2\s*\^?\s*\{?\s*[-−]\s*10|1/1024", r"log\s*2"],
            expected="Delta = 4 log Z(n,n) - log Z(2n,2n); Delta(2n) <= 64 Delta(n)^2; eta_* = 2^-10 gives m_lat >= log 2 / n"),
        "b_extract_long": dict(
            prompt="Below is the full text of a mathematics paper (plain text extracted from a PDF).\n\n" + n1
            + "\n\nAnswer from the paper only, in two numbered lines:\n"
            "1. For which values of n does the main theorem on exponential decay hold, and at which temperatures?\n"
            "2. What lower bound does the paper prove for the ratio of the twisted to the untwisted partition function on an annulus of inner radius N?",
            keys=[r"n\s*(≥|>=|\\ge|\\geq)\s*3|three or more|at least 3", r"log\s*N"],
            expected="n >= 3, every temperature (all beta); ratio >= 1 - O(1/log N)"),
        "c_lean_quantifiers": dict(
            prompt="Here is a Lean 4 definition.\n\n" + joint + "\n\nRead only the quantifier structure. Answer in exactly three lines:\n"
            "C: chosen BEFORE or AFTER the sequence p?\nq: chosen BEFORE or AFTER p?\nT: may it depend on q? YES or NO",
            keys=[r"C:.*BEFORE", r"q:.*AFTER", r"T:.*YES"],
            expected="C before p; q after p; T may depend on q"),
        "d_referee_wrong_proof": dict(prompt=PROOF_WRONG + REVIEW_QUESTION,
            keys=[r"VERDICT:\s*INCORRECT", r"FIRST WRONG STEP:\s*2"], expected="INCORRECT, step 2 (r^2 <= r on [0,1])"),
        "e_referee_correct_proof": dict(prompt=PROOF_CORRECT + REVIEW_QUESTION,
            keys=[r"VERDICT:\s*CORRECT", r"FIRST WRONG STEP:\s*none"], expected="CORRECT, none"),
        "h_referee_telescope_wrong": dict(prompt=TELESCOPE_WRONG + REVIEW_QUESTION,
            keys=[r"VERDICT:\s*INCORRECT", r"FIRST WRONG STEP:\s*(1|4)"], expected="INCORRECT: P1 must be raised to the 8th power (step 1 gives Z(2n,n,n,n)^4, step 2 needs ^8 to cancel); accept step 1 or 4"),
        "i_referee_exponential_wrong": dict(prompt=EXPONENTIAL_WRONG + REVIEW_QUESTION,
            keys=[r"VERDICT:\s*INCORRECT", r"FIRST WRONG STEP:\s*2"], expected="INCORRECT, step 2 (e^x - 1 >= x; convexity gives e^x - 1 <= (e-1)x)"),
        "j_referee_trace_wrong": dict(prompt=TRACE_WRONG + REVIEW_QUESTION,
            keys=[r"VERDICT:\s*INCORRECT", r"FIRST WRONG STEP:\s*3"], expected="INCORRECT, step 3 (sum_{i,j} A_ii B_jj = tr(A) tr(B))"),
        "k_referee_layer_correct": dict(prompt=LAYER_CORRECT + REVIEW_QUESTION,
            keys=[r"VERDICT:\s*CORRECT", r"FIRST WRONG STEP:\s*none"], expected="CORRECT, none"),
        "f_compute": dict(
            prompt="For su(N) with inner product <X,Y> = -(1/N) Re tr(XY), the Killing form is K(X,Y) = 2N tr(XY) and the Casimir constant C is defined by K = -C <.,.>.\n"
            "1. Compute C for N = 4.\n2. With b = 11 C / (24 pi^2), give b for N = 2 as a simplified fraction times 1/pi^2.\nShow one line of work each, then 'ANSWERS: C=<value>, b=<value>'.",
            keys=[r"C\s*=\s*32", r"11\s*/\s*\(?3|\\frac\{11\}\{3"], expected="C = 32; b = 11/(3 pi^2)"),
        "g_rust_test": dict(
            prompt="Rust. A function `sphere_monomial_average(exponents: &[u32]) -> Option<BigRational>` returns the average of n_1^e_1 ... n_d^e_d over the unit sphere of R^d, d = exponents.len(). "
            "Formula: zero if some exponent is odd; otherwise with e_i = 2 a_i and k = sum a_i, the value is prod_i (2 a_i - 1)!! / (d (d+2) ... (d+2k-2)). "
            "A helper `q(numerator: i64, denominator: i64) -> BigRational` exists.\n"
            "Write exactly one #[test] function that asserts the value for exponents [2, 2, 2, 0, 0, 0]. Compute the expected value by hand in a comment. Output only Rust code.",
            keys=[r"480", r"#\[test\]", r"\[\s*2\s*,\s*2\s*,\s*2\s*,\s*0\s*,\s*0\s*,\s*0\s*\]"], expected="q(1, 480): 1!!^3 / (6*8*10)"),
    }

def check_endpoint(endpoint):
    host = urllib.parse.urlparse(endpoint).hostname or ""
    if host == "localhost":
        return
    address = ipaddress.ip_address(host)
    if not (address.is_loopback or address.is_private):
        raise SystemExit("refusing a non-local endpoint")

def main():
    parser = argparse.ArgumentParser()
    for name in ("n1", "n2", "lean", "out"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--endpoint", default="http://127.0.0.1:8080/v1")
    parser.add_argument("--only", default="")
    parser.add_argument("--max-tokens", type=int, default=1500)
    args = parser.parse_args()
    check_endpoint(args.endpoint)
    out = pathlib.Path(args.out); out.mkdir(parents=True, exist_ok=True)
    models = json.load(urllib.request.urlopen(args.endpoint + "/models", timeout=10))
    model = models["data"][0]["id"]
    tasks = build_tasks(args)
    wanted = [name for name in args.only.split(",") if name] or list(tasks)
    summary = []
    for name in wanted:
        task = tasks[name]
        body = json.dumps({"model": model, "temperature": 0.2, "max_tokens": args.max_tokens,
                           "messages": [{"role": "user", "content": task["prompt"]}]}).encode()
        request = urllib.request.Request(args.endpoint + "/chat/completions", data=body, headers={"Content-Type": "application/json"})
        started = time.time()
        try:
            reply = json.load(urllib.request.urlopen(request, timeout=1800))
            message = reply["choices"][0]["message"]
            content = message.get("content") or ""
            usage = reply.get("usage", {})
            finish = reply["choices"][0].get("finish_reason")
        except Exception as error:  # record and continue
            content, usage, finish, message = f"ERROR: {error}", {}, "error", {}
        seconds = round(time.time() - started, 1)
        hits = [bool(re.search(key, content, re.I | re.S)) for key in task["keys"]]
        (out / f"{model}__{name}.json").write_text(json.dumps({"task": name, "model": model, "expected": task["expected"],
            "content": content, "reasoning": message.get("reasoning_content"), "usage": usage, "finish": finish,
            "seconds": seconds, "keyword_hits": hits}, indent=1))
        summary.append((name, sum(hits), len(hits), usage.get("prompt_tokens"), usage.get("completion_tokens"), seconds, finish))
        print(f"{name:26s} keys {sum(hits)}/{len(hits)}  prompt {usage.get('prompt_tokens')}  completion {usage.get('completion_tokens')}  {seconds}s  {finish}", flush=True)
    (out / f"{model}__summary.json").write_text(json.dumps(summary, indent=1))

if __name__ == "__main__":
    main()
