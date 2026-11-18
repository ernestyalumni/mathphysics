"""Offline tests for bin/ym-micro. No model server is needed.

Run: python3 -I -m unittest docs/plans/yang-mills-mass-gap/tests/test_ym_micro.py
"""
import contextlib, importlib.machinery, importlib.util, io, json, pathlib, re, tempfile, unittest

PLAN = pathlib.Path(__file__).resolve().parents[1]


def load(name, path):
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


core = load("ym_micro", PLAN / "bin" / "ym-micro")

PROOF = r"""
\begin{proposition}\label{p:one}
If $0\le r\le1$ then Eq.~\eqref{e:target} holds.
\end{proposition}
\begin{proof}
Since $0\le r\le1$, $r^2\le r$. Next,
\begin{equation}\label{e:mid}
 q=\sum_j r_j^2\le\sum_j r_j=s .
\end{equation}
Hence $P\le1$ by Eq.~\eqref{e:target}. This uses the value $0.5$ only once.
\end{proof}
\begin{equation}\label{e:target}
 P=\frac{1+q}{(1+s)^2}
\end{equation}
"""


class Steps(unittest.TestCase):
    def test_display_stays_in_its_sentence(self):
        steps = core.split_steps(r"Since $a\le b$, $c\le d$. Next, \begin{equation} x=y . \end{equation} Hence $z$ follows from it. A decimal $0.5$ stays.")
        self.assertEqual(len(steps), 4)
        self.assertIn(r"\begin{equation}", steps[1])
        self.assertTrue(steps[1].startswith("Next,"))
        self.assertIn("0.5", steps[3])

    def test_abbreviation_does_not_split(self):
        steps = core.split_steps(r"By Prop. \ref{p:one} the bound holds. Then we are done here.")
        self.assertEqual(len(steps), 2)

    def test_proof_units_carry_statement_and_references(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory, "t.tex")
            path.write_text(PROOF)
            units = core.proof_units(path, "t.tex", "t")
        self.assertEqual(len(units), 1)
        unit = units[0]
        self.assertIn("p-one", unit["name"])
        self.assertIn(r"\begin{proposition}", unit["context"])
        self.assertTrue(any("e:target" in block for block in unit["refs"]))
        self.assertEqual(len(unit["steps"]), 4)

    def test_prompt_uses_the_probe_question(self):
        probe = (PLAN / "ledger/evidence/local-model-probe/probe.py").read_text()
        self.assertIn("VERDICT: CORRECT or INCORRECT", probe)
        self.assertIn("FIRST WRONG STEP: <number or none>", core.REVIEW_QUESTION)
        unit = {"context": "S", "refs": [], "steps": ["a", "b"]}
        prompt = core.referee_prompt(unit, (1, 2))
        self.assertIn("Step 2: b", prompt)
        self.assertIn("Check only steps 1 to 2", prompt)

    def test_windows(self):
        self.assertEqual(core.windows_of(14), [None])
        self.assertEqual(core.windows_of(25), [(1, 10), (11, 20), (21, 25)])


class Mutants(unittest.TestCase):
    def test_one_change_at_a_known_step(self):
        unit = {"name": "u", "truth": "unknown", "steps": [r"$a\le b$ and nothing else.", r"so $x^{8}=y$ follows.", r"The \left( bracket \right) is not an inequality."]}
        found = core.mutants(unit, limit=6)
        self.assertEqual(len(found), 2)
        for mutant in found:
            changed = [number for number, (old, new) in enumerate(zip(unit["steps"], mutant["steps"]), 1) if old != new]
            self.assertEqual(changed, [mutant["mutated_step"]])
            self.assertEqual(mutant["truth"], "mutant")
        self.assertTrue(any(r"\ge" in mutant["steps"][0] for mutant in found))
        self.assertTrue(all(mutant["steps"][2] == unit["steps"][2] for mutant in found))

    def test_no_exponent_change_on_a_power_compared_with_zero(self):
        # ||E psi||^2 = 0 and ||E psi||^3 = 0 are the same statement: not a planted error.
        harmless = [r"Hence $\|E_H((0,\Delta'])\psi\|^2=0$ on the dense set.", r"so $\left|f\right|^{2}\ge 0$.",
                    r"hence $0 < \|v\|^{4}$.", r"and $(a-b)^2 = 0$."]
        for step in harmless:
            found = core.mutants({"name": "u", "truth": "unknown", "steps": [step]})
            self.assertFalse([m for m in found if m["mutation"].startswith("exponent")], step)
        kept = [r"$\|v\|^2 = 2\|w\|$.", r"$\Delta_0 = c\,x^2$.", r"$x^{2} \le 10$.", r"and $(\lambda-\mu)^2>0$.",
                r"$0=\bigl\langle H_2(g,g)+X((df)g)^2\bigr\rangle$.", r"$1 - x^2 = 0$."]
        for step in kept:
            self.assertTrue(any(m["mutation"].startswith("exponent") for m in
                                core.mutants({"name": "u", "truth": "unknown", "steps": [step]})), step)


class Verdicts(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(core.parse_verdict("VERDICT: INCORRECT\nFIRST WRONG STEP: 3\nREASON: sign.")["step"], "3")
        self.assertEqual(core.parse_verdict("**VERDICT:** CORRECT\nFIRST WRONG STEP: none\nREASON: fine")["verdict"], "CORRECT")
        self.assertEqual(core.parse_verdict("I think it is fine.")["verdict"], "UNPARSED")


class LeanReading(unittest.TestCase):
    SNIPPET = "/-- The sites. -/\ndef Site (n : Nat) : Type := Fin n -> Fin 4\ntheorem site_card (n : Nat) : Fintype.card (Site n) = n ^ 4 := by"

    def good(self, text):
        accepted, why, fields = core.check_lean_reading(text, self.SNIPPET)
        self.assertTrue(accepted, why)
        return fields

    def bad(self, text):
        self.assertFalse(core.check_lean_reading(text, self.SNIPPET)[0])

    def test_three_good(self):
        self.assertEqual(self.good("KIND: definition\nSTATEMENT: Site n is the type of functions.\nIDENTIFIERS: Site, Fin")["kind"], "definition")
        self.good("KIND: theorem\nSTATEMENT: The number of sites is n to the fourth.\nIDENTIFIERS: site_card, Fintype.card, Site")
        self.good("**KIND:** theorem\n**STATEMENT:** Counts sites.\n**IDENTIFIERS:** `site_card`, `Site`")

    def test_three_bad(self):
        self.bad("KIND: theorem\nSTATEMENT: Counts sites.")  # header missing
        self.bad("KIND: theorem\nSTATEMENT: Counts sites.\nIDENTIFIERS: site_card, Lattice.volume")  # invented identifier
        self.bad("KIND: theorem\nSTATEMENT: " + "word " * 90 + "\nIDENTIFIERS: Site")  # too long


class LeanSnippet(unittest.TestCase):
    def test_starts_at_target_with_its_doc_comment_only(self):
        source = "\n".join([
            "theorem helper (h : a = b) : b = a := by",  # 1
            "  rw [h]",                                  # 2
            "",                                          # 3
            "/-- The answer. -/",                        # 4
            "@[simp]",                                   # 5
            "theorem answer (x : Nat) : x + 0 = x := by",  # 6
            "  simp",                                    # 7
        ])
        with tempfile.TemporaryDirectory() as directory:
            pathlib.Path(directory, "A.lean").write_text(source)
            got = core.lean_snippet(directory, "https://x/blob/abc123/A.lean#L6")
            plain = core.lean_snippet(directory, "https://x/blob/abc123/A.lean#L1")
        self.assertEqual(got["text"].splitlines(), ["/-- The answer. -/", "@[simp]", "theorem answer (x : Nat) : x + 0 = x := by"])
        self.assertNotIn("helper", got["text"])
        self.assertEqual(plain["text"], "theorem helper (h : a = b) : b = a := by")


class Quotes(unittest.TestCase):
    def test_only_verbatim_quotes_get_a_locator(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory, "s.tex")
            path.write_text("First line of text.\n\\begin{equation}\n \\Delta_\\beta(n)=4\\log Z(n,n)-\\log Z(2n,2n)\n\\end{equation}\nLast line.\n")
            text, squashed, origin = core.load_source([path])
        answer = ("ITEM: the defect\nQUOTE: \\Delta_\\beta(n) = 4\\log Z(n,n) - \\log Z(2n,2n)\n"
                  "ITEM: invented\nQUOTE: \\Delta_\\beta(n) = 16\\log Z(n,n,n,n) - something else entirely")
        items = core.check_quotes(answer, squashed, origin)
        self.assertEqual([item["verbatim"] for item in items], [True, False])
        self.assertEqual((items[0]["file"], items[0]["line"]), ("s.tex", 3))
        self.assertNotIn("file", items[1])


class Fragments(unittest.TestCase):
    GOOD = "\\subsection{A title}\n" + "Some text with $x^2$ in it. " * 60 + "\n\\begin{equation}\n a=b\n\\end{equation}\n"

    def test_checks(self):
        self.assertEqual(core.check_fragment(self.GOOD), "")
        self.assertIn("forbidden", core.check_fragment(self.GOOD + "\\input{/etc/passwd}"))
        self.assertIn("forbidden", core.check_fragment(self.GOOD + "\\immediate\\write18{ls}"))
        self.assertIn("unbalanced", core.check_fragment(self.GOOD + "\\begin{proof} x"))
        self.assertIn("short", core.check_fragment("\\subsection{A}\n short"))
        self.assertIn("sectioning", core.check_fragment("Here is your fragment. " * 80))

    def test_extract_from_fence(self):
        self.assertTrue(core.extract_fragment("Sure:\n```latex\n\\subsection{T}\nbody\n```\nDone.").startswith("\\subsection{T}"))


class Endpoints(unittest.TestCase):
    def test_refuses_public_hosts(self):
        for endpoint in ("http://8.8.8.8/v1", "https://api.example.com/v1"):
            with self.assertRaises(SystemExit):
                core.check_endpoint(endpoint)
        core.check_endpoint("http://127.0.0.1:8080/v1")
        core.check_endpoint("http://10.0.0.7:8080/v1")


class DryRun(unittest.TestCase):
    def test_queue_runs_without_a_server_and_resumes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            unit = {"name": "u1", "source": "t.tex:1", "truth": "unknown", "context": "", "refs": [], "steps": ["a.", "b."]}
            (root / "u1.json").write_text(json.dumps(unit))
            queue = root / "q.jsonl"
            queue.write_text(json.dumps({"id": "u1", "type": "referee", "unit": str(root / "u1.json"), "seeds": [0, 1], "repeatable": False}) + "\n")
            arguments = type("A", (), dict(endpoint="http://127.0.0.1:9/v1", out=str(root / "out"), queue=str(queue), hours=0.01,
                                           model=None, dry_run=True, timeout=5, max_passes=2, down_retries=1))
            self.assertEqual(core.run_queue(arguments), 0)
            result = json.loads((root / "out/referee/u1.p0.json").read_text())
            self.assertEqual(result["runs"], 2)
            self.assertEqual(result["flags"], 0)
            self.assertIn("unreviewed", result["label"])
            before = (root / "out/referee/u1.p0.json").stat().st_mtime_ns
            self.assertEqual(core.run_queue(arguments), 0)
            self.assertEqual((root / "out/referee/u1.p0.json").stat().st_mtime_ns, before)


class ExportCases(unittest.TestCase):
    def test_split_is_by_unit_and_mutants_stay_with_their_proof(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            (root / "state/units").mkdir(parents=True)
            for k in range(4):
                unit = {"name": "u%d" % k, "truth": "unknown", "context": "S", "refs": [],
                        "steps": [r"$a\le b$ here.", r"so $x^{2}=y$.", "done."]}
                (root / "state/units" / ("u%d.json" % k)).write_text(json.dumps(unit))
            (root / "state/units/long.json").write_text(json.dumps(
                {"name": "long", "truth": "unknown", "steps": ["s."] * 30}))
            (root / "state/units/u0--mut1.json").write_text(json.dumps({"name": "u0--mut1", "truth": "mutant", "steps": ["x"]}))
            arguments = type("A", (), dict(state=str(root / "state"), out=str(root / "out"), max_steps=25, mutants=4))
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(core.export_cases(arguments), 0)
            parts = {name: json.loads((root / "out" / ("referee-%s.json" % name)).read_text()) for name in ("discovery", "test")}
        groups = {name: {case["group"] for case in cases} for name, cases in parts.items()}
        self.assertFalse(groups["discovery"] & groups["test"])
        self.assertEqual(groups["discovery"] | groups["test"], {"u0", "u1", "u2", "u3"})
        for cases in parts.values():
            for case in cases:
                if case["expected"] == "INCORRECT":
                    self.assertTrue(case["id"].startswith(case["group"] + "--mut"))
                    self.assertIn("Step %d:" % case["step"], case["prompt"])
                else:
                    self.assertIsNone(case["step"])
                self.assertNotIn("#", case["id"])
        self.assertEqual(sum(c["expected"] == "INCORRECT" for cs in parts.values() for c in cs), 8)


if __name__ == "__main__":
    unittest.main()
