# TheoPhysEval: harness, verifier, and worked problems

Source: *TheoPhysEval: Assessing AI agents at autonomous theoretical physics problem solving*,
arXiv:2610.16733 (local source: `Data/Public/papers/AI/arXiv-2610.16733v1/`). The 50
problems are unpublished research-level "turn-the-crank" calculations with numerical answers.
The official prompts, data files and answer keys are in the gated Hugging Face dataset
`RameshChandra1/TheoPhysEval`.

## What the paper found

- Realistic track (web allowed): GPT-6-Astra in Codex solved 49/50 and Opus 5.5 in Claude
  Code 46/50. Gemini 3.8 Flash solved 34, DeepSeek-V4.1-Flash 27, GLM 5.3 21.
- Isolated track (no web): Astra and Opus 46/50 each. Accuracy barely changes without the
  web, even on problems whose notation is defined only in a paper that was not provided.
  The authors read this as memorization of much of the research literature.
- Answers are graded only as numbers. Symbolic answers were turned into numbers by fixing
  the physical parameters, including irrelevant ones, so no LLM judge is needed.
- The one problem nobody solved, `autonomous-de`, is an inverse problem: find the equation
  behind a plotted curve. The answer is a fractional *delay* equation, a play on
  "differential equation".
- The paper lists Opus's failures:
  - It used initial data from a cited paper instead of the prompt's.
  - It made a wrong convergence claim about an infinite series.
  - It repeatedly hit the output-token limit.
  - It set up a modular differential equation incorrectly.

## Files

| file | role |
|---|---|
| `extract_problems.py` | builds `problems.jsonl` (50 rows: statement, theme, the paper's per-agent grades) from the arXiv source |
| `verifier.py` | numeric grader: parses `1.2e-3`, `1.2\times10^{-3}`, `a/b`, multiples of π; grades the final `\boxed{}` / `KEY=` region; PASS / FAIL / INCONCLUSIVE plus `admitted_unknown`, `confident_wrong` |
| `harness.py` | runs problems through `claude -p` (realistic / isolated tracks, `--max-budget-usd` cap) or an OpenAI-compatible endpoint; dry run unless `--execute` |
| `make_references.py` | runs the solvers and writes `references.local.json` |
| `solutions/maple_leaf.py` | `maple-leaf-dirac-nodes`, `maple-leaf-chern-number`, `maple-leaf-exact-diagonalization` |
| `solutions/theta_logistic.py` | `theta-logistics` |
| `solutions/fractional_delay.py` | `autonomous-de`: forward solve of the published answer, figure digitizer, coefficient fit, identifiability |
| `_test/test_verifier.py` | known-answer tests for the verifier |

Untracked by design: `problems.jsonl`, `references.local.json`, `RESULTS.local.md`,
`results/`, `runs/`, `out/` (see the local `.gitignore`). The benchmark keeps its answers off
the public web to avoid training-data contamination. This repository has a public remote, so
computed answers stay local. Running the solvers regenerates them.

## Worked problems and how each answer was checked

The official keys are gated, so each answer here is **self-checked, not graded**.

- **maple-leaf-dirac-nodes.**
  - Method: a 400² grid search for the minima of ε₄−ε₃, then Nelder–Mead refinement and a
    Brent line search on the Γ–M line.
  - Hermiticity and C₆ invariance of the spectrum are checked to 10⁻¹⁵, with C₆ acting on
    fractional momenta as (u,v) ↦ (−v, u+v), from R b₁ = b₂ and R b₂ = −b₁ + b₂.
  - The nodes form a single C₆ orbit.
  - The Berry phase is a discrete Wilson loop around a node, stable under changes of radius.
  - Symmetry question: the chiral deformation of the Chern problem keeps C₆ but breaks the
    mirrors. It opens a gap linear in φ, so p6 alone does not pin the nodes.
- **maple-leaf-chern-number.**
  - The Chern number comes from the prompt's link-variable formula and is stable from a 24²
    to a 192² mesh. A Kubo-formula integral gives the same magnitude with the opposite sign.
    The prompt's formula measures −∫Ω for the connection A = i⟨n|∂n⟩, which explains the sign.
    The prompt fixes the sign convention through its formula.
  - The gaps come from grid minima refined by two optimizers. All six C₆ images agree to
    10⁻¹⁵.
- **maple-leaf-exact-diagonalization.**
  - Lanczos in the Sᶻ = 0 and Sᶻ = 1 sectors (dimensions 48620 and 43758).
  - Checks: bond counts 18/18/9; coordination (2, 2, 1) at every site; E₀ reproduced in
    the full 2¹⁸ space; every Sᶻ = 1 level reappears in Sᶻ = 0 (SU(2)); the first
    excitation is threefold degenerate.
  - The translation eigenvalue is ⟨ψ₀|T|ψ₀⟩ for the nondegenerate ground state.
- **theta-logistics.**
  - Closed form: the stationary densities are a gamma density (θ = 1) and a generalized
    gamma density (θ = 2), so V₁ = 19 exactly.
  - V₂ is a root of a Laplace transform. The z^(s−1) endpoint singularity silently fooled
    both scipy `quad` and mpmath at 15–20 digits, in the fourth digit. Substituting
    z = w^(1/s) removes the singularity. Three methods now agree.
  - Check: a 20,000-path Euler–Maruyama simulation in ln X gives P_detect = 0.949–0.950.
- **autonomous-de.**
  - The published fractional delay equation is solved by implicit Euler with
    Grünwald–Letnikov weights. Converged against `fig1.png`: 0.012 decades rms, about
    0.7 pixel.
  - The human reading of the figure supplies the delay and constant terms:
    - kinks at t = 1, 2, 3, 4 mean τ = 1;
    - the flat drift on [0, 1) means c₃h₀ + c₄ = 0.
  - With those two facts, least squares recovers c₂ and c₃ within a few percent.
  - c₁ is not identifiable from the picture. The weakest Jacobian direction is almost pure
    c₁, and a 10% change moves the curve by 0.0005 decades, against 0.017 for one pixel.
    The prompt's "round numbers" condition supplies that missing information.

## Running

Use the mathphysics uv venv (numpy, scipy, matplotlib, pillow, mpmath, pytest):

```bash
PY=/home/propdev/.openclaw/workspace/workspace2/repos/mathphysics/.venv/bin/python
cd evals/theophyseval
$PY extract_problems.py                    # problems.jsonl
$PY make_references.py                     # references.local.json (~10 s)
$PY solutions/maple_leaf.py all            # each solver prints its checks
$PY solutions/theta_logistic.py            # ~30 s with the Monte Carlo
$PY solutions/fractional_delay.py --fit --plot out/fractional_delay.png
$PY -m pytest -q _test
$PY harness.py --list                      # paper grades per agent, REF = reference available
$PY harness.py                             # dry run on the referenced problems
$PY harness.py --problems theta-logistics --execute --model sonnet --budget-usd 2
```

`--execute` with the claude backend gives the agent Bash on this machine. The paper ran each
agent in a Docker sandbox with 2 vCPUs, 4 GB and an 8 h cutoff. Here the isolated track only
removes the web tools.

To use the official prompts and keys:
1. Accept the dataset terms on Hugging Face.
2. `huggingface-cli login`.
3. Download the dataset.
4. Pass `--task-dir <bundle>`.

The bundle's `tests/grading_spec.json` format has not been inspected yet, so an adapter from
it to the reference format in `verifier.py` is still to be written.

## Lessons for TMPAI evals

1. **Numeric conversion beats symbolic judging** for agent evals: fix every parameter and
   grade numbers. The registry in `agent/EVAL-AND-TOOLING.md` asks for symbolic ground
   truth. Pairing each symbolic item with a numeric instance gives a judge-free check.
2. **Sign conventions are where confident errors live.** The Chern number has opposite signs
   from the link-variable and Kubo formulas. An agent that "knows" the Kubo sign fails a
   prompt that defines its own convention. `confident_wrong` is the right headline metric.
3. **Silent numerical drift.** Generic quadrature lost the fourth digit on an integrable
   endpoint singularity, and the only signal was a warning. Remove the singularity
   analytically and require two methods to agree.
4. **Prompt over prior.** Opus's realistic-track miss used a cited paper's initial data. The
   harness preamble now says to use the prompt's parameters even when a known reference
   differs.
5. **Inverse problems need priors.** The curve fixes c₂, c₃, c₄ but not c₁. "Round numbers"
   is information, and an agent has to use it deliberately, by searching a lattice of round
   values rather than fitting continuously.

The geometry behind the maple-leaf answers:
- the Berry phase is the holonomy of the Berry connection on the band-3 line bundle around a
  node;
- the Chern number is the first Chern class of the rank-3 occupied bundle over the torus;
- the link variables are lattice parallel transports, as in lattice gauge theory.
