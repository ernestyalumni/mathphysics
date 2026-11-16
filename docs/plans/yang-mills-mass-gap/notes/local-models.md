# Local models for this program: what Ernest has, what fits, what each can do

Author: strategist (Fable 5.2), 2026-11-15, session 3. **Working note, written incrementally.** If the last section is "Probe: pending", the session was cut off before measurements; resume at "Resume here".

## 1. Inventory (read from `claw-work/Deployments`, not from memory)

**Desktop (canonical host).** Xeon E5-1650 v3 (12 threads), 31 GB RAM, RTX 3060 12 GB (the llama.cpp server's GPU), GTX 980 Ti 6 GB. Docker llama.cpp server, profiles in `claw-work/Deployments/Scripts/llama-cpp-server/profiles/`.

| Profile | Weights on disk | Context in profile | VRAM measured | tok/s | Workflow checks (their BENCHMARKS.md) |
|---|---|---:|---:|---:|---|
| `qwythos-9b-q5` (Qwythos 9B, "Claude Mythos 5" distill) | yes | **65,536** | 6.2 GB at 32k | 40.9 | **6/6** (their recommended model) |
| `qwythos-9b-q6` | yes | 32,768 | 6.9 GB | 35.7 | 6/6 |
| `qwen38-9b-distill-q5` / `q4` | on the second (NVMe) model volume per profile; not verified today | 32,768 | 6.2 / 5.5 GB | 41 / 45 | 4/6, 5/6 |
| `qwen35-9b-distilled` Q4 / Q5 / Q8 | yes | 32,768 / 32,768 / 16,384 | 5.5 / 6.2 / 8.8 GB | 45 / 41 / 31 | 3/6, 4/6, 5/6 |
| `gemma-4-e4b-it` Q5 / Q6 | not verified today | 32,768 | — | — | not in the table |

Also on disk, not in GGUF form: Qwen3.5-2B, Nemotron-3-Nano-4B, Kearm/Meta8byo (safetensors).

**Mac (M4 Pro, 64 GB unified memory, per `llama-cpp-server-macos/ADDING-A-MODEL.md`).** Native Metal server. Profiles for Qwen3.5 9B and Qwen3.8 9B at Q4–Q8 with context **262,144** (and one at 1,048,576 split over 4 slots). Measured there for Qwen3.8-9B Q8: decode 24 tok/s at short context, 16 tok/s at 60k, 11 tok/s at 131k; prefill 300 → 158 tok/s, so a 131k-token prompt takes about 14 minutes before the first output token.

## 2. Sizes of the things we want read (characters / 4 ≈ tokens)

| Source | Size | Approx. tokens | Fits in 32k | Fits in 65k | Fits in 262k (Mac) |
|---|---:|---:|---|---|---|
| N1 (2D O(n) decay, 20 pp.) | 74 KB | 19k | yes | yes | yes |
| N2 (O(4) mass bounds, 102 pp.) | 390 KB | 100k | by section | by section | **whole** |
| Assessment (review group) | 116 KB | 29k | barely | yes | yes |
| `ARM/statement.md` | 132 KB | 34k | no | yes | yes |
| `ARM/summary.md` | 74 KB | 19k | yes | yes | yes |
| One Lean declaration with context | 1–3 KB | under 1k | yes | yes | yes |

**So context is not the binding constraint.** With `qwythos-9b-q5` at 65k on the desktop, or any 9B at 262k on the Mac, every source fits. The question is whether a 9B model's reading is *correct*. That needs measurement, not opinion: section 4.

## 3. Correction to what the strategist said earlier

Session 2 said "a 9B model with 32k context cannot do the reading packets". That overstated the context problem (the 65k and 262k profiles exist) and asserted the capability limit without testing it. The claim to test is narrower: *a 9B model cannot be trusted, unreviewed, on packets whose output is a mathematical digest.* Section 4 tests it.

## 3a. GPU containers on the desktop: toolkit missing, workaround in place (2026-11-16)

`./run-recommended.sh` fails with `could not select device driver "" with capabilities: [[gpu]]`. Cause, from `/var/log/dpkg.log`: `nvidia-container-toolkit` was **purged on 2026-11-11 19:16**; `apt-cache policy` shows no candidate, so the NVIDIA apt repository is not configured either. `daemon.json` still names `nvidia-container-runtime`, which no longer exists. The host driver is 580.178.04 and `nvidia-smi` works. Reinstalling needs sudo; the strategist has none and did not try.

**Workaround, no sudo, no system change:** `bin/ym-local-llm` starts the same llama.cpp image with the GPU device files and five host driver libraries bind-mounted by hand, pinned to the RTX 3060 by UUID. Verified 2026-11-16: Qwythos 9B Q5 with 65,536 context loads in 16 s, 6.6 GB on the card, 47 tok/s generation, 870 tok/s prompt processing; the whole N1 paper (17k tokens) is read and answered in 22 s.

    bin/ym-local-llm start empero-ai/Qwythos-9B-Claude-Mythos-5-1M-GGUF/Qwythos-9B-Claude-Mythos-5-1M-Q5_K_M.gguf qwythos-9b-q5-k-m 65536
    bin/ym-local-llm status
    bin/ym-local-llm stop

It serves `http://127.0.0.1:8081/v1` (not 8080, so it never collides with claw-work's own launcher). Other GPU containers (Isaac, CuLLM images, …) are still broken until Ernest reinstalls the toolkit.

## 4. Probe results: Qwythos 9B Q5_K_M (claw-work's benchmark winner)

Four runs per task (one on CPU, three on the GPU), temperature 0.2, reasoning budget 512, tasks whose answers the strategist knows. Script: `ledger/evidence/local-model-probe/probe.py`. Raw answers are outside git. Small numbers: read the pattern, not the percentages.

| Task | What it tests | Right / 4 | What went wrong when it did |
|---|---|---:|---|
| a. Extract Eq. (1.6) and two inequalities from a 1,500-token excerpt of N2 | Verbatim extraction | 3 | Once the last line was cut short; never a wrong formula |
| b. Whole paper N1 in context (17k tokens), two questions | Long-context retrieval | **4** | — ("n ≥ 3, every β"; "Z(π)/Z(0) ≥ 1 − C_β/log N") |
| c. Quantifier order of `JointExistence` (3 lines of Lean) | Reading formal syntax | 3 | Once answered off-template |
| d. Referee: planted error r² ≥ r on [0,1] | Easy one-line error | **4** | — |
| i. Referee: "convexity gives e^x − 1 ≤ x" | Wrong inequality | **4** | — |
| e, k. Referee two correct proofs | False alarms | **4, 4** | — |
| j. Referee: "tr(Z ↦ AZB) = tr(AB)", error in step 3 | Wrong identity | **1** | Twice blamed a correct step; once said CORRECT |
| h. Referee: telescoping product with P₁⁴ instead of P₁⁸ | Bookkeeping across four steps | **0** | Said CORRECT every time |
| f. C for su(4), b_* for su(2) from given formulas | Two-step arithmetic | 1 | C = 4; b = 11/48π², 11/24π², 11/16π² |
| g. One Rust test, expected value computed from a given formula | Applying a formula | 1 | Denominators 1,451,520; 5,760; 112,400,064,000 instead of 480 |

**Reading.**

- It copies accurately from a short excerpt. It follows a fixed answer template most of the time, not always.
- As a referee it reliably caught the two errors that are visible in a single line, never caught the one that needs tracking exponents through several lines, and was mostly wrong on the one that needs an identity it has to recall. It raised no false alarm on two correct proofs (8 of 8).
- It cannot be trusted to compute: the two tasks that required applying a formula were right once in four tries each.
- It reads a whole 20-page paper and answers a factual question about it correctly, in 22 seconds on the GPU. Context and speed are not the problem; judgment is.
- So: **an INCORRECT verdict from this model is a cheap flag worth a stronger model's minute. A CORRECT verdict carries no weight at all.** Its own calculations are never evidence.

## 6. What to give which local model

| Use | Model and where | Why | Guard |
|---|---|---|---|
| **Flagging referee.** Feed each numbered step group of our proofs and of every digest, several seeds each; collect INCORRECT flags for a T2+ reviewer | `qwythos-9b-q5` on the desktop GPU (after the fix in 3a) | Probe rows d, e, i, k; costs nothing; runs unattended | Only flags are used. Silence means nothing (row h). |
| **Verbatim indexing.** For one excerpt at a time: list every numbered equation, definition and theorem statement exactly as printed, with its number | same; or any 9B at 262k context on the Mac for whole papers | Probe row a | Accept a line only if it occurs verbatim in the source (string match by script) |
| **Lean signature cards** (packet YM-T14): kind, identifiers, one-sentence reading of one declaration | same | Row c shows it reads quantifiers but drifts from the format | Mechanical checks in YM-T14; label `unreviewed` |
| **Rust tests with the expected value supplied by the packet** | same, through Pi or Hermes | Row g: it cannot derive the value | `cargo test` decides; never let it choose the oracle |
| Reruns, formatting, receipts, watch tasks | any 9B | No judgment needed | — |
| **Not for any local 9B:** digests of papers (YM-R01–R05), derivations, statement red-teaming (YM-A03), anything in workstream C, any review whose "accept" would close a task | — | Rows f, g, h, j | — |

Which of the existing profiles:

1. **`qwythos-9b-q5`** is the one to use on the desktop: best in claw-work's own benchmark (6/6), 65k context, 6.2 GB, about 41 tok/s on the GPU.
2. `qwen35-9b-distilled-q8` or `qwen38-9b-distill-q5` as a **second voice** for the flagging referee (a flag raised by two different fine-tunes is worth more). They are fine-tunes of the same 9B base family, so this is weak independence. `gemma-4-e4b-it` is a different family and would be better for that, but it is smaller and unmeasured here.
3. Not worth using here: Q8 of a model already run at Q5 (same parameters; claw-work measured no gain), the 4B Nemotron (failed claw-work's tool checks), the 2B Qwen.

What would actually raise the ceiling: a **27B-class model on the Mac** (section 5). Before giving any new model a packet, run the admission test:

    python3 -I ledger/evidence/local-model-probe/probe.py --n1 <N1 text> --n2 <N2 text> --lean <LEAN checkout> --out <dir outside git> --endpoint http://<host>:<port>/v1

(`pdftotext` of the two preprints named in `SOURCES.md` gives the text files.) A model that gets rows f, g, h and j right is a candidate for T2 review work; one that does not stays at T0/T1.


## 5. Larger models that could run locally (web search 2026-11-16; **secondary sources, unverified**; check the Hugging Face model cards before any download)

What the search says exists: Qwen3.5 comes as 27B dense, 35B-A3B (mixture of experts, about 3B parameters active per token) and larger; a Qwen3.6-35B-A3B exists; the official Qwen3.8 release is a **27B**. The "Qwen3.8-9B-Distill" and "Qwythos-9B" that Ernest has are community distillations onto a 9B model, not the official Qwen3.8. Reported GGUF sizes for Qwen3.8-27B: about 16.5 GB at Q4_K_M, 29 GB at Q8_0; a 3.7-bit community build is 12.6 GB. Qwen3.5-35B-A3B at Q4_K_M: about 21 GB.

| Machine | Budget | Candidate | Expected behaviour (estimates from the bandwidth rule in ADDING-A-MODEL.md; not measured) |
|---|---|---|---|
| Mac M4 Pro, 64 GB | about 38 GB (60% of RAM) for weights + KV + buffers | **Qwen3.8-27B Q6 or Q8** (22–29 GB), context 32k–64k | Decode ceiling ≈ 273 GB/s ÷ file size: about 9 tok/s at Q8, 12 at Q6, 16 at Q4. Three times the parameters of the 9B models. This is the first thing to try. |
| Mac M4 Pro, 64 GB | same | Qwen3.5/3.6-35B-A3B Q4–Q6 | Much faster decode (few active parameters), weaker than a dense 27B per the usual pattern. Good for bulk T1 work. |
| Desktop, RTX 3060 12 GB | 12 GB VRAM | 9B at Q8, or a 12–14B at Q4 with 8k–16k context | Full GPU. A step up from 9B Q5 only in quantization or modestly in size. |
| Desktop, RTX 3060 + 31 GB RAM | split | 27B Q4 with about 60% of layers on the GPU | Runs, but slowly (single-digit tok/s) and RAM is tight (15 GB was in use by other things today). Not recommended as a worker. |
| Desktop, GTX 980 Ti as second GPU | +6 GB | — | Do not plan on it: CUDA 13 dropped that GPU generation; whether the server image can use it is unchecked. |

None of this helps until two things are done by Ernest: (1) reinstall `nvidia-container-toolkit` on the desktop (section 3a), (2) decide whether the Mac serves models to the LAN (bind address in its profile).

## Resume here

Done: inventory, GPU workaround, probe of Qwythos Q5 (four runs per task), assignments, candidates. In progress when this line was written: the same probe on `qwen35-9b-distilled-q8` and `qwen38-9b-distill-q5` for comparison (section 4a, below, if present). If a container `ym-local-llm` or `ym-probe-llama-gpu` is running, stop it. Open for others: reinstall the GPU container toolkit (Ernest); admission test of a 27B model on the Mac; build YM-T14 (`referee`, `index` first).
