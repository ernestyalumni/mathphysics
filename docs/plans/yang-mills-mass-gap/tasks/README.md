# Task packets

One file per task. A packet plus [../06-AGENT-PROTOCOL.md](../06-AGENT-PROTOCOL.md) is everything a worker needs. Status lives in [../BOARD.md](../BOARD.md) and on the tracker thread.

ID prefixes: `R` reading and digests, `A` audit, `G` gap research, `T` Rust tooling, `N` numerics, `O` operations, `W` watch.

Path shorthands used in every packet:

- `REPOS` = `/media/propdev/Expansion/openclaw/.openclaw/workspace/repos`
- `LEAN` = `REPOS/yang-mills-group/OpenAI-YM-20261028` (pin `3456e5e`)
- `YMG` = `REPOS/yang-mills-group/YM` (pin `55ba9fc`)
- `ARM` = `REPOS/scottnarmstrong/yang-mills-informalization` (pin `daf1db4`)
- `HAI` = `REPOS/hairer/humanise-yang-mills` (pin `16b070a9`)
- `OAI` = `REPOS/openai/math/preprints` (pin `4eaee87`)
- `N1` = `OAI/Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026`
- `N2` = `OAI/Sharp-mass-bounds-for-the-two-dimensional-O4-model-September-23-2026`
- `N3` = `OAI/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026`
- `PLAN` = this directory's parent

Digest template (all `R` packets): for each section read, give (1) the statements proved, with hypotheses and quantifiers in full; (2) the mechanism in your own words, with every displayed identity you rely on reproduced and its equation number; (3) what it uses from earlier sections; (4) what is specific to the model or dimension and why; (5) what you could not follow, marked `GAP:`. Quote theorem statements exactly. Length follows content; no padding.
