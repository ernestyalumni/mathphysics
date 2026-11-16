# Symbolic Wilson gauge words without lost conventions

As Ernest learning Wilson lattice gauge theory, I want reusable Rust symbolic operations that preserve group-word order, inverse orientation, gauge endpoints, trace meaning and normalization, so I can check the derivations in my notes and build toward exact lattice calculations without hidden convention changes.

Source: request-2026-11-16-wilson-reading.txt; first bounded proposal YM-T15. Draft independent acceptance source; no generated contract or implementation/semantic PDD validation claimed.

Acceptance: (1) Caller-supplied symbols and explicit assumptions, not chapter IDs, drive calculations. (2) Cyclic trace equivalence holds while an exact3-factor matrix witness rejects inappropriate swaps. (3) Path inversion reverses order and orientations; complex trace and normalized real trace remain distinct. (4) Positive dimension, scalar-domain assumptions, resource bounds and unsupported operations are explicit; no silent numerical/proof fallback. (5) Geometry and later Haar/partition/CUDA layers consume one owning algebra library; maintain independent oracles. (6) Existing users remain compatible and a small maintained example explains what was and was not checked. Finite calculations do not certify a continuum gap.
