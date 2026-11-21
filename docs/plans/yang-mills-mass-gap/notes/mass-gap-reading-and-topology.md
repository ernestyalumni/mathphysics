# Mass gap: reading list, the Higgs question, and topology

Label: unreviewed (claude-opus55-strat, 2026-11-20). Written for Ernest, a differential geometer, as background for the 4D Yang–Mills mass-gap program. Citations are from memory; check each one before quoting it.

## Reading list

**The problem**
- A. Jaffe and E. Witten, *Quantum Yang–Mills Theory*: the official Clay problem description.
- M. R. Douglas, *Report on the status of the Yang–Mills Millennium Prize problem* (2004).

**A mathematician's entry point**
- S. Chatterjee, *Yang–Mills for probabilists*, arXiv:1803.01950.
- S. Chatterjee, *A probabilistic mechanism for quark confinement*, arXiv:2006.16229.

**Differential geometry, already in our code**
- E. Witten, *On quantum gauge theories in two dimensions*, Commun. Math. Phys. 141 (1991).
- M. F. Atiyah and R. Bott, *The Yang–Mills equations over Riemann surfaces*, Phil. Trans. R. Soc. A 308 (1983).

The character-expansion formula Z = Σ_j d_j^χ (c_j/d_j)^V that `Rust/lattice-gauge` checks against Monte Carlo in 2D is Witten's result on the lattice. See `character::torus_plaquette_2d`, `torus_twist_2d` and `square_defect_2d`.

**Rigorous lattice gauge theory**
- K. Osterwalder and E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. 110 (1978): confinement and a mass gap at strong coupling.
- E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, Springer Lecture Notes in Physics 159 (1982).
- J. Glimm and A. Jaffe, *Quantum Physics: A Functional Integral Point of View*: cluster expansions.

**The 4D continuum frontier**
- T. Bałaban, the 1980s series in Commun. Math. Phys. on ultraviolet stability of lattice gauge theories.
- J. Magnen, V. Rivasseau and R. Sénéor, *Construction of YM₄ with an infrared cutoff*, Commun. Math. Phys. 155 (1993).

**Physics of confinement**
- J. Greensite, *An Introduction to the Confinement Problem*, Springer Lecture Notes in Physics 821 (2011).

Suggested first three: Chatterjee's survey, Göpfert–Mack (below) and 't Hooft 1978 (below). Together they give what the mass gap is, the one rigorous proof of a gap from topology, and the twist criterion our tools measure.

## Is the Higgs mechanism relevant?

**Not as the source of the gap.** Pure Yang–Mills has no scalar field to break a symmetry, and local gauge symmetry cannot break spontaneously:
- S. Elitzur, *Impossibility of spontaneously breaking local symmetries*, Phys. Rev. D 12 (1975).

The pure Yang–Mills gap is expected to be dynamical: dimensional transmutation turns the dimensionless coupling into a mass scale Λ.

Higgs models are still worth reading, for two reasons.

1. **Complementarity.** For a Higgs field in the fundamental representation, the "Higgs" and "confinement" regions are continuously connected: one phase, with no transition between them.
   - E. Fradkin and S. Shenker, *Phase diagrams of lattice gauge theories with Higgs fields*, Phys. Rev. D 19 (1979).
   - Osterwalder–Seiler (1978).

   Gapped behaviour with a Higgs field is easier to analyze, and in this sense it is the same phase as confinement.
2. **A positive control.** The counter-model filter rejects any argument that also works for U(1) in 4D, which has a massless Coulomb phase. A gauge–Higgs model is the opposite check: a known-gapped case where a good argument must succeed. A Higgs field is a small addition to `Rust/lattice-gauge`.

## Topology: where a gap mechanism might come from

**Centre vortices and twists**
- G. 't Hooft, *On the phase transition towards permanent quark confinement*, Nucl. Phys. B138 (1978).

This paper introduced twisted boundary conditions and vortex free energies. Z_tw/Z → 1 at large volume is the confinement criterion. YM-N02 measures it, and the 2D twist code (`Lattice::set_twist`, `examples/twist_ratio_2d.rs`) is its groundwork.

Plain thermodynamic integration is too noisy for ratios near 1: on the 4×4 torus at β = 6 it gave ln(Z_tw/Z) = −0.054(18) against exact −0.026. YM-N02 therefore needs a chain of intermediate ensembles.

E. T. Tomboulis claimed a proof of confinement along these lines (arXiv:0707.2179). The claim is disputed, and understanding why is instructive.

**Monopoles: the dual superconductor**
- A. M. Polyakov, *Quark confinement and topology of gauge groups*, Nucl. Phys. B120 (1977): monopoles give confinement and a mass gap in 3D compact U(1).
- M. Göpfert and G. Mack, *Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant*, Commun. Math. Phys. 82 (1982).

Göpfert–Mack is the rigorous version of Polyakov's result, and the only rigorous theorem deriving a gap from topology in a setting close to a continuum limit. It also shows why 4D is harder. In 4D, monopoles are worldlines, not points, and 4D U(1) has a massless Coulomb phase.

**Seiberg–Witten**
- N. Seiberg and E. Witten, Nucl. Phys. B426 (1994).

Confinement by monopole condensation in N = 2 supersymmetric Yang–Mills, broken softly to N = 1, is the one 4D case where a gap mechanism is understood exactly.

**Instantons**

Instantons matter for the θ-vacuum, the η′ mass and the topological susceptibility χ_t. A dilute instanton gas does not confine in 4D, so they are probably not the mechanism of the gap. They are still measurable on the lattice, by gradient flow and a field-theoretic charge Q, and χ_t is tabulated in Teper's review (hep-th/9812187). That makes them a good topology tool to build and validate first.

## Topology experiments queued

All in `Rust/lattice-gauge`, CPU path first, CUDA C++ where measured to pay:

1. Gradient (Wilson) flow and the clover topological charge Q. Check that Q is near an integer after flow, and compare χ_t with Teper's table.
2. Centre projection: maximal centre gauge, P-vortex density, and string tension from projected links (centre dominance).
3. Abelian projection: maximal abelian gauge, monopole currents, monopole density.
4. YM-N02 twist free energies in 4D, by a chain of ensembles.

Each must also be run on U(1) in 4D (no gap in the Coulomb phase). A diagnostic that looks the same there cannot be what produces the gap.
