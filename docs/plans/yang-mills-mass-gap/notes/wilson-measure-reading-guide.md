# One book and the first Wilson-measure study session

Recommendation for Ernest: Christof Gattringer and Christian B. Lang, **Quantum Chromodynamics on the Lattice: An Introductory Presentation**, Lecture Notes in Physics788, Springer2010. DOI: https://doi.org/10.1007/978-3-642-01850-3 . The publisher describes an introduction for newcomers with a progression from first principles to calculations. Its preface assumes quantum mechanics, some continuum QFT and basic particle physics, and promises detailed derivations of selected topics. That is a good fit for the requested Srednicki-like working style; it is a judgment of fit, not a claim the two authors write identically.

## What was checked

The filename inventory under the supplied Data/Public/books found Srednicki and the Frankel corpus, but no Gattringer/Lang title or other dedicated lattice-gauge textbook under the searched author/title names. This is not a content search of anonymously named book files. The repository's Frankel Ch20 notes already discuss connections/gauge potentials, which provide a bridge; they do not supply a complete Wilson-measure treatment. No book was purchased, downloaded wholesale or OCR'd.

Publisher book/chapter records and the publicly available publisher front matter/contents were opened. Full chapter text is not claimed read. Verified reading locators below are printed pages from that contents, not guessed PDF offsets.

- Book: https://link.springer.com/book/10.1007/978-3-642-01850-3
- Publisher front matter/contents: https://link.springer.com/content/pdf/bfm:978-3-642-01850-3/1
- Chapter2 record: https://link.springer.com/chapter/10.1007/978-3-642-01850-3_2
- Chapter3 record: https://link.springer.com/chapter/10.1007/978-3-642-01850-3_3

## Read in this order

| Reading | Printed pages | Stop when you can explain |
|---|---|---|
| §2.2.2–§2.2.3 |33–35| Why a gauge field is assigned to a link and how its transformation law fixes multiplication order |
| §2.3–§2.4 |36–41| Plaquette action, configuration integral and normalization |
| §3.1 |44–48| Normalized Haar integration; why compact-group gauge transformations preserve the reference measure |
| §3.3 |54–58| Closed-path/Wilson-loop observables |
| Later: §3.5.3–§3.5.4 |67–69| Running coupling and what taking a continuum limit requires |

If Euclidean correlators/path integrals need refreshing, first use Chapter1, especially §1.1.3 p4, §1.3.5–§1.3.6 pp16–18 and §1.4 pp19–23. Remain in the pure-gauge sector for this session; no need to read all fermion/spectroscopy chapters now. The book is background for these objects, not a solution or replacement for constructive continuum analysis.

## Today’s dictionary and learning checks

- ℓ counts sites per direction; a is physical lattice spacing; L=aℓ. The lattice coordinate vector e_mu is not the internal Lie-algebra basis e_a.
- G⊂U(m); m is matrix size. One group-valued variable per positive link:4ℓ^4 independent matrices. A reverse edge is a shifted inverse of an existing link.
- A plaquette is an **ordered** product around one square; its matrix is gauge covariant and its trace is invariant.
- Haar measure is the normalized invariant group measure. Product Haar is the reference law; the Wilson factor reweights it and generally couples links.
- Z normalizes the law. For beta>=0, exp(-12betaℓ^4)<=Z<=1; normalized Haar gives Z(0)=1.
- An observable is a function of the configuration; its expectation is the weighted group integral. A Wilson loop is one example, not the whole theory.
- For defining SU(m), Tr(TaTb)=delta_ab/2 implies g²=2m g0², beta=1/g0²=2m/g². Compare actions and trace conventions when comparing book/Lean/notes symbols.

Do these checks by hand before coding: identify all four plaquette corners and signed edges; derive cancellation of the three intermediate gauge transformations; derive Haar and Wilson-law invariance; show <1>=1 and Z(0)=1; translate SU2 beta4/g² and SU3 beta6/g² into the notes' g0. The expanded mathematical derivations are in the existing lattice topic, not repeated as code prose in the manuscript.

## Normalization finding and independent review

The earlier paragraph after ym:running-inverse claimed C'/g0'^2=C_g/g0². It was false: the invariant is the **product** g0'^2 C'=g0² C_g. With g²=c g0² and C'=C_g/c, beta=c/g² gives d beta/dlog(mu)=c*2b0'=11C_g/(24pi²)+O(g²). The final slope already written in the notes was correct; its intermediate justification is now corrected explicitly. Independent reviewer /root/ym_progress_audit accepted the first-measure lesson and this correction, source SHA256 a4b50005d75b4b529e994dcfe4bb4c771261ecd4d0244d01f17a7245bf5b6380. This review does not certify the later defect/mass-gap propositions.
