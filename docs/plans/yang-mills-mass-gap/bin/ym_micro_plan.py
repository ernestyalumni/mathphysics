"""Task plan for `ym-micro prepare`: referee units, quote questions, LaTeX briefs and the two queues.

Written by the strategist on 2026-11-16 for an unattended ten-hour run of two local models:
a fast 9B model on the desktop GPU and a 27B model on the Mac. Everything a queue produces is
unreviewed raw material. Oracles in the briefs are hidden from the model and checked by the runner.
"""
import csv, json, pathlib, re

N1 = "openai/math/preprints/Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026/build"
N2 = "openai/math/preprints/Sharp-mass-bounds-for-the-two-dimensional-O4-model-September-23-2026/build"
N3 = "openai/math/preprints/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026/build"
LEAN = "yang-mills-group/OpenAI-YM-20261028"
KEYS = "scottnarmstrong/yang-mills-informalization/data/key-statements.tsv"
TOPIC = "documents/mono-notes/topics/yang-mills-lattice.tex"

QUESTIONS = {
    "YM-R01": (["main.tex", "sections/rotation.tex", "sections/annulus.tex", "sections/clusters.tex"], N1, [
        "State the main theorem on exponential decay with all its hypotheses. Quote it.",
        "What are the definitions of the quantities H(A) and W(A)? Quote the defining equations.",
        "Quote the statement of the main proposition of the section on the fourth-order rotation estimate.",
        "In the rotation estimate, which terms form a nonnegative square, and which terms are local terms that are controlled separately? Quote the sentences that say so.",
        "Describe the annulus: its geometry, what is pinned or fixed, and the family of test functions. Quote each.",
        "Where do the factors log N and (log N)^2 come from? Quote the lines where they first appear.",
        "What is the final inequality for the ratio of the twisted to the untwisted partition function? Quote it.",
        "What identity relates a crossing probability to the ratio of twisted to untwisted partition functions? Quote it.",
        "Which correlation inequalities are used (for example Ginibre, FKG, Campbell-Chayes), and which property of the model does each one need? Quote each place.",
        "Which steps use that the lattice is two-dimensional or planar? Quote each.",
        "Which steps use that the spins lie on a sphere, or reduce n components to three? Quote each.",
        "Which steps use ferromagnetism, that is nonnegative couplings? Quote each.",
    ]),
    "YM-R02": (["sections/introduction.tex", "sections/finite-volume.tex"], N2, [
        "Quote the definition of the defect quantity Delta_beta(n) and the theorem or theorems of the finite-volume section, with all constants.",
        "Which properties of the model does the finite-volume argument use: reflection positivity and in which directions, symmetry between the time and space periods, positivity of the transfer operator, compactness, anything else? Quote each place.",
        "Quote the inequality that bounds the weight of the excited states, a sum of powers of ratios of transfer eigenvalues, by the defect.",
        "How does the argument pass from one spatial size to a larger spatial volume? Quote the key inequalities.",
    ]),
    "YM-R03a": (["sections/preliminary.tex"], N2, [
        "For each subsection, quote one sentence that states what the subsection establishes.",
        "Where does control pass from spin correlations to rotation-invariant observables? Quote the statement and the key step.",
        "Quote the construction of the Hilbert space and of the transfer operator.",
        "Quote every statement about mixing or decay of correlations for all observables, with its constants.",
    ]),
    "YM-R03b": (["sections/comparison.tex", "sections/assembly.tex"], N2, [
        "What exactly is compared between different cutoffs, and why is the comparison uniform? Quote.",
        "What is the fixed reference box? Quote its definition and one place where it is used.",
        "Quote the passage that goes from relative densities to partition functions.",
        "Which earlier results does the assembly section cite as inputs? Quote each citation with its sentence.",
    ]),
    "YM-R04": (["sections/introduction.tex", "sections/trace.tex", "sections/continuum.tex", "sections/uniqueness.tex"], N3, [
        "Quote the statements of the two main theorems of the introduction.",
        "What do 'canonical' and 'without selecting subsequences' mean here? Quote the sentences that define or explain them.",
        "Quote the fixed-box estimate of the trace section and the statement that the physical mass gap is uniform in the cutoff.",
        "Quote the proposition that compares tori and passes test functions to the limit.",
        "How is the bounded bare-offset ambiguity removed? Quote the key statements, including anything about two incommensurable periods.",
    ]),
}

WHOLE_PAPER_QUESTIONS = [  # for the Mac only: the whole of N2 in one context
    ("YM-R02", "Across the whole paper, list every place where the defect quantity Delta_beta(n) or its recurrence is used after it is defined. Quote each."),
    ("YM-R02", "Which results proved in later sections are needed to verify the smallness hypothesis on the defect at the starting scale? Quote the statements."),
    ("YM-R03", "Trace the chain from the renormalization map to the final mass bound: quote, in order, the statement that closes each link of the chain."),
    ("YM-R03", "Which statements hold uniformly in the cutoff or in the volume? Quote each with the word that expresses the uniformity."),
    ("YM-G02", "Does the paper state that its transfer operator is nonnegative or positive, and where is that proved or cited? Quote."),
    ("YM-G06", "Is any step of the paper computer-assisted or does any step rely on explicit numerical constants? Quote each such place."),
]

BRIEFS = [
    ("b01-compact-self-adjoint", [], "transfer",
     "Title: Integral operators with continuous kernel on a compact space.\n"
     "Let X be a compact metric space with a finite Borel measure mu of full support and K a continuous real symmetric kernel on X x X. "
     "Define the integral operator T on L^2(X,mu). State and prove that T is bounded and self-adjoint, with the kernel computations displayed. "
     "State (proof omitted) that T is Hilbert-Schmidt, hence compact, and the spectral theorem for compact self-adjoint operators. "
     "State the theorem that gives Tr T as the integral of K over the diagonal when T is nonnegative, and say which hypothesis it needs. "
     "Then derive Tr T^n = sum of lambda_j^n and the kernel formula for T^n by composing kernels in displayed steps.",
     [r"Hilbert--?Schmidt", r"K\(x,\s*x\)", r"Mercer", r"\\lambda_j\^\{?n"]),
    ("b02-perron-jentzsch", [], "transfer",
     "Title: Positive kernels: a simple largest eigenvalue.\n"
     "State the Perron-Frobenius theorem for a matrix with strictly positive entries, and prove that the largest eigenvalue is simple and its "
     "eigenvector can be chosen strictly positive, with each inequality displayed. Then state the analogue for an integral operator with a continuous "
     "strictly positive symmetric kernel on a compact space (name the theorem). Deduce that such an operator has a gap log(lambda_0/lambda_1) > 0. "
     "End with one remark: why this finite-volume gap says nothing about a gap that is uniform in the volume.",
     [r"Jentzsch", r"simple", r"\\log", r"\\lambda_0", r"uniform"]),
    ("b03-peter-weyl-convolution", [], None,
     "Title: Peter-Weyl and convolution by a class function.\n"
     "Let G be a compact group with normalized Haar measure. State the Schur orthogonality relations for matrix coefficients of irreducible unitary "
     "representations, with the factor involving the dimension d_lambda, and the Peter-Weyl theorem. Define the character chi_lambda and a class "
     "function. Let f = sum over lambda of a_lambda chi_lambda be a continuous class function and (C_f phi)(g) = integral of f(g h^{-1}) phi(h) dh. "
     "Compute C_f on a matrix coefficient of the representation lambda in displayed steps, using Schur orthogonality, and find the scalar by which "
     "C_f acts on the block of lambda. State when C_f is a nonnegative operator.",
     [r"d_\{?\\lambda", r"\\frac\{a_\{?\\lambda\}?\}\{d_\{?\\lambda\}?\}|a_\{?\\lambda\}?\s*/\s*d_\{?\\lambda", r"\\delta", r"\\chi_\{?\\lambda"]),
    ("b04-su2-characters-wilson-weight", [], None,
     "Title: The Wilson weight of SU(2) in characters.\n"
     "For g in SU(2) with eigenvalues e^{i theta/2}, e^{-i theta/2}, 0 <= theta <= 2 pi, the character of the spin-j representation is "
     "chi_j(theta) = sin((2j+1) theta/2)/sin(theta/2), j = 0, 1/2, 1, .... State the Weyl integration formula for class functions on SU(2) in the "
     "variable theta with the normalized Haar measure, and verify the orthonormality of chi_0 and chi_{1/2} by explicit integrals. Then compute the "
     "coefficients a_j in exp((beta/2) tr g) = sum over j of a_j chi_j(g), for beta > 0, in displayed steps: write the integral for a_j, reduce it with "
     "product-to-sum identities, and identify the result with a standard special function using its integral representation. Give the final formula "
     "for a_j and say why a_j > 0.",
     [r"I_\{2j\s*\+\s*1\}", r"2j\s*\+\s*1", r"\\sin\^\{?2", r"\\frac\{2\s*\(2j\+1\)\}\{\\beta\}|2\(2j\+1\)\s*\\?,?\s*\\frac\{I_\{2j\+1\}\(\\beta\)\}\{\\beta\}|\(2j\+1\)\s*I_\{2j\+1\}\(\\beta\)\s*/\s*\\beta"]),
    ("b05-unbounded-spectral-theorem", [], None,
     "Title: Unbounded self-adjoint operators and the definition of a mass gap.\n"
     "Define a densely defined operator on a Hilbert space, its adjoint with its domain, symmetric versus self-adjoint. State the spectral theorem "
     "for a self-adjoint operator H in terms of a projection-valued measure E on the real line, and the formulas for f(H) and for the quadratic form "
     "<psi, f(H) psi> as an integral against the measure <psi, E(d lambda) psi>. Define the spectrum, point spectrum and continuous spectrum. "
     "For H >= 0 with a unique vector Omega such that H Omega = 0, define the mass gap Delta as a statement about spec(H), and derive in displayed "
     "steps the bound |<psi, e^{-tH} psi> - |<Omega,psi>|^2| <= e^{-t Delta} ||psi||^2 for unit Omega.",
     [r"projection-valued|projection valued", r"\[\s*\\Delta\s*,\s*\\infty", r"e\^\{-\s*t\s*\\Delta\}", r"\\Omega"]),
    ("b06-reflection-positivity-wilson", ["wilson-measure"], None,
     "Title: Reflection positivity of the Wilson measure.\n"
     "Use the notation of the supplied material. On a periodic lattice with an even number of sites in the time direction, define the reflection theta "
     "through the hyperplane midway between two time slices (link reflection), acting on links and on functions of the links. State what reflection "
     "positivity of the Wilson measure means: for every function F of the links in the positive half, the expectation of conj(theta F) F is >= 0. "
     "Give the standard proof in displayed steps: fix the temporal gauge on the links that cross the hyperplane, or integrate them explicitly; split the "
     "action into a part in the positive half, its reflection, and the plaquettes that cross; expand the crossing factor so that the integrand becomes a "
     "sum of terms conj(theta G) G with nonnegative coefficients. Mark with GAP: any step you cannot justify. End with the definition of the Hilbert "
     "space obtained from this positive semidefinite form.",
     [r"\\theta", r"\\ge(q)?\s*0", r"temporal|time", r"GAP:|nonnegative"]),
    ("b07-haar-su2-moments", [], None,
     "Title: Haar integrals on SU(2) by hand.\n"
     "Write U in SU(2) as U = q_0 1 + i (q_1 sigma_1 + q_2 sigma_2 + q_3 sigma_3) with q on the unit sphere S^3 in R^4, so that normalized Haar measure "
     "is the normalized round measure on S^3. Use hyperspherical coordinates q_0 = cos chi, with the measure (1/(2 pi^2)) sin^2 chi sin theta d chi d theta d phi, "
     "and check that it has total mass one. Compute, each in displayed steps with the final value as a fraction: (a) the integral of q_alpha q_beta; "
     "(b) the integral of q_0^4; (c) the integral of q_0^2 q_1^2; (d) the integral of |tr U|^2; (e) the integral of |U_11|^2; (f) the integral of U_11 U_22; "
     "(g) the integral of U_12 U_21. Use symmetry of the sphere wherever it shortens a computation, and say so.",
     [r"\\delta_\{\\alpha\\beta\}\s*/\s*4|\\frac\{\\delta_\{\\alpha\\beta\}\}\{4\}|\\t?frac\{?1\}?\{?4\}?\s*\\?,?\s*\\delta_\{\\alpha\\beta\}",
      r"\\t?frac\{?1\}?\{?8\}?|1/8", r"\\t?frac\{?1\}?\{?24\}?|1/24", r"-\s*\\t?frac\{?1\}?\{?2\}?|-\s*1/2"]),
    ("b08-matrix-integration-basics", [], None,
     "Title: Three meanings of 'integrating a matrix'.\n"
     "Distinguish, with one displayed example each: (1) the entrywise integral of a matrix-valued function of a real variable; (2) an integral over "
     "the vector space of real m x m matrices with Lebesgue measure on the entries; (3) an integral over a compact matrix group with Haar measure. "
     "For (1): show that the trace commutes with the integral, and that the determinant does not, by computing both the determinant of the integral "
     "and the integral of the determinant for A(t) = t 1_2 on 0 <= t <= 1. For (2): compute the Gaussian integral of exp(-(1/2) tr(A^T A)) over real "
     "2 x 2 matrices, reducing it to one-dimensional Gaussian integrals in displayed steps. For (3): state left and right invariance and normalization, "
     "and show that the integral of U over U(1) vanishes.",
     [r"\\t?frac\{?1\}?\{?4\}?|1/4", r"\\t?frac\{?1\}?\{?3\}?|1/3", r"\(2\\pi\)\^\{?2\}?|4\s*\\pi\^\{?2", r"invarian"]),
    ("b09-laplacian-compact-group", [], None,
     "Title: The Laplacian on a compact Lie group and the heat-kernel weight.\n"
     "Let G be a compact connected Lie group with a bi-invariant metric and Laplace-Beltrami operator. State that on the Peter-Weyl block of an "
     "irreducible representation lambda the Laplacian acts by a scalar given by the Casimir eigenvalue c_lambda >= 0, and write the heat kernel as a "
     "sum over lambda of d_lambda e^{-t c_lambda} chi_lambda. Specialize to SU(2) identified with the unit sphere S^3: give the eigenvalues of minus "
     "the Laplacian in terms of ell = 2j, their multiplicities, and check the multiplicity against the dimension of the Peter-Weyl block. "
     "Conclude that convolution by the heat kernel is a nonnegative operator with all character coefficients positive.",
     [r"\\ell\s*\(\s*\\ell\s*\+\s*2\s*\)", r"\(\s*\\ell\s*\+\s*1\s*\)\^\{?2", r"e\^\{-\s*t", r"d_\{?\\lambda"]),
    ("b10-transfer-matrix-correlators", [], None,
     "Title: Masses from Euclidean correlators.\n"
     "Let T = e^{-aH} be a positive self-adjoint transfer operator on a Hilbert space with discrete spectrum, eigenvalues e^{-a E_n}, E_0 < E_1 <= ..., "
     "on a lattice with N_t time slices and periodic boundary conditions. Starting from Z = Tr T^{N_t}, derive in displayed steps the spectral "
     "decomposition of the two-point function <O(t) O(0)> as a double sum over eigenstates, take N_t to infinity, and show that the connected correlator "
     "decays as e^{-t (E_1 - E_0)} times a matrix element squared. Define the effective mass from the ratio of the correlator at neighbouring times and "
     "state what it converges to. State exactly which assumption fails when the spectrum above the vacuum is continuous.",
     [r"e\^\{-", r"\\log", r"E_1\s*-\s*E_0", r"\\operatorname\{Tr\}|\\mathrm\{Tr\}|\\text\{Tr\}|\\Tr"]),
    ("b11-u1-two-dimensions", ["wilson-measure"], None,
     "Title: Exact solution of U(1) lattice gauge theory in two dimensions with free boundary conditions.\n"
     "Use the notation of the supplied material with G = U(1), plaquette weight exp(beta cos theta_p) up to normalization. On a rectangle with free "
     "boundary conditions, show by a change of variables (a gauge fixing along a maximal tree, displayed) that the plaquette angles are independent and "
     "uniformly distributed under Haar measure. Compute the expectation of e^{i theta_p} for one plaquette as a ratio of two standard special functions "
     "from their integral representations. Deduce the expectation of a rectangular Wilson loop enclosing A plaquettes and the string tension. "
     "State the behaviour of the string tension as beta tends to zero and to infinity.",
     [r"I_1\(\\beta\)", r"I_0\(\\beta\)", r"\^\{?A\}?", r"\\log"]),
    ("b12-exercises-wilson-measure", ["wilson-measure"], None,
     "Title: Exercises on the Wilson measure.\n"
     "Write six exercises with complete solutions on the supplied material, each in an exercise environment followed by a proof environment headed "
     "Solution. Topics, in this order: (1) the trace of a plaquette holonomy is gauge invariant; (2) reversing the orientation of a plaquette conjugates "
     "its trace for unitary groups; (3) the expectation of a single link variable vanishes; (4) the expectation of a gauge-invariant function is unchanged "
     "by a gauge transformation of the integration variables, using invariance of Haar measure; (5) for G = U(1) the plaquette energy reduces to "
     "1 - cos of the lattice curl; (6) the partition function with normalized Haar measure lies between e^{-2 beta |P|} and 1 for the energy "
     "1 - (1/m) Re tr. Every algebraic step displayed.",
     [r"\\begin\{exercise\}", r"gauge", r"Haar", r"\\cos"]),
    ("b13-z2-gauge-two-dimensions", [], None,
     "Title: Exact partition function of Z_2 lattice gauge theory on a two-dimensional torus.\n"
     "Links carry sigma = +1 or -1, the plaquette variable sigma_p is the product of the four links of p, the weight is the product over plaquettes of "
     "exp(beta sigma_p), and the measure on links is normalized counting measure. The torus has N plaquettes. Expand each factor as "
     "cosh(beta) (1 + sigma_p tanh(beta)), multiply out, and determine which sets of plaquettes survive the average over links, using that every link "
     "lies in exactly two plaquettes. Give the exact partition function in closed form. Then explain in one displayed line the constraint on the "
     "product of all plaquette variables that your answer reflects.",
     [r"\\cosh", r"\\tanh", r"1\s*\+\s*\\tanh\s*(\(\\beta\)|\\beta)?\s*\)?\^\{?N|\\sinh\s*(\(\\beta\)|\\beta)?\s*\)?\^\{?N", r"\\prod"]),
    ("b14-su2-two-dimensions-wilson-loop", ["wilson-measure"], None,
     "Title: The Wilson loop of SU(2) lattice gauge theory in two dimensions.\n"
     "Use the notation of the supplied material with G = SU(2) and plaquette weight exp((beta/2) tr U_p). With free boundary conditions the plaquette "
     "variables are independent Haar-distributed group elements. Using the character expansion of the weight and the orthogonality of characters, "
     "compute the expectation of (1/2) tr of the holonomy around a rectangle enclosing A plaquettes, by integrating out one plaquette at a time in "
     "displayed steps. Express the answer as a ratio of two special functions raised to the power A, and give its leading behaviour for small beta.",
     [r"I_2\(\\beta\)", r"I_1\(\\beta\)", r"\\frac\{\\beta\}\{4\}|\\beta\s*/\s*4", r"\^\{?A\}?"]),
]

PI_CONCEPTS = [
    ["the path integral for a free scalar field", "the LSZ reduction formula", "gauge invariance and the photon path integral", "nonabelian gauge theory: the Yang-Mills Lagrangian", "the Faddeev-Popov gauge fixing"],
    ["the beta function of nonabelian gauge theory", "asymptotic freedom", "Wilson loops and lattice gauge theory", "the renormalization group", "spontaneous breaking of gauge symmetry"],
]


def subsection_text(tex, title_fragment):
    heads = list(re.finditer(r"\\subsection\{([^}]*)\}", tex))
    for number, head in enumerate(heads):
        if title_fragment.lower() in head.group(1).lower():
            end = heads[number + 1].start() if number + 1 < len(heads) else len(tex)
            return tex[head.start():end]
    raise SystemExit("subsection not found: " + title_fragment)


def c23_unit(plan_dir):
    """The transfer-positivity sketch (claim C23) from the study-route note, as five steps."""
    text = (plan_dir / "notes/spectral-theory-study-route.md").read_text()
    section = text[text.index("## 3."):text.index("## 4.")]
    intro = section[:section.index("\n1. **")]
    steps = [" ".join(step.split()) for step in re.findall(r"^\d\. (\*\*.*?)(?=^\d\. \*\*|\Z)", section, re.M | re.S)]
    steps[-1] = steps[-1].split("The lesson of step 5")[0].strip()
    return {"name": "c23-transfer-positivity", "source": "docs/plans/yang-mills-mass-gap/notes/spectral-theory-study-route.md:section-3",
            "truth": "unknown", "context": " ".join(intro.split())[-1400:], "refs": [], "steps": steps}


def build(core, state, repos):
    plan_dir = core.PLAN
    state.mkdir(parents=True, exist_ok=True)
    units_dir, material_dir = state / "units", state / "material"
    units_dir.mkdir(exist_ok=True)
    material_dir.mkdir(exist_ok=True)
    topic_path = core.MP / TOPIC
    topic = core.strip_comments(topic_path.read_text())

    material = {"wilson-measure": material_dir / "wilson-measure.tex", "transfer": material_dir / "transfer-operator.tex"}
    core.atomic_write(material["wilson-measure"], subsection_text(topic, "Links, plaquettes, Wilson measure"))
    core.atomic_write(material["transfer"], subsection_text(topic, "Transfer operator and the purity"))

    def save(units):
        for unit in units:
            core.atomic_write(units_dir / (unit["name"] + ".json"), json.dumps(unit, indent=1))
        return [unit["name"] for unit in units]

    our_proofs = core.proof_units(topic_path, TOPIC, "ym") + [c23_unit(plan_dir)]
    our_sections = core.section_units(topic_path, TOPIC, "ym")
    n1_units = []
    for name in ("rotation", "annulus", "clusters"):
        n1_units += core.proof_units(repos / N1 / "sections" / (name + ".tex"), "N1/sections/%s.tex" % name, "n1-" + name, truth="source (N1 is formalized)")
    n2_units = core.proof_units(repos / N2 / "sections/finite-volume.tex", "N2/sections/finite-volume.tex", "n2-finite-volume", truth="source")
    proof_mutants = [mutant for unit in our_proofs for mutant in core.mutants(unit, limit=6)]
    section_mutants = [mutant for unit in our_sections for mutant in core.mutants(unit, limit=3)]
    names = {key: save(units) for key, units in (("our_proofs", our_proofs), ("our_sections", our_sections), ("n1", n1_units),
                                                 ("n2", n2_units), ("proof_mutants", proof_mutants), ("section_mutants", section_mutants))}

    def referee(group, seeds, repeatable, prefix="", max_tokens=5000):
        return [{"id": prefix + name, "type": "referee", "unit": str(units_dir / (name + ".json")), "seeds": seeds, "repeatable": repeatable,
                 "max_tokens": max_tokens} for name in names[group]]

    def quotes(packet, limit=None):
        files, base, questions = QUESTIONS[packet]
        return [{"id": "%s-q%02d" % (packet.lower(), number), "type": "quote", "packet": packet, "question": question,
                 "sources": [str(repos / base / name) for name in files]} for number, question in enumerate(questions[:limit], 1)]

    def drafts(max_tokens):
        return [{"id": name, "type": "tex-draft", "brief": brief, "oracle": oracle, "max_tokens": max_tokens,
                 "material": [str(material[key]) for key in keys] + ([str(material[extra])] if extra else [])}
                for name, keys, extra, brief, oracle in BRIEFS]

    with open(repos / KEYS, newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("link")]
    rank = {"main result": 0, "major result": 1, "key result": 2, "numbered result": 3}
    rows.sort(key=lambda row: rank.get(row["importance"], 4))

    def lean(kind, selection):
        return [{"id": row["number"], "type": kind, "k_id": row["number"], "lean_name": row["lean_name"], "link": row["link"],
                 "summary": row["summary"], "lean_root": str(repos / LEAN)} for row in selection]

    drafts_folders = [str(state / machine / "tex-draft" / verdict) for machine in ("desktop", "mac") for verdict in ("accepted", "rejected")]

    def draft_referee(tag, seeds):
        return [{"id": "drafts-" + tag, "type": "draft-referee", "folders": drafts_folders, "seeds": seeds, "repeatable": True, "max_tokens": 5000}]

    profile = state / "pi-profile"
    profile.mkdir(exist_ok=True)
    example = repos / "claw-work/Deployments/Scripts/llama-cpp-server/examples/pi"
    for name in ("models.json", "settings.json"):
        if (example / name).is_file():
            core.atomic_write(profile / name, (example / name).read_text())
    srednicki = repos.parent / "Data/Public/books/Physics/Srednicki-QuantumFieldTheory/Srednicki_QFT.tex"
    canonical = repos / "mathphysics"
    pi_tasks = []
    for number, concepts in enumerate(PI_CONCEPTS, 1):
        prompt = ("You have read-only tools: read, grep, find, ls. Task: build a reading map. For each concept below, use grep to find (1) where the textbook "
                  "source %s treats it and (2) where the notes under %s treat it. Concepts:\n%s\n\n"
                  "Answer with one markdown table and nothing else. Columns: concept | textbook locator | notes locator. A locator is a path with a line "
                  "number in backticks, like `documents/mono-notes/topics/example.tex:12`, taken from grep output. Write 'not found' if grep finds nothing. "
                  "Do not guess a line number." % (srednicki, canonical / "documents/mono-notes/topics", "\n".join("- " + c for c in concepts)))
        pi_tasks.append({"id": "r09-map-%d" % number, "type": "pi", "model": "llama-cpp-local/qwythos-9b-q5-k-m", "profile": str(profile),
                         "cwd": str(canonical), "prompt": prompt, "output": str(state / "desktop" / "pi" / ("r09-map-%d.md" % number)),
                         "tools": "read,grep,find,ls", "search_roots": [str(srednicki.parent), str(srednicki.parent.parent)], "must_match": [r"\|.*\|.*\|"], "minimum_paths": 3, "timeout": 1500})

    desktop = (referee("our_proofs", [0, 1, 2], True)
               + drafts(12000)
               + referee("proof_mutants", [0, 1, 2], False)
               + quotes("YM-R01") + quotes("YM-R02")
               + referee("our_sections", [0, 1, 2], True)
               + draft_referee("a", [0, 1, 2])
               + referee("n1", [0, 1, 2], True) + referee("n2", [0, 1, 2], True)
               + lean("lean-reading", rows[:150]) + lean("lean-compare", rows[:150])
               + pi_tasks
               + referee("section_mutants", [0, 1, 2], False)
               + quotes("YM-R03a") + quotes("YM-R03b") + quotes("YM-R04")
               + draft_referee("b", [0, 1, 2])
               + lean("lean-reading", rows[150:]) + lean("lean-compare", rows[150:])
               + draft_referee("c", [0, 1, 2]))

    whole_n2 = [str(repos / N2 / "sections" / (name + ".tex")) for name in (
        "introduction", "finite-volume", "preliminary", "linear", "renormalization", "calibration", "trajectories", "comparison", "assembly", "upper", "physical")]
    whole = [{"id": "n2-whole-q%02d" % number, "type": "quote", "packet": packet, "question": question, "sources": whole_n2, "max_tokens": 4000}
             for number, (packet, question) in enumerate(WHOLE_PAPER_QUESTIONS, 1)]
    probe = {"type": "probe", "n1": str(state / "sources/N1.txt"), "n2": str(state / "sources/N2.txt"), "lean_root": str(repos / LEAN)}
    mac_mutants = [task for task in referee("proof_mutants", [0, 1], False, max_tokens=4500) if re.search(r"--mut\d+-(flip|exponent)", task["id"])][::3]
    desktop_drafts = [{"id": "desktop-drafts", "type": "draft-referee", "folders": drafts_folders[:2], "seeds": [0], "repeatable": True, "max_tokens": 4500}]
    hard_sources = [task for task in referee("n2", [0], True, max_tokens=4500) + referee("n1", [0], True, max_tokens=4500)
                    if task["id"].startswith(("n2-", "n1-rotation-02"))]
    mac = (referee("our_proofs", [0, 1], True, max_tokens=4500)
           + [dict(probe, id="admission-2")]
           + drafts(10000)[:7]
           + mac_mutants
           + referee("our_sections", [0], True, max_tokens=4500)
           + drafts(10000)[7:]
           + hard_sources
           + quotes("YM-R02")
           + whole
           + quotes("YM-R04")
           + [dict(probe, id="admission-3")]
           + desktop_drafts)

    for name, tasks in (("desktop", desktop), ("mac", mac)):
        identifiers = [task["type"] + "/" + task["id"] for task in tasks]
        assert len(identifiers) == len(set(identifiers)), "duplicate task id in " + name
        core.atomic_write(state / ("queue-%s.jsonl" % name), "".join(json.dumps(task) + "\n" for task in tasks))
    summary = {key: len(value) for key, value in names.items()}
    summary.update(briefs=len(BRIEFS), desktop_tasks=len(desktop), mac_tasks=len(mac), lean_rows=len(rows))
    core.atomic_write(state / "plan-summary.json", json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    return 0
