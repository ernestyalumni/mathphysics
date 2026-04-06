import katex from "katex";
import "katex/dist/katex.min.css";

interface KeyEquation {
  label: string;
  latex: string;
  description?: string;
}

interface Paper {
  id: string;
  title: string;
  authors: string;
  year: number;
  arxivId: string;
  arxivUrl: string;
  priority: "Critical" | "High" | "Medium" | "Low";
  status: "Active" | "Backlog" | "Completed";
  description: string;
  keyEquations: KeyEquation[];
  latexFile: string;
}

function renderLatex(latex: string): string {
  try {
    return katex.renderToString(latex, {
      displayMode: true,
      throwOnError: false,
      strict: false,
    });
  } catch {
    return `<span style="color:#ff6b6b">LaTeX error</span>`;
  }
}

const PAPERS: Paper[] = [
  {
    id: "smga2025",
    title: "Single-minus gluon tree amplitudes are nonzero",
    authors: "Guevara, Lupsasca, Skinner, Strominger, Weil (OpenAI)",
    year: 2025,
    arxivId: "2602.12176",
    arxivUrl: "https://arxiv.org/abs/2602.12176",
    priority: "Critical",
    status: "Active",
    description:
      "OUR BENCHMARK. GPT-5.2 Pro conjectured a closed-form formula for single-minus gluon amplitudes in the half-collinear regime (Klein signature). A new OpenAI model proved it via Berends-Giele recursion. This is the discovery we want to match or exceed.",
    keyEquations: [
      {
        label: "Main Result (Eq. 16)",
        latex:
          "A_{1\\cdots n}\\Big|_{\\mathcal{R}_1} = \\frac{1}{2^{n-2}} \\prod_{m=2}^{n-1} \\left(\\mathrm{sg}_{m,m+1} + \\mathrm{sg}_{1,2\\cdots m}\\right)",
        description: "Closed-form single-minus stripped amplitude in the first region",
      },
      {
        label: "Stripped Amplitude Ansatz",
        latex:
          "\\mathcal{A}_n = i^{2-n} A_{1\\cdots n} \\prod_{a=2}^{n} \\delta(\\langle 1a\\rangle) \\cdot \\delta^2\\!\\left(\\sum_{i=1}^{n} \\langle ri\\rangle \\tilde\\lambda_i\\right)",
        description: "Full amplitude decomposition with half-collinear δ-functions",
      },
      {
        label: "Berends-Giele Recursion (SDYM)",
        latex:
          "\\mathcal{F}_{1\\cdots m} = \\frac{1}{p_{1\\cdots m}^2 + i\\varepsilon} \\sum_{j=1}^{m-1} [\\tilde\\lambda_{1\\cdots j}\\, \\tilde\\lambda_{j+1\\cdots m}]\\, \\mathcal{F}_{1\\cdots j}\\, \\mathcal{F}_{j+1\\cdots m}",
        description: "Off-shell current recursion for self-dual Yang-Mills (App. B)",
      },
      {
        label: "Weinberg Soft Theorem",
        latex:
          "\\lim_{\\omega_n \\to 0} A_{1\\cdots n} = \\tfrac{1}{2}\\left(\\mathrm{sg}_{n-1,n} + \\mathrm{sg}_{n,1}\\right) A_{1\\cdots n-1}",
        description: "Consistency check: soft limit of stripped amplitudes",
      },
    ],
    latexFile: "08-literature-review.tex",
  },
  {
    id: "dixon2013",
    title: "A brief introduction to modern amplitude methods",
    authors: "Lance J. Dixon",
    year: 2013,
    arxivId: "1310.5353",
    arxivUrl: "https://arxiv.org/abs/1310.5353",
    priority: "High",
    status: "Active",
    description:
      "Pedagogical review of spinor-helicity formalism, color ordering, MHV amplitudes via Parke-Taylor, and BCFW recursion. The bridge between Srednicki and modern methods.",
    keyEquations: [
      {
        label: "Parke-Taylor MHV Formula",
        latex:
          "A_n(1^-,2^-,3^+,\\ldots,n^+) = \\frac{\\langle 12\\rangle^4}{\\langle 12\\rangle\\langle 23\\rangle\\cdots\\langle n1\\rangle}",
        description: "Color-ordered MHV amplitude, all massless gluons (Dixon §3)",
      },
      {
        label: "BCFW Recursion Relation",
        latex:
          "A_n = \\sum_{\\text{diagrams}} A_L(\\hat{z}_P)\\,\\frac{1}{\\hat{P}^2}\\,A_R(\\hat{z}_P)",
        description: "On-shell recursion via complex momentum shift (Dixon §4)",
      },
      {
        label: "Gluon Polarization Vectors",
        latex:
          "\\varepsilon^+_\\mu(k;q) = \\frac{\\langle q|\\gamma_\\mu|k]}{\\sqrt{2}\\langle qk\\rangle}, \\qquad \\varepsilon^-_\\mu(k;q) = \\frac{[q|\\gamma_\\mu|k\\rangle}{\\sqrt{2}[qk]}",
        description: "Helicity polarization vectors in spinor-helicity notation (Dixon §2)",
      },
    ],
    latexFile: "08-literature-review.tex",
  },
  {
    id: "elvang2013",
    title: "Scattering Amplitudes in Gauge Theory and Gravity",
    authors: "Henriette Elvang & Yu-tin Huang",
    year: 2013,
    arxivId: "1308.1697",
    arxivUrl: "https://arxiv.org/abs/1308.1697",
    priority: "High",
    status: "Backlog",
    description:
      "THE modern amplitudes textbook. Covers spinor-helicity, BCFW, on-shell methods, superamplitudes, Grassmannian, and gravity amplitudes. Replaces needing many individual papers. Essential for going beyond Dixon.",
    keyEquations: [
      {
        label: "Spinor Momentum Decomposition",
        latex:
          "p_{\\alpha\\dot\\alpha} = \\lambda_\\alpha \\tilde\\lambda_{\\dot\\alpha}, \\quad p^2 = 0",
        description: "Massless momentum as spinor outer product",
      },
      {
        label: "Schouten Identity",
        latex:
          "\\langle ij\\rangle\\langle kl\\rangle + \\langle ik\\rangle\\langle lj\\rangle + \\langle il\\rangle\\langle jk\\rangle = 0",
        description: "Fundamental 2D identity for manipulating brackets",
      },
    ],
    latexFile: "08-literature-review.tex",
  },
  {
    id: "witten2003",
    title: "Perturbative Gauge Theory As A String Theory In Twistor Space",
    authors: "Edward Witten",
    year: 2003,
    arxivId: "hep-th/0312171",
    arxivUrl: "https://arxiv.org/abs/hep-th/0312171",
    priority: "High",
    status: "Backlog",
    description:
      "Started the modern amplitudes revolution. MHV amplitudes localize on curves in twistor space. SMGA explicitly references this paper — Witten noted single-minus tree amplitudes are supported at a point in twistor space.",
    keyEquations: [
      {
        label: "Twistor Variables",
        latex:
          "Z^I = (\\lambda_\\alpha,\\, \\mu^{\\dot\\alpha}), \\quad \\mu^{\\dot\\alpha} = x^{\\alpha\\dot\\alpha}\\lambda_\\alpha",
        description: "Penrose twistor transform: spacetime → twistor space",
      },
    ],
    latexFile: "08-literature-review.tex",
  },
  {
    id: "bcfw2005",
    title: "Direct Proof Of Tree-Level Recursion Relation in Yang-Mills Theory",
    authors: "Britto, Cachazo, Feng, Witten",
    year: 2005,
    arxivId: "hep-th/0501052",
    arxivUrl: "https://arxiv.org/abs/hep-th/0501052",
    priority: "Medium",
    status: "Backlog",
    description:
      "The original BCFW paper. Short and elegant. Proves that tree amplitudes can be reconstructed from their poles under a complex momentum shift. One of the most important results in modern amplitudes.",
    keyEquations: [
      {
        label: "BCFW Shift",
        latex:
          "\\hat\\lambda_i(z) = \\lambda_i + z\\lambda_j, \\quad \\hat{\\tilde\\lambda}_j(z) = \\tilde\\lambda_j - z\\tilde\\lambda_i",
        description: "Complex deformation preserving on-shell conditions",
      },
    ],
    latexFile: "08-literature-review.tex",
  },
  {
    id: "grassmannian2012",
    title: "Scattering Amplitudes and the Positive Grassmannian",
    authors: "Arkani-Hamed, Bourjaily, Cachazo, Goncharov, Postnikov, Trnka",
    year: 2012,
    arxivId: "1212.5605",
    arxivUrl: "https://arxiv.org/abs/1212.5605",
    priority: "Medium",
    status: "Backlog",
    description:
      "The amplituhedron and geometric approach to scattering amplitudes. Amplitudes as volumes of positive geometries. A potential source of new discoveries if we can connect SMGA-type results to this framework.",
    keyEquations: [
      {
        label: "Grassmannian Integral",
        latex:
          "\\mathcal{A}_{n,k} = \\int \\frac{d^{k \\times n} C}{\\text{GL}(k)} \\frac{\\delta^{k \\times 4}(C \\cdot \\tilde\\eta)\\, \\delta^{k \\times 2}(C \\cdot \\tilde\\lambda)}{(1\\cdots k)(2\\cdots k{+}1)\\cdots(n\\cdots k{-}1)}",
        description: "N^{k-2}MHV amplitude as Grassmannian contour integral",
      },
    ],
    latexFile: "08-literature-review.tex",
  },
];

function PriorityBadge({ priority }: { priority: Paper["priority"] }) {
  const cls =
    priority === "Critical"
      ? "badge badge--priority-critical"
      : priority === "High"
      ? "badge badge--priority-high"
      : priority === "Medium"
      ? "badge badge--priority-medium"
      : "badge badge--priority-low";
  return <span className={cls}>{priority}</span>;
}

function StatusBadge({ status }: { status: Paper["status"] }) {
  const cls =
    status === "Active"
      ? "badge badge--status-active"
      : status === "Completed"
      ? "badge badge--status-completed"
      : "badge badge--status-backlog";
  return <span className={cls}>{status}</span>;
}

function PaperCard({ paper }: { paper: Paper }) {
  return (
    <div className="lit-card">
      <div className="lit-card-header">
        <div className="lit-card-meta">
          <PriorityBadge priority={paper.priority} />
          <StatusBadge status={paper.status} />
          <span className="lit-card-year">{paper.year}</span>
        </div>
        <h2 className="lit-card-title">{paper.title}</h2>
        <div className="lit-card-authors">{paper.authors}</div>
        <div className="lit-card-links">
          <a
            href={paper.arxivUrl}
            target="_blank"
            rel="noopener noreferrer"
            className="lit-arxiv-link"
          >
            arXiv:{paper.arxivId}
          </a>
          <span className="lit-latex-ref">
            Notes: <code>{paper.latexFile}</code>
          </span>
        </div>
      </div>

      <div className="lit-card-description">{paper.description}</div>

      <div className="lit-equations-section">
        <div className="lit-equations-label">Key Equations</div>
        <div className="lit-equations-list">
          {paper.keyEquations.map((eq) => (
            <div key={eq.label} className="lit-equation-item">
              <div className="lit-equation-label">{eq.label}</div>
              <div
                className="lit-equation-body"
                dangerouslySetInnerHTML={{ __html: renderLatex(eq.latex) }}
              />
              {eq.description && (
                <div className="lit-equation-desc">{eq.description}</div>
              )}
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

export function LiteratureTab() {
  return (
    <div className="lit-tab">
      <div className="lit-tab-header">
        <h1 className="lit-tab-title">Literature & Reading List</h1>
        <p className="lit-tab-subtitle">
          Papers we're actively reading and incorporating into our amplitude notes.
          Detailed write-ups live in the <code>amplitudes/</code> directory.
        </p>
      </div>
      <div className="lit-cards-grid">
        {PAPERS.map((paper) => (
          <PaperCard key={paper.id} paper={paper} />
        ))}
      </div>
    </div>
  );
}
