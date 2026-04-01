import { useState } from "react";
import { type Chapter } from "../data/chapters";

interface NotesPanelProps {
  chapter: Chapter;
}

const CONVENTIONS = [
  { label: "Metric", value: "g_μν = diag(+1, −1, −1, −1)  [mostly minus]" },
  { label: "Epsilon", value: "ε¹² = ε₂₁ = +1  (spinor index raising/lowering)" },
  { label: "Sigma", value: "σ⁰ = 1,  σⁱ = Pauli matrices" },
  { label: "Natural units", value: "ℏ = c = 1" },
];

export function NotesPanel({ chapter }: NotesPanelProps) {
  const [open, setOpen] = useState(true);

  return (
    <div className="notes-panel">
      <button
        className="notes-toggle"
        onClick={() => setOpen((v) => !v)}
        aria-expanded={open}
      >
        <span className="notes-toggle-icon">{open ? "▼" : "▶"}</span>
        <span className="notes-toggle-label">NOTES & CONVENTIONS</span>
      </button>

      {open && (
        <div className="notes-body">
          <div className="notes-description">
            <h3 className="notes-section-title">
              {typeof chapter.number === "number" ? `Ch.${chapter.number}` : chapter.number} — {chapter.title}
            </h3>
            <p className="notes-text">{chapter.description}</p>
          </div>

          <div className="notes-conventions">
            <h3 className="notes-section-title">Conventions</h3>
            <table className="conventions-table">
              <tbody>
                {CONVENTIONS.map((conv) => (
                  <tr key={conv.label} className="conventions-row">
                    <td className="conventions-label">{conv.label}</td>
                    <td className="conventions-value">
                      <code>{conv.value}</code>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
}
