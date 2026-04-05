import { type Chapter } from "../data/chapters";
import { KatexEquation } from "./KatexEquation";

interface EquationPanelProps {
  chapter: Chapter;
}

export function EquationPanel({ chapter }: EquationPanelProps) {
  return (
    <section className="equation-panel">
      <div className="equation-panel-header">
        <h2 className="equation-panel-title">
          <span className="eq-chapter-num">
            {typeof chapter.number === "number" ? `Ch.${chapter.number}` : chapter.number}
          </span>
          <span className="eq-chapter-title">{chapter.title}</span>
        </h2>
        <div className="equation-panel-divider" />
      </div>

      <div className="equations-grid">
        {chapter.equations.map((eq, i) => (
          <div key={i} className="equation-card">
            <div className="equation-label">{eq.label}</div>
            <div className="equation-body">
              <KatexEquation latex={eq.latex} displayMode={true} />
            </div>
            {eq.description && (
              <div className="equation-description">{eq.description}</div>
            )}
          </div>
        ))}
      </div>
    </section>
  );
}
