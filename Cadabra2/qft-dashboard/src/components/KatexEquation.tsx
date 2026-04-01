import katex from "katex";
import "katex/dist/katex.min.css";

interface KatexEquationProps {
  latex: string;
  displayMode?: boolean;
  className?: string;
}

export function KatexEquation({ latex, displayMode = true, className = "" }: KatexEquationProps) {
  let html = "";
  let error = false;

  try {
    html = katex.renderToString(latex, {
      displayMode,
      throwOnError: false,
      strict: false,
    });
  } catch {
    error = true;
    html = `<span style="color:#ff6b6b">LaTeX error</span>`;
  }

  return (
    <span
      className={className}
      style={error ? { color: "#ff6b6b" } : undefined}
      dangerouslySetInnerHTML={{ __html: html }}
    />
  );
}
