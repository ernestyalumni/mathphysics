export function Header() {
  return (
    <header className="header">
      <div className="header-left">
        <span className="header-title">QFT Lab</span>
        <span className="header-subtitle">Srednicki · Spinors · MHV · AdS/CFT</span>
      </div>
      <div className="header-right">
        <span className="convention-badge">
          <span className="convention-label">metric</span>
          <span className="convention-value">(−,+,+,+)</span>
        </span>
        <a
          href="#"
          className="header-link"
          title="Srednicki QFT Textbook (placeholder)"
        >
          Srednicki PDF
        </a>
        <a
          href="https://github.com/InServiceOfX/Monoclaw"
          target="_blank"
          rel="noopener noreferrer"
          className="header-link"
        >
          GitHub
        </a>
      </div>
    </header>
  );
}
