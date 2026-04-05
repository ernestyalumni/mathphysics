import { chapters, roadmap, type Chapter, type ChapterStatus } from "../data/chapters";

const STATUS_ICON: Record<ChapterStatus, string> = {
  done: "✅",
  "in-progress": "🔄",
  "not-started": "⬜",
};

const STATUS_CLASS: Record<ChapterStatus, string> = {
  done: "status-done",
  "in-progress": "status-progress",
  "not-started": "status-todo",
};

interface ChapterSidebarProps {
  selectedId: string;
  onSelect: (id: string) => void;
}

function ChapterItem({
  chapter,
  isSelected,
  onSelect,
}: {
  chapter: Chapter;
  isSelected: boolean;
  onSelect: () => void;
}) {
  return (
    <button
      className={`chapter-item ${STATUS_CLASS[chapter.status]} ${isSelected ? "chapter-item--selected" : ""}`}
      onClick={onSelect}
    >
      <span className="chapter-icon">{STATUS_ICON[chapter.status]}</span>
      <span className="chapter-info">
        <span className="chapter-number">
          {typeof chapter.number === "number" ? `Ch.${chapter.number}` : chapter.number}
        </span>
        <span className="chapter-title-text">{chapter.title}</span>
      </span>
    </button>
  );
}

export function ChapterSidebar({ selectedId, onSelect }: ChapterSidebarProps) {
  const ordered = roadmap
    .map((id) => chapters.find((c) => c.id === id))
    .filter((c): c is Chapter => c !== undefined);

  return (
    <aside className="sidebar">
      <div className="sidebar-header">
        <span className="sidebar-title">CHAPTERS</span>
        <span className="sidebar-subtitle">Srednicki Roadmap</span>
      </div>
      <nav className="chapter-list">
        {ordered.map((chapter) => (
          <ChapterItem
            key={chapter.id}
            chapter={chapter}
            isSelected={chapter.id === selectedId}
            onSelect={() => onSelect(chapter.id)}
          />
        ))}
      </nav>
      <div className="sidebar-legend">
        <div className="legend-item">
          <span>✅</span> <span>Done</span>
        </div>
        <div className="legend-item">
          <span>🔄</span> <span>In Progress</span>
        </div>
        <div className="legend-item">
          <span>⬜</span> <span>Not Started</span>
        </div>
      </div>
    </aside>
  );
}
