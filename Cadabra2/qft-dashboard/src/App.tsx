import { useState } from "react";
import { Header } from "./components/Header";
import { ChapterSidebar } from "./components/ChapterSidebar";
import { EquationPanel } from "./components/EquationPanel";
import { CadabraLauncher } from "./components/CadabraLauncher";
import { NotesPanel } from "./components/NotesPanel";
import { chapters } from "./data/chapters";

const DEFAULT_CHAPTER = "ch34";

export function App() {
  const [selectedId, setSelectedId] = useState(DEFAULT_CHAPTER);

  const chapter = chapters.find((c) => c.id === selectedId) ?? chapters[0]!;

  return (
    <div className="app-root">
      <Header />
      <div className="app-body">
        <ChapterSidebar selectedId={selectedId} onSelect={setSelectedId} />
        <main className="main-content">
          <EquationPanel chapter={chapter} />
          <NotesPanel chapter={chapter} />
        </main>
        <CadabraLauncher chapter={chapter} />
      </div>
    </div>
  );
}
