import { useState } from "react";
import { Header } from "./components/Header";
import { ChapterSidebar } from "./components/ChapterSidebar";
import { EquationPanel } from "./components/EquationPanel";
import { CadabraLauncher } from "./components/CadabraLauncher";
import { NotesPanel } from "./components/NotesPanel";
import { LiteratureTab } from "./components/LiteratureTab";
import { chapters } from "./data/chapters";

type Tab = "srednicki" | "literature" | "tools";

const DEFAULT_TAB: Tab = "srednicki";
const DEFAULT_CHAPTER = "ch34";

export function App() {
  const [activeTab, setActiveTab] = useState<Tab>(DEFAULT_TAB);
  const [selectedChapterId, setSelectedChapterId] = useState(DEFAULT_CHAPTER);

  const selectedChapter = chapters.find((c) => c.id === selectedChapterId) ?? chapters[0]!;

  return (
    <div className="app-root">
      <Header />

      {/* Tab Navigation */}
      <div className="tab-bar">
        <button
          className={`tab-button ${activeTab === "srednicki" ? "tab-button--active" : ""}`}
          onClick={() => setActiveTab("srednicki")}
        >
          Srednicki Chapters
        </button>
        <button
          className={`tab-button ${activeTab === "literature" ? "tab-button--active" : ""}`}
          onClick={() => setActiveTab("literature")}
        >
          Literature & Reading List
        </button>
        <button
          className={`tab-button ${activeTab === "tools" ? "tab-button--active" : ""}`}
          onClick={() => setActiveTab("tools")}
        >
          Computations & Tools
        </button>
      </div>

      <div className="app-body">
        {activeTab === "srednicki" && (
          <>
            <ChapterSidebar
              selectedId={selectedChapterId}
              onSelect={setSelectedChapterId}
            />
            <main className="main-content">
              <EquationPanel chapter={selectedChapter} />
              <NotesPanel chapter={selectedChapter} />
            </main>
            <CadabraLauncher chapter={selectedChapter} />
          </>
        )}

        {activeTab === "literature" && (
          <div className="tab-full-content">
            <LiteratureTab />
          </div>
        )}

        {activeTab === "tools" && (
          <div className="tab-full-content">
            <div className="tools-placeholder">
              <div className="tools-placeholder-inner">
                <div className="tools-placeholder-icon">⚗️</div>
                <h2 className="tools-placeholder-title">Computations & Tools</h2>
                <p className="tools-placeholder-desc">
                  Cadabra2 utilities, Docker launchers, and symbolic computation tools
                  will be organized here.
                </p>
              </div>
              <CadabraLauncher chapter={selectedChapter} />
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
