import { useState } from "react";
import { type Chapter } from "../data/chapters";

interface CadabraLauncherProps {
  chapter: Chapter;
}

export function CadabraLauncher({ chapter }: CadabraLauncherProps) {
  const [copied, setCopied] = useState(false);

  const handleCopy = async () => {
    try {
      await navigator.clipboard.writeText(chapter.dockerCmd);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch {
      // fallback: select all text in the textarea
      const el = document.getElementById("docker-cmd-text");
      if (el instanceof HTMLTextAreaElement) {
        el.select();
        document.execCommand("copy");
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      }
    }
  };

  const scriptPath = `Python/Cadabra2/Srednicki/${chapter.scriptFile}`;

  return (
    <aside className="cadabra-panel">
      <div className="cadabra-header">
        <span className="cadabra-title">CADABRA2 LAUNCHER</span>
        <span className="cadabra-subtitle">Docker · CAS</span>
      </div>

      <div className="cadabra-script-path">
        <span className="cadabra-script-label">Script</span>
        <code className="cadabra-script-file">{scriptPath}</code>
      </div>

      <div className="cadabra-cmd-block">
        <label className="cadabra-cmd-label" htmlFor="docker-cmd-text">
          Docker command
        </label>
        <textarea
          id="docker-cmd-text"
          className="cadabra-cmd-textarea"
          readOnly
          value={chapter.dockerCmd}
          rows={6}
          spellCheck={false}
        />
      </div>

      <div className="cadabra-actions">
        <button
          className={`cadabra-btn cadabra-btn--copy ${copied ? "cadabra-btn--copied" : ""}`}
          onClick={handleCopy}
        >
          {copied ? "✓ Copied!" : "Copy command"}
        </button>
        <button
          className="cadabra-btn cadabra-btn--run"
          onClick={handleCopy}
          title="Copies the command to clipboard (backend execution out of scope for MVP)"
        >
          Run script
        </button>
      </div>

      <p className="cadabra-note">
        Run script copies the command to clipboard. Paste into a terminal with Docker running.
      </p>
    </aside>
  );
}
