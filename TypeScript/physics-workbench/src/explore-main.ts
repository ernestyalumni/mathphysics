import { mountExplore, type SourceRef } from "./explore";
import "./explore-standalone.css";
import "./explore-scene-layout.css";
const root = document.querySelector<HTMLDivElement>("#game")!;
root.innerHTML = `<header class="game-bar"><a class="game-brand" href="/explore.html">✦ Explore Physics</a><span>WORLDS OF MATHEMATICAL PHYSICS</span><button id="guide-toggle" type="button" aria-expanded="false">Field guide</button><label><input id="game-sound" type="checkbox"> Sound</label><button id="save-journey" type="button">Export journey</button></header><dialog id="field-guide"><button id="guide-close" type="button">Close guide</button><h2>Field guide</h2><p>Attach a source to a discovery without leaving your journey.</p><label>Collection <select id="source-root" aria-label="Reference collection"></select></label><label>Reference path <input id="source-path" aria-label="Reference path" placeholder="Path within collection, e.g. LaTeX_and_pdfs/Srednicki_QFT.tex"></label><button id="open-reference" type="button">Read reference</button><p id="guide-status" role="status"></p><pre id="guide-text"></pre></dialog><main id="journey"></main><footer class="game-footer">Explore freely · inspect the mathematics · keep your discoveries</footer>`;
const $ = <T extends HTMLElement>(id: string) =>
  document.getElementById(id) as T;
let sound = false;
try {
  sound = localStorage.getItem("mathphysics.explore.sound") === "on";
} catch {}
const checkbox = $<HTMLInputElement>("game-sound");
checkbox.checked = sound;
checkbox.onchange = () => {
  sound = checkbox.checked;
  try {
    localStorage.setItem("mathphysics.explore.sound", sound ? "on" : "off");
  } catch {}
};
let audio: AudioContext | undefined;
function cue(kind: "select" | "record" = "select") {
  if (!sound) return;
  try {
    audio ??= new AudioContext();
    void audio.resume();
    const notes = kind === "record" ? [523.25, 659.25, 783.99] : [554.37];
    notes.forEach((frequency, index) => {
      const o = audio!.createOscillator(),
        g = audio!.createGain(),
        at = audio!.currentTime + index * 0.065;
      o.type = "triangle";
      o.frequency.value = frequency;
      g.gain.setValueAtTime(0, at);
      g.gain.linearRampToValueAtTime(0.035, at + 0.008);
      g.gain.exponentialRampToValueAtTime(0.001, at + 0.1);
      o.connect(g).connect(audio!.destination);
      o.start(at);
      o.stop(at + 0.11);
    });
  } catch {}
}
const selector = $<HTMLSelectElement>("source-root"),
  path = $<HTMLInputElement>("source-path");
let roots = new Set<string>();
function source(): SourceRef | undefined {
  const value = path.value.trim();
  return roots.has(selector.value) && value
    ? {
        root: selector.value,
        path: value,
        name: value.split("/").at(-1) || value,
        kind: value.toLowerCase().endsWith(".pdf") ? "pdf" : "text",
      }
    : undefined;
}
async function openSource(s: SourceRef) {
  if (!$<HTMLDialogElement>("field-guide").open)
    $<HTMLDialogElement>("field-guide").showModal();
  $("guide-toggle").setAttribute("aria-expanded", "true");
  selector.value = s.root;
  path.value = s.path;
  const status = $("guide-status");
  try {
    if (!roots.has(s.root))
      throw Error("This reference collection is unavailable.");
    if (s.kind === "pdf" || s.path.toLowerCase().endsWith(".pdf")) {
      window.open(
        "/api/pdf?" + new URLSearchParams({ root: s.root, path: s.path }),
        "_blank",
        "noopener",
      );
      status.textContent = "PDF opened in a separate reader.";
      return;
    }
    status.textContent = "Opening reference…";
    const response = await fetch(
      "/api/text?" + new URLSearchParams({ root: s.root, path: s.path }),
    );
    const data = await response.json();
    if (!response.ok)
      throw Error(data.error?.message || "Reference unavailable");
    $("guide-text").textContent = data.text;
    status.textContent = s.name;
  } catch (e) {
    status.textContent = (e as Error).message;
  }
}
$("guide-toggle").onclick = () => {
  $<HTMLDialogElement>("field-guide").showModal();
  $("guide-toggle").setAttribute("aria-expanded", "true");
};
$("guide-close").onclick = () => $<HTMLDialogElement>("field-guide").close();
$<HTMLDialogElement>("field-guide").addEventListener("close", () =>
  $("guide-toggle").setAttribute("aria-expanded", "false"),
);
$("open-reference").onclick = () => {
  const s = source();
  if (s) void openSource(s);
  else
    $("guide-status").textContent =
      "Choose a collection and enter a reference path.";
};
void fetch("/api/roots")
  .then(async (response) => {
    if (!response.ok) throw Error("Reference collections unavailable");
    const data = await response.json();
    for (const r of data.roots) {
      roots.add(r.id);
      const o = document.createElement("option");
      o.value = r.id;
      o.textContent = r.label;
      selector.append(o);
    }
  })
  .catch((e) => ($("guide-status").textContent = (e as Error).message));
let view = mountExplore($("journey"), { sound: cue, source, openSource });
$("save-journey").onclick = () => view.exportState();
window.addEventListener("keydown", (e) => {
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "s") {
    e.preventDefault();
    view.exportState();
  }
});
window.addEventListener("pagehide", () => view.dispose());
window.addEventListener("pageshow", (e) => {
  if (e.persisted)
    view = mountExplore($("journey"), { sound: cue, source, openSource });
});
