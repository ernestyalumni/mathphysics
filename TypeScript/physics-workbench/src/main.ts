import { mountPdf, type PdfView } from "./pdf";
import "./style.css";
import "katex/dist/katex.min.css";
import { renderDocument, renderMath } from "./renderer";
import {
  restoreWorkspace,
  saveWorkspace,
  rememberPlace,
  type Cell,
  type Entry,
  type FileTab,
  type Place,
  type Root,
} from "./model";
import {
  getRoots,
  getFiles,
  getText,
  pdfUrl,
  runTrace,
  request,
  type Listing,
} from "./api";
const state = restoreWorkspace();
let roots: Root[] = [];
let listing: Listing | null = null;
let listGeneration = 0;
let listAbort: AbortController | undefined;
let searchTimer: number | undefined;
let audio: AudioContext | undefined;
let back: Place[] = [];
let forward: Place[] = [];
const blobs = new Map<string, string>();
const running = new Map<string, number>();
let runSequence = 0;
let requestedDocument = 0;
let activePdf: PdfView | undefined;
const $ = <T extends HTMLElement = HTMLElement>(selector: string) =>
  document.querySelector<T>(selector)!;
function el<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  className = "",
  text = "",
): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  e.className = className;
  e.textContent = text;
  return e;
}
function button(
  label: string,
  className: string,
  fn: () => void,
  title = label,
) {
  const b = el("button", className, label);
  b.type = "button";
  b.title = title;
  b.addEventListener("click", fn);
  return b;
}
function notify(message: string) {
  $("#toast").textContent = message;
  $("#toast").hidden = false;
  window.setTimeout(() => ($("#toast").hidden = true), 5500);
}
function persist() {
  if (!saveWorkspace(state))
    notify("Browser storage is full. Export your notebook to keep a copy.");
}
function sound(kind: "select" | "done" | "error" = "select") {
  if (!state.sound || state.style !== "pixel") return;
  try {
    audio ??= new AudioContext();
    void audio.resume();
    const t = audio.currentTime;
    const notes =
      kind === "done" ? [523, 784] : kind === "error" ? [196, 147] : [440, 554];
    notes.forEach((f, i) => {
      const o = audio!.createOscillator();
      const g = audio!.createGain();
      o.type = "triangle";
      o.frequency.value = f;
      g.gain.setValueAtTime(0, t + i * 0.07);
      g.gain.linearRampToValueAtTime(state.volume * 0.25, t + i * 0.07 + 0.008);
      g.gain.exponentialRampToValueAtTime(0.001, t + i * 0.07 + 0.075);
      o.connect(g).connect(audio!.destination);
      o.start(t + i * 0.07);
      o.stop(t + i * 0.07 + 0.08);
    });
  } catch {
    /* Visible feedback is independent of audio. */
  }
}
$("#app").innerHTML =
  `<div class="app-shell"><header class="topbar"><div class="brand"><span class="brand-mark" aria-hidden="true">∑</span><span><strong class="brand-name">mathphysics</strong><small>RESEARCH WORKBENCH</small></span></div><div class="top-actions"><span class="connection-badge" id="connection">Connecting to workspace…</span><div class="style-toggle" role="group" aria-label="Workspace style"><button id="style-notebook" type="button">Notebook</button><button id="style-pixel" type="button">Pixel</button></div><label class="sound-control"><input id="sound" type="checkbox"> Sound</label><input id="volume" type="range" min="0" max="1" step="0.05" aria-label="Sound volume" title="Sound volume"><button id="mobile-files" class="icon-button" aria-label="Toggle files">☰</button></div></header><div class="workspace"><aside class="sidebar"><div class="section-nav"><a id="launch-explore" class="nav-item" href="/explore.html" target="_blank" rel="noopener">✦ &nbsp; Launch Explore Physics ↗</a><button id="show-notebook" class="nav-item">▤ &nbsp; Research notebook</button><button id="open-local" class="nav-item">↥ &nbsp; Open a local file</button></div><div class="places" id="places"></div><section class="files-panel"><div class="sidebar-heading"><span>FILES</span><button id="favorite" class="icon-button" aria-label="Favorite this folder" title="Favorite this folder">☆</button></div><select id="root" aria-label="Workspace root"></select><label class="search-field"><span aria-hidden="true">⌕</span><input id="search" type="search" placeholder="Find files…" aria-label="Find files"><kbd>⌘ K</kbd></label><div class="files-toolbar"><button id="back" class="icon-button" aria-label="Back">←</button><button id="forward" class="icon-button" aria-label="Forward">→</button><button id="parent" class="icon-button" aria-label="Parent folder">↑</button><button id="refresh" class="icon-button" aria-label="Refresh files">↻</button><button id="copy-path" class="icon-button" aria-label="Copy folder path">⧉</button></div><form id="path-form"><input id="path" aria-label="Folder or file path" spellcheck="false" placeholder="Folder or file path"></form><nav id="breadcrumbs" class="breadcrumbs" aria-label="File location"></nav><div id="file-list" class="file-list" role="list" aria-label="Files"></div><p id="file-summary" class="progress-note" role="status"></p></section></aside><main class="document-area"><nav id="tabs" class="document-tabs" aria-label="Open documents"></nav><div id="document-toolbar" class="document-toolbar"></div><div id="document-body" class="document-body"></div></main><aside class="inspector"><div class="panel-title">WORKSPACE</div><div class="detail-row"><span>Engine</span><strong>Rust · local</strong></div><div class="detail-row"><span>Notebook</span><strong id="save-state">Saved locally</strong></div><div class="detail-row"><span>Style</span><strong id="style-label">Notebook</strong></div><div class="panel-title">CONNECTIONS</div><p class="progress-note">Available harnesses are listed below. A connection is shown only after an adapter is working.</p><div id="agents" class="agent-list"></div><button id="refresh-agents" class="small-button">Refresh availability</button><div class="panel-title">READING & SOURCES</div><p class="progress-note">Open PDFs and TeX side by side with your notebook. Use matching-file actions to move between source and compiled pages.</p><div class="welcome-card"><strong>Keep the assumptions.</strong><p>Calculated expressions, source notes and checked proofs carry different evidence.</p></div><div class="panel-title">SHORTCUTS</div><div class="detail-row"><span>Evaluate cell</span><kbd>Shift ↵</kbd></div><div class="detail-row"><span>Find a file</span><kbd>Ctrl/⌘ K</kbd></div><div class="detail-row"><span>Export document</span><kbd>Ctrl/⌘ S</kbd></div></aside></div><footer class="statusbar"><span id="status" role="status">Local workspace</span><span id="document-status">Notebook · draft autosave</span></footer><div id="toast" class="toast" role="status" hidden></div><input id="local-input" type="file" accept=".pdf,.tex,.md,.mmd,.txt,.json,.mathphysics.json" hidden></div>`;
function applyStyle() {
  document.documentElement.dataset.style = state.style;
  $("#style-notebook").setAttribute(
    "aria-pressed",
    String(state.style === "notebook"),
  );
  $("#style-pixel").setAttribute(
    "aria-pressed",
    String(state.style === "pixel"),
  );
  $("#style-label").textContent =
    state.style === "pixel" ? "Pixel" : "Notebook";
  $("#sound").toggleAttribute("checked", state.sound);
  ($("#sound") as HTMLInputElement).checked = state.sound;
  ($("#volume") as HTMLInputElement).value = String(state.volume);
  document.documentElement.style.setProperty(
    "--notebook-zoom",
    String(state.zoom / 100),
  );
}
for (const style of ["notebook", "pixel"] as const)
  $("#style-" + style).onclick = () => {
    state.style = style;
    applyStyle();
    persist();
    sound();
  };
$("#sound").addEventListener("change", () => {
  state.sound = ($("#sound") as HTMLInputElement).checked;
  persist();
  sound();
});
$("#volume").addEventListener("input", () => {
  state.volume = Number(($("#volume") as HTMLInputElement).value);
  persist();
});
$("#mobile-files").onclick = () => {
  const open = document.documentElement.dataset.sidebar !== "open";
  document.documentElement.dataset.sidebar = open ? "open" : "closed";
  $("#mobile-files").setAttribute("aria-expanded", String(open));
  if (open) $("#search").focus();
};
function activate(id: string) {
  requestedDocument++;
  state.activeTab = id;
  persist();
  renderTabs();
  renderDocumentArea();
  sound();
}
function renderTabs() {
  const tabs = $("#tabs");
  tabs.replaceChildren();
  const notebook = button("Research notebook", "document-tab", () =>
    activate("notebook"),
  );
  notebook.setAttribute(
    "aria-selected",
    String(state.activeTab === "notebook"),
  );
  notebook.setAttribute("role", "tab");
  tabs.append(notebook);
  for (const tab of state.tabs) {
    const wrap = el("div", "document-tab-wrap");
    const b = button(
      (tab.kind === "pdf" ? "▧ " : "▤ ") + tab.name,
      "document-tab",
      () => activate(tab.id),
    );
    b.setAttribute("role", "tab");
    b.setAttribute("aria-selected", String(state.activeTab === tab.id));
    b.title = tab.path;
    const close = button(
      "×",
      "tab-close",
      () => {
        const i = state.tabs.findIndex((t) => t.id === tab.id);
        state.tabs.splice(i, 1);
        if (blobs.has(tab.id)) {
          URL.revokeObjectURL(blobs.get(tab.id)!);
          blobs.delete(tab.id);
        }
        if (state.activeTab === tab.id) state.activeTab = "notebook";
        persist();
        renderTabs();
        renderDocumentArea();
      },
      "Close " + tab.name,
    );
    close.setAttribute("aria-label", "Close " + tab.name);
    wrap.append(b, close);
    tabs.append(wrap);
  }
}
function exportText(
  text: string,
  name: string,
  type = "text/plain;charset=utf-8",
) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = el("a");
  a.href = url;
  a.download = name;
  a.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 2000);
  notify("Exported " + name);
}
function exportCurrent() {
  const tab = state.tabs.find((t) => t.id === state.activeTab);
  if (!tab) {
    exportText(
      JSON.stringify(
        { format: "mathphysics.notebook.v1", cells: state.cells },
        null,
        2,
      ),
      "research.mathphysics.json",
      "application/json",
    );
    return;
  }
  if (tab.kind === "pdf") {
    const a = el("a");
    a.href = blobs.get(tab.id) || pdfUrl(tab.root, tab.path);
    a.download = tab.name;
    a.target = "_blank";
    a.rel = "noopener";
    a.click();
  } else exportText(tab.draft ?? tab.text ?? "", tab.name);
}
function renderDocumentArea() {
  activePdf?.dispose();
  activePdf = undefined;
  const toolbar = $("#document-toolbar");
  const body = $("#document-body");
  toolbar.replaceChildren();
  body.replaceChildren();
  body.style.setProperty("--notebook-zoom", String(state.zoom / 100));
  const tab = state.tabs.find((t) => t.id === state.activeTab);
  if (!tab) {
    state.activeTab = "notebook";
    toolbar.append(
      el("span", "toolbar-title", "NOTEBOOK"),
      button("+ Math", "small-button", () => addCell("math")),
      button("+ Text", "small-button", () => addCell("text")),
      button("+ Trace", "small-button", () => addCell("trace")),
      button("Export notebook", "primary-button", exportCurrent),
    );
    const zoom = el("select");
    zoom.setAttribute("aria-label", "Notebook zoom");
    for (const n of [100, 125, 150, 200]) {
      const option = el("option", "", n + "%");
      option.value = String(n);
      option.selected = n === state.zoom;
      zoom.append(option);
    }
    zoom.onchange = () => {
      state.zoom = Number(zoom.value);
      applyStyle();
      persist();
    };
    toolbar.append(zoom);
    const notebook = el("section", "notebook");
    notebook.append(
      el("div", "notebook-eyebrow", "MATHEMATICAL PHYSICS / WORKING NOTES"),
      el("h1", "document-heading", "Quantum field theory"),
      el(
        "p",
        "notebook-subtitle",
        "Working notes · Srednicki and mathematical physics",
      ),
    );
    for (const cell of state.cells) notebook.append(renderCell(cell));
    body.append(notebook);
    $("#document-status").textContent =
      state.cells.length + " cells · draft autosave";
    return;
  }
  toolbar.append(
    el(
      "span",
      "toolbar-title",
      tab.kind === "pdf" ? "PDF READER" : "SOURCE & MATH",
    ),
    el("span", "toolbar-path", tab.name),
  );
  if (tab.paired) {
    toolbar.append(
      button(
        tab.paired.kind === "pdf" ? "Open compiled PDF" : "Open TeX source",
        "small-button",
        () =>
          void openFile(
            {
              name: tab.paired!.path.split("/").pop()!,
              path: tab.paired!.path,
              kind: tab.paired!.kind,
              size: 0,
              modified_secs: 0,
            },
            tab.root,
          ),
      ),
    );
  }
  toolbar.append(
    button(
      tab.kind === "pdf" ? "Download PDF" : "Export draft",
      "small-button",
      exportCurrent,
    ),
  );
  if (tab.kind === "pdf") {
    const url =
      blobs.get(tab.id) ||
      (tab.root !== "local" ? pdfUrl(tab.root, tab.path) : "");
    if (!url) {
      body.append(
        el(
          "div",
          "empty-state",
          "Reopen this local PDF to restore access after a browser restart.",
        ),
      );
      return;
    }
    const page = el("input");
    page.type = "number";
    page.min = "1";
    page.value = String(tab.page || 1);
    page.setAttribute("aria-label", "PDF page");
    page.style.width = "64px";
    toolbar.append(el("label", "field-label", "Page"), page);
    const pdfHost = el("div", "pdf-view");
    const count = el("span", "progress-note", "Loading pages…");
    toolbar.append(
      button(
        "←",
        "icon-button",
        () => activePdf?.goTo((tab.page || 1) - 1),
        "Previous PDF page",
      ),
      button(
        "→",
        "icon-button",
        () => activePdf?.goTo((tab.page || 1) + 1),
        "Next PDF page",
      ),
    );
    page.onchange = () => activePdf?.goTo(Number(page.value));
    const zoom = el("select");
    zoom.setAttribute("aria-label", "PDF zoom");
    for (const [value, label] of [
      [1, "Fit width"],
      [0.75, "75% width"],
      [1.25, "125% width"],
      [1.5, "150% width"],
    ] as const) {
      const option = el("option", "", label);
      option.value = String(value);
      zoom.append(option);
    }
    zoom.onchange = () => activePdf?.setZoom(Number(zoom.value));
    toolbar.append(zoom);
    const native = el("a", "small-button", "Open in tab ↗");
    native.href = url;
    native.target = "_blank";
    native.rel = "noopener";
    toolbar.append(native);
    body.append(pdfHost);
    activePdf = mountPdf(pdfHost, url, tab.page || 1, (current, total) => {
      tab.page = current;
      page.value = String(current);
      page.max = String(total);
      count.textContent = current + " / " + total + " pages";
      $("#document-status").textContent =
        "PDF · page " + current + " of " + total + " · " + tab.name;
      persist();
    });
    toolbar.append(count);
    $("#document-status").textContent = "PDF · " + tab.path;
    return;
  }
  const split = el("div", "reader-split");
  const edit = el("textarea", "source-editor");
  edit.setAttribute("aria-label", "LaTeX or text source");
  edit.spellcheck = false;
  edit.value = tab.draft ?? tab.text ?? "";
  const preview = el("article", "latex-preview");
  const draw = () => {
    if (tab.kind === "tex" || tab.kind === "markdown")
      preview.innerHTML = renderDocument(edit.value, tab.kind);
    else preview.replaceChildren();
    if (tab.kind !== "tex" && tab.kind !== "markdown")
      preview.textContent = edit.value;
  };
  draw();
  let timer: number | undefined;
  edit.addEventListener("input", () => {
    tab.draft = edit.value;
    persist();
    window.clearTimeout(timer);
    timer = window.setTimeout(draw, 200);
  });
  split.append(edit, preview);
  body.append(split);
  $("#document-status").textContent =
    (tab.kind === "tex" ? "LaTeX source + math preview" : "Text preview") +
    " · local draft";
  if (!tab.text && tab.root !== "local") {
    preview.textContent = "Loading source…";
    void getText(tab.root, tab.path)
      .then((file) => {
        tab.text = file.text;
        tab.paired = file.paired;
        persist();
        if (state.activeTab === tab.id) renderDocumentArea();
      })
      .catch((e) => {
        preview.textContent = String(e.message || e);
      });
  }
}
function addCell(kind: Cell["kind"]) {
  if (state.cells.length >= 200) {
    notify(
      "This notebook has reached its 200-cell limit. Export it before starting another.",
    );
    return;
  }
  const cell: Cell = {
    id: crypto.randomUUID(),
    kind,
    source:
      kind === "trace"
        ? "mu, nu"
        : kind === "math"
          ? "E=mc^2"
          : "## New thought",
  };
  state.cells.push(cell);
  state.activeTab = "notebook";
  persist();
  renderTabs();
  renderDocumentArea();
  const editor = document.querySelector<HTMLTextAreaElement>(
    `[data-cell-id="${cell.id}"] textarea`,
  );
  editor?.focus();
  sound();
}
function renderCell(cell: Cell) {
  const box = el("section", "cell");
  box.dataset.cellId = cell.id;
  const label = el(
    "div",
    "cell-label",
    (cell.kind === "trace" ? "In" : "Note") +
      "[" +
      (state.cells.indexOf(cell) + 1) +
      "]:=",
  );
  const content = el("div", "cell-content");
  const controls = el("div", "command-row");
  controls.append(
    el(
      "span",
      "cell-kind",
      cell.kind === "trace"
        ? "Clifford trace"
        : cell.kind === "math"
          ? "Mathematical note"
          : "Text note",
    ),
  );
  const run = button(
    cell.kind === "trace" ? "Run cell ↵" : "Preview",
    "small-button",
    () => void evaluateCell(cell),
  );
  run.disabled = running.has(cell.id);
  if (cell.kind !== "trace") {
    const editButton = button("Edit source", "small-button", () => {
      const area = box.querySelector<HTMLTextAreaElement>("textarea")!;
      area.hidden = !area.hidden;
      editButton.textContent = area.hidden ? "Edit source" : "Hide source";
      if (!area.hidden) area.focus();
    });
    controls.append(editButton);
  }
  if (cell.kind === "trace") controls.append(run);
  controls.append(
    button(
      "×",
      "icon-button",
      () => {
        if (state.cells.length === 1) return;
        state.cells = state.cells.filter((c) => c.id !== cell.id);
        running.delete(cell.id);
        persist();
        renderDocumentArea();
      },
      "Remove cell",
    ),
  );
  const editor = el("textarea", "cell-editor");
  editor.setAttribute(
    "aria-label",
    (cell.kind === "trace" ? "Gamma labels" : "Cell source") +
      " " +
      (state.cells.indexOf(cell) + 1),
  );
  editor.spellcheck = false;
  editor.rows = cell.kind === "trace" ? 2 : 3;
  editor.hidden = cell.kind !== "trace";
  editor.value = cell.source;
  const output = el("div", "cell-output");
  output.dataset.output = cell.id;
  const paint = () => {
    output.replaceChildren();
    if (cell.kind === "math") output.innerHTML = renderMath(cell.source, true);
    else if (cell.kind === "text")
      output.innerHTML = renderDocument(cell.source, "markdown");
    else if (cell.output) paintTrace(output, cell.output);
    else
      output.append(
        el(
          "p",
          "progress-note",
          "Enter free labels separated by commas. Shift+Enter evaluates with the local Rust kernel.",
        ),
      );
  };
  paint();
  editor.oninput = () => {
    cell.source = editor.value;
    if (cell.kind === "trace") cell.output = undefined;
    persist();
    paint();
  };
  editor.onkeydown = (e) => {
    if (e.shiftKey && e.key === "Enter") {
      e.preventDefault();
      void evaluateCell(cell);
    }
  };
  content.append(controls, editor, output);
  box.append(label, content, el("div", "cell-bracket"));
  return box;
}
function paintTrace(target: HTMLElement, value: unknown) {
  target.replaceChildren();
  const data = value as {
    status?: string;
    result?: {
      polynomial?: {
        display: string;
        terms: { coefficient: string; factors: string[] }[];
      };
      assumptions?: string[];
    };
    error?: { message: string };
    evidence_kind?: string;
  };
  const status = data.status || "unknown";
  const pill = el(
    "span",
    "status-pill",
    status === "computed"
      ? "Computed · exact algebra"
      : status.replaceAll("_", " "),
  );
  pill.dataset.status = status;
  target.append(pill);
  if (status === "computed" && data.result?.polynomial) {
    const terms = data.result.polynomial.terms;
    const tex = terms.length
      ? terms
          .map((term, i) => {
            const n = Number(term.coefficient);
            const coeff =
              i > 0 && n >= 0 ? "+" + term.coefficient : term.coefficient;
            return (
              coeff +
              term.factors
                .map((f) => {
                  const m = f.match(
                    /^h\(\d+:([A-Za-z0-9_]+);\d+:([A-Za-z0-9_]+)\)$/,
                  );
                  return m
                    ? "\\,h^{" + labelTex(m[1]) + labelTex(m[2]) + "}"
                    : "\\,\\mathrm{metric}";
                })
                .join("")
            );
          })
          .join(" ")
      : "0";
    const math = el("div", "trace-math");
    math.innerHTML = renderMath(tex, true);
    target.append(math);
    const details = el("details");
    details.append(el("summary", "", "Assumptions & result details"));
    const assumptions = el(
      "p",
      "progress-note",
      data.result.assumptions?.join(" · ") ||
        "Formal symmetric metric; free labels do not imply sums.",
    );
    const code = el("pre", "result-json", JSON.stringify(data, null, 2));
    details.append(assumptions, code);
    target.append(details);
  } else
    target.append(
      el(
        "p",
        "progress-note",
        data.error?.message || "No computed result was returned.",
      ),
    );
}
function labelTex(label: string) {
  const greek: Record<string, string> = {
    mu: "\\mu",
    nu: "\\nu",
    rho: "\\rho",
    sigma: "\\sigma",
    alpha: "\\alpha",
    beta: "\\beta",
  };
  return greek[label] || "\\mathrm{" + label.replaceAll("_", "\\_") + "}";
}
async function evaluateCell(cell: Cell) {
  if (cell.kind !== "trace") {
    const target = document.querySelector<HTMLElement>(
      `[data-output="${cell.id}"]`,
    );
    if (target)
      target.innerHTML =
        cell.kind === "math"
          ? renderMath(cell.source, true)
          : renderDocument(cell.source, "markdown");
    sound();
    return;
  }
  const source = cell.source;
  const generation = ++runSequence;
  running.set(cell.id, generation);
  const target = document.querySelector<HTMLElement>(
    `[data-output="${cell.id}"]`,
  );
  if (target) target.textContent = "Calculating with Rust…";
  try {
    const labels = source.trim() ? source.split(",").map((s) => s.trim()) : [];
    const result = await runTrace(labels);
    if (running.get(cell.id) !== generation || cell.source !== source) {
      if (target && cell.source !== source)
        target.textContent =
          "Input changed. Evaluate again for a fresh result.";
      return;
    }
    cell.output = result;
    persist();
    if (target) paintTrace(target, result);
    sound(result.status === "computed" ? "done" : "error");
    $("#status").textContent =
      result.status === "computed"
        ? "Calculation complete · bounded exact algebra"
        : "Calculation returned " + String(result.status);
  } catch (e) {
    if (target)
      target.textContent =
        "Cannot reach the calculation service: " + String((e as Error).message);
    sound("error");
  } finally {
    if (running.get(cell.id) === generation) running.delete(cell.id);
  }
}
async function loadDirectory(place: Place, remember = true, history = true) {
  if (
    history &&
    (place.root !== state.location.root || place.path !== state.location.path)
  ) {
    back.push({ ...state.location });
    forward = [];
  }
  state.location = place;
  if (remember) rememberPlace(state, place);
  persist();
  renderPlaces();
  ($("#root") as HTMLSelectElement).value = place.root;
  ($("#path") as HTMLInputElement).value = place.path;
  await refreshFiles();
}
async function refreshFiles() {
  const generation = ++listGeneration;
  listAbort?.abort();
  listAbort = new AbortController();
  const root = state.location.root,
    path = state.location.path,
    q = ($("#search") as HTMLInputElement).value;
  $("#file-summary").textContent = q
    ? "Searching filenames…"
    : "Loading folder…";
  try {
    const result = await getFiles(root, path, q, listAbort.signal);
    if (generation !== listGeneration) return;
    listing = result;
    renderFileList();
    renderBreadcrumbs();
    $("#file-summary").textContent =
      result.entries.length +
      " " +
      (q ? "matches" : "items") +
      (result.truncated ? " · limited results; narrow your search" : "");
    $("#status").textContent =
      "Workspace connected · " +
      (roots.find((r) => r.id === root)?.label || root);
  } catch (e) {
    if ((e as Error).name === "AbortError") return;
    if (generation !== listGeneration) return;
    $("#file-list").replaceChildren(
      el("p", "empty-state", (e as Error).message),
    );
    $("#file-summary").textContent = "Folder unavailable";
  }
}
function renderFileList() {
  const list = $("#file-list");
  list.replaceChildren();
  for (const file of listing?.entries || []) {
    const b = button(
      "",
      "file-entry",
      () => {
        if (file.kind === "directory")
          void loadDirectory({
            root: state.location.root,
            path: file.path,
            name: file.name,
          });
        else void openFile(file, state.location.root);
      },
      file.path,
    );
    b.dataset.kind = file.kind;
    b.setAttribute("role", "listitem");
    b.append(
      el(
        "span",
        "file-icon",
        file.kind === "directory"
          ? "▰"
          : file.kind === "pdf"
            ? "PDF"
            : file.kind === "tex"
              ? "TᴇX"
              : "▤",
      ),
    );
    const text = el("span", "file-name", file.name);
    text.append(
      el(
        "small",
        "file-meta",
        ($("#search") as HTMLInputElement).value
          ? file.path
          : file.kind === "directory"
            ? "Folder"
            : file.kind.toUpperCase() + " · " + size(file.size),
      ),
    );
    b.append(text);
    list.append(b);
  }
  if (!listing?.entries.length)
    list.append(el("p", "empty-state", "No matching files in this folder."));
}
function size(n: number) {
  return n > 1024 * 1024
    ? (n / 1024 / 1024).toFixed(1) + " MB"
    : n > 1024
      ? Math.round(n / 1024) + " KB"
      : n + " B";
}
function renderBreadcrumbs() {
  const nav = $("#breadcrumbs");
  nav.replaceChildren();
  const root = roots.find((r) => r.id === state.location.root);
  nav.append(
    button(
      root?.label || state.location.root,
      "breadcrumb",
      () =>
        void loadDirectory({
          root: state.location.root,
          path: "",
          name: root?.label || "Workspace",
        }),
    ),
  );
  const pieces = state.location.path.split("/").filter(Boolean);
  pieces.forEach((part, i) =>
    nav.append(
      el("span", "", "/"),
      button(
        part,
        "breadcrumb",
        () =>
          void loadDirectory({
            root: state.location.root,
            path: pieces.slice(0, i + 1).join("/"),
            name: part,
          }),
      ),
    ),
  );
  $("#favorite").textContent = state.favorites.some(
    (p) => p.root === state.location.root && p.path === state.location.path,
  )
    ? "★"
    : "☆";
}
function renderPlaces() {
  const node = $("#places");
  node.replaceChildren();
  const group = (label: string, places: Place[]) => {
    if (!places.length) return;
    node.append(el("div", "sidebar-heading", label));
    for (const p of places) {
      const b = button(p.name, "place", () => void loadDirectory(p));
      b.title = p.root + ":/" + p.path;
      node.append(b);
    }
  };
  const defaults = roots.map((r) => ({
    root: r.id,
    path: r.id === "mathphysics" ? "LaTeX_and_pdfs" : "",
    name: r.id === "mathphysics" ? "Notes & derivations" : r.label,
  }));
  group("PLACES", [
    ...state.favorites,
    ...defaults.filter(
      (p) =>
        !state.favorites.some((f) => f.root === p.root && f.path === p.path),
    ),
  ]);
  group("RECENT", state.recent.slice(0, 4));
  const frequent = Object.entries(state.visits)
    .filter(([key, n]) => Number.isFinite(n) && n > 1 && key.includes(":"))
    .sort((a, b) => b[1] - a[1])
    .map(([key]) => {
      const split = key.indexOf(":");
      const root = key.slice(0, split),
        path = key.slice(split + 1);
      return {
        root,
        path,
        name:
          path.split("/").pop() ||
          roots.find((r) => r.id === root)?.label ||
          root,
      };
    })
    .filter((p) => roots.some((r) => r.id === p.root))
    .slice(0, 3);
  group("FREQUENT", frequent);
  if (state.recent.length)
    node.append(
      button(
        "Clear visit history",
        "small-button",
        () => {
          state.recent = [];
          state.visits = {};
          persist();
          renderPlaces();
        },
        "Clear recent and frequent folders; keep favorites",
      ),
    );
}
async function openFile(file: Entry, root: string) {
  const generation = ++requestedDocument;
  const existing = state.tabs.find(
    (t) => t.root === root && t.path === file.path,
  );
  if (existing) {
    activate(existing.id);
    return;
  }
  if (!["pdf", "tex", "markdown", "text"].includes(file.kind)) {
    notify("Preview supports PDF, LaTeX, Markdown and UTF-8 text files.");
    return;
  }
  if (state.tabs.length >= 30) {
    notify("Close a document tab before opening another.");
    return;
  }
  const tab: FileTab = {
    id: crypto.randomUUID(),
    root,
    path: file.path,
    name: file.name,
    kind: file.kind,
    page: 1,
  };
  try {
    if (file.kind !== "pdf") {
      const data = await getText(root, file.path);
      tab.text = data.text;
      tab.paired = data.paired;
    } else {
      const parent = file.path.split("/").slice(0, -1).join("/");
      const data = await getFiles(root, parent, "");
      const expected = file.name.replace(/\.pdf$/i, ".tex").toLowerCase();
      const pair = data.entries.find(
        (f) => f.name.toLowerCase() === expected && f.kind === "tex",
      );
      if (pair) tab.paired = { path: pair.path, kind: "tex" };
    }
    const duplicate = state.tabs.find(
      (t) => t.root === root && t.path === file.path,
    );
    if (duplicate) {
      if (generation === requestedDocument) activate(duplicate.id);
      return;
    }
    state.tabs.push(tab);
    if (generation === requestedDocument) {
      activate(tab.id);
      document.documentElement.dataset.sidebar = "closed";
    } else {
      persist();
      renderTabs();
    }
  } catch (e) {
    notify("Could not open " + file.name + ": " + (e as Error).message);
  }
}
$("#show-notebook").onclick = () => activate("notebook");
$("#root").onchange = () =>
  void loadDirectory({
    root: ($("#root") as HTMLSelectElement).value,
    path: "",
    name: "Workspace",
  });
$("#search").oninput = () => {
  clearTimeout(searchTimer);
  searchTimer = window.setTimeout(() => void refreshFiles(), 180);
};
$("#refresh").onclick = () => void refreshFiles();
$("#parent").onclick = () => {
  const path = state.location.path.split("/").slice(0, -1).join("/");
  void loadDirectory({
    root: state.location.root,
    path,
    name: path.split("/").pop() || "Workspace",
  });
};
$("#back").onclick = () => {
  const p = back.pop();
  if (p) {
    forward.push({ ...state.location });
    void loadDirectory(p, false, false);
  }
};
$("#forward").onclick = () => {
  const p = forward.pop();
  if (p) {
    back.push({ ...state.location });
    void loadDirectory(p, false, false);
  }
};
$("#favorite").onclick = () => {
  const i = state.favorites.findIndex(
    (p) => p.root === state.location.root && p.path === state.location.path,
  );
  if (i >= 0) state.favorites.splice(i, 1);
  else
    state.favorites.push({
      ...state.location,
      name:
        state.location.path.split("/").pop() ||
        roots.find((r) => r.id === state.location.root)?.label ||
        "Workspace",
    });
  persist();
  renderPlaces();
  renderBreadcrumbs();
  sound();
};
$("#copy-path").onclick = () => {
  const r = roots.find((r) => r.id === state.location.root);
  void navigator.clipboard
    .writeText(
      (r?.path || state.location.root) +
        (state.location.path ? "/" + state.location.path : ""),
    )
    .then(() => notify("Folder path copied."))
    .catch(() => notify("Select the folder path to copy it."));
};
$("#path-form").onsubmit = (e) => {
  e.preventDefault();
  let value = ($("#path") as HTMLInputElement).value.trim();
  let root = state.location.root;
  if (value.startsWith("/")) {
    const match = roots
      .filter((r) => value === r.path || value.startsWith(r.path + "/"))
      .sort((a, b) => b.path.length - a.path.length)[0];
    if (!match) {
      notify("Choose a path inside one of the configured workspace roots.");
      return;
    }
    root = match.id;
    value = value.slice(match.path.length).replace(/^\//, "");
  }
  const ext = value.split(".").pop()?.toLowerCase();
  const kind: Entry["kind"] =
    ext === "pdf"
      ? "pdf"
      : ext === "tex"
        ? "tex"
        : ext === "md" || ext === "mmd"
          ? "markdown"
          : [
                "txt",
                "json",
                "jsonl",
                "yaml",
                "yml",
                "csv",
                "rs",
                "toml",
                "prompt",
              ].includes(ext || "")
            ? "text"
            : "directory";
  if (kind === "directory")
    void loadDirectory({
      root,
      path: value,
      name: value.split("/").pop() || "Workspace",
    });
  else
    void openFile(
      {
        name: value.split("/").pop()!,
        path: value,
        kind,
        size: 0,
        modified_secs: 0,
      },
      root,
    );
};
$("#open-local").onclick = () => $("#local-input").click();
$("#local-input").onchange = async () => {
  const input = $("#local-input") as HTMLInputElement;
  const file = input.files?.[0];
  if (!file) return;
  try {
    if (file.name.endsWith(".mathphysics.json")) {
      const data = JSON.parse(await file.text());
      if (
        data.format !== "mathphysics.notebook.v1" ||
        !Array.isArray(data.cells) ||
        !data.cells.every(
          (c: Cell) =>
            c &&
            ["math", "text", "trace"].includes(c.kind) &&
            typeof c.source === "string",
        )
      )
        throw new Error("Unrecognized notebook file");
      const cells: Cell[] = data.cells.slice(0, 200).map((c: Cell) => ({
        id: crypto.randomUUID(),
        kind: c.kind,
        source: c.source,
      }));
      if (
        !confirm(
          "Replace this local notebook with the imported cells? Export the current notebook first if you need a separate copy.",
        )
      )
        return;
      state.cells = cells;
      activate("notebook");
      return;
    }
    const ext = file.name.split(".").pop()?.toLowerCase();
    const kind =
      ext === "pdf"
        ? "pdf"
        : ext === "tex"
          ? "tex"
          : ext === "md" || ext === "mmd"
            ? "markdown"
            : "text";
    if (kind !== "pdf" && file.size > 2 * 1024 * 1024)
      throw new Error("Text preview is limited to 2 MiB.");
    const tab: FileTab = {
      id: crypto.randomUUID(),
      root: "local",
      path: file.name,
      name: file.name,
      kind,
      page: 1,
    };
    if (kind === "pdf") blobs.set(tab.id, URL.createObjectURL(file));
    else tab.text = await file.text();
    state.tabs.push(tab);
    activate(tab.id);
  } catch (e) {
    notify((e as Error).message);
  } finally {
    input.value = "";
  }
};
async function refreshAgents() {
  try {
    const data = await request<{
      harnesses: {
        id: string;
        name: string;
        available: boolean;
        status: string;
        detail: string;
      }[];
    }>("/api/harnesses");
    const node = $("#agents");
    node.replaceChildren();
    for (const h of data.harnesses) {
      const row = el("div", "agent-row");
      row.append(
        el("span", "agent-dot", "●"),
        el("strong", "", h.name),
        el(
          "small",
          "",
          h.available ? "Available · not connected" : "Not found",
        ),
      );
      row.title = h.detail;
      node.append(row);
    }
  } catch {
    $("#agents").replaceChildren(
      el("p", "progress-note", "Connection service unavailable"),
    );
  }
}
$("#refresh-agents").onclick = () => void refreshAgents();
document.addEventListener("keydown", (e) => {
  if (
    e.key === "Escape" &&
    document.documentElement.dataset.sidebar === "open"
  ) {
    document.documentElement.dataset.sidebar = "closed";
    $("#mobile-files").setAttribute("aria-expanded", "false");
    $("#mobile-files").focus();
  }
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "k") {
    e.preventDefault();
    document.documentElement.dataset.sidebar = "open";
    $("#search").focus();
  }
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "s") {
    e.preventDefault();
    exportCurrent();
  }
});
window.addEventListener("pagehide", persist);
window.addEventListener("beforeunload", () => {
  for (const url of blobs.values()) URL.revokeObjectURL(url);
});
applyStyle();
renderTabs();
renderDocumentArea();
void (async () => {
  try {
    roots = (await getRoots()).roots;
    const select = $("#root") as HTMLSelectElement;
    for (const root of roots) {
      const option = el("option", "", root.label);
      option.value = root.id;
      select.append(option);
    }
    if (!roots.some((r) => r.id === state.location.root))
      state.location = {
        root: roots[0]?.id || "",
        path: "",
        name: roots[0]?.label || "Workspace",
      };
    $("#connection").textContent = "● Local workspace";
    await loadDirectory(state.location, false, false);
    void refreshAgents();
  } catch (e) {
    $("#connection").textContent = "Workspace offline";
    $("#file-summary").textContent =
      "Start the Rust workbench service to browse project files.";
    notify("Workspace service unavailable: " + (e as Error).message);
  }
})();
