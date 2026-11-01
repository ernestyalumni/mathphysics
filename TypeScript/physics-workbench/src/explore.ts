import { mountScene } from "./explore-scene";
import { mountWorld } from "./explore-world";
import { mountPlots, type PlotSnapshot } from "./explore-plots";
import "./explore.css";

type Vec3 = [number, number, number];
export interface SourceRef {
  root: string;
  path: string;
  name: string;
  kind: string;
}
interface Input {
  axis: Vec3;
  t: number;
  slice: number;
}
interface Snapshot extends PlotSnapshot {
  convention: string;
  t: number;
  axis: Vec3;
  matrix: number[][][];
  checks: { unitarity_error: number; determinant_error: number };
}
interface Entry {
  id: string;
  title: string;
  created: string;
  input: Input;
  source?: SourceRef;
  note: string;
}
interface Saved {
  version: 1;
  convention: "pauli-plus-i-v1";
  location: "world" | "lab";
  input: Input;
  journal: Entry[];
}
export interface ExploreOptions {
  sound(kind?: "select" | "record"): void;
  source?(): SourceRef | undefined;
  openSource?(source: SourceRef): void;
}
const KEY = "mathphysics.explore.v1";
const defaults = (): Saved => ({
  version: 1,
  convention: "pauli-plus-i-v1",
  location: "world",
  input: { axis: [0, 0, 1], t: 0, slice: 0 },
  journal: [],
});
const finite = (v: unknown): v is number =>
  typeof v === "number" && Number.isFinite(v);
function validInput(x: unknown): x is Input {
  const a = x as Input;
  return (
    !!a &&
    Array.isArray(a.axis) &&
    a.axis.length === 3 &&
    a.axis.every(finite) &&
    Math.abs(Math.hypot(...a.axis) - 1) <= 1e-12 &&
    finite(a.t) &&
    a.t >= 0 &&
    a.t <= 4 * Math.PI &&
    finite(a.slice) &&
    Math.abs(a.slice) <= 1.2
  );
}
function validSource(s: unknown): s is SourceRef {
  return (
    !!s &&
    typeof s === "object" &&
    ["root", "path", "name", "kind"].every(
      (k) =>
        typeof (s as Record<string, unknown>)[k] === "string" &&
        (s as Record<string, string>)[k].length <= 2000,
    )
  );
}
function decode(raw: string): Saved {
  if (raw.length > 2_000_000) throw Error("Saved exploration is too large.");
  const s = JSON.parse(raw) as Saved;
  if (
    s.version !== 1 ||
    s.convention !== "pauli-plus-i-v1" ||
    !["world", "lab"].includes(s.location) ||
    !validInput(s.input) ||
    !Array.isArray(s.journal) ||
    s.journal.length > 100 ||
    !s.journal.every(
      (e) =>
        e &&
        typeof e.id === "string" &&
        typeof e.title === "string" &&
        e.title.length <= 200 &&
        typeof e.created === "string" &&
        validInput(e.input) &&
        typeof e.note === "string" &&
        e.note.length <= 2000 &&
        (!e.source || validSource(e.source)),
    )
  )
    throw Error("Unsupported or invalid exploration data.");
  return s;
}
const clone = <T>(v: T): T => JSON.parse(JSON.stringify(v)) as T;
function node<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  cls = "",
  text = "",
) {
  const e = document.createElement(tag);
  e.className = cls;
  e.textContent = text;
  return e;
}
function btn(text: string, fn: () => void, cls = "small-button") {
  const b = node("button", cls, text);
  b.type = "button";
  b.onclick = fn;
  return b;
}
function download(text: string, name: string) {
  const url = URL.createObjectURL(
    new Blob([text], { type: "application/json" }),
  );
  const a = node("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
function matrixText(s: Snapshot) {
  return s.matrix
    .map((row) =>
      row
        .map(
          ([re, im]) =>
            `${Math.abs(re) < 1e-12 ? 0 : re.toFixed(5)} ${im < 0 ? "−" : "+"} ${Math.abs(im) < 1e-12 ? "0" : Math.abs(im).toFixed(5)}i`,
        )
        .join("    "),
    )
    .join("\n");
}
function inputKey(x: Input) {
  return JSON.stringify(x);
}
export function mountExplore(
  host: HTMLElement,
  options: ExploreOptions,
): { dispose(): void; exportState(): void } {
  let saved = defaults(),
    blocked = false,
    raw = "";
  try {
    raw = localStorage.getItem(KEY) || "";
    if (raw) saved = decode(raw);
  } catch {
    blocked = true;
  }
  let disposed = false;
  let sequence = 0;
  let abort: AbortController | undefined;
  let timer: number | undefined;
  let frame = 0;
  let playing = false;
  let start = 0;
  let startT = 0;
  let lastRequest = 0;
  let snapshot: Snapshot | undefined;
  let snapshotKey = "";
  let inFlight = false;
  let compare: Entry | undefined;
  let scene: ReturnType<typeof mountScene> | undefined;
  let interactionTimer: number | undefined;
  let dockPlay: HTMLButtonElement | undefined;
  let axisCycle: HTMLButtonElement | undefined;
  let keepRequestedKey: string | undefined;
  let feedbackTimer: number | undefined;
  const shell = node("section", "explore-shell");
  shell.setAttribute("aria-label", "Mathematical physics exploration");
  const header = node("header", "explore-header");
  header.append(
    node("div", "explore-eyebrow", "FIELD STATION 01 · LIE GROUPS"),
    node("h1", "", "The symmetry observatory"),
    node(
      "p",
      "",
      "Explore a world. Change a transformation. Keep what you discover.",
    ),
  );
  const nav = node("div", "explore-actions");
  const worldButton = btn("Explore the grounds", () => show("world"));
  const labButton = btn("Open SU(2) lab", () => show("lab"), "primary-button");
  nav.append(worldButton, labButton);
  header.append(nav);
  shell.append(header);
  const notice = node("p", "progress-note");
  notice.setAttribute("role", "status");
  shell.append(notice);
  const recovery = node("div", "explore-recovery");
  if (blocked) {
    recovery.append(
      node(
        "p",
        "",
        "Saved exploration data could not be read. It is preserved; new changes will not overwrite it.",
      ),
      btn("Export saved data", () =>
        download(raw, "exploration-recovery.json"),
      ),
      btn("Reset exploration storage", () => {
        if (
          !window.confirm(
            "Reset only exploration data? Export your saved data first if needed.",
          )
        )
          return;
        keepRequestedKey = undefined;
        mark.removeAttribute("aria-busy");
        compare = undefined;
        comparison.textContent = "";
        blocked = false;
        saved = defaults();
        persist();
        recovery.remove();
        show("world");
      }),
    );
    shell.append(recovery);
  }
  const world = node("div", "explore-world-host");
  const lab = node("section", "explore-lab");
  lab.hidden = true;
  const worldView = mountWorld(
    world,
    () => show("lab"),
    () => options.sound(),
  );
  shell.append(world, lab);
  host.append(shell);
  const labHead = node("div", "explore-lab-heading");
  labHead.append(
    node("h2", "", "SU(2) laboratory"),
    node(
      "p",
      "",
      "Turn the dial through 4π. The group element changes sign after 2π and returns after 4π.",
    ),
  );
  lab.append(labHead);
  const layout = node("div", "explore-layout");
  const controls = node("section", "explore-controls");
  controls.setAttribute("aria-label", "SU(2) controls");
  layout.append(controls);
  const panels = node("div", "explore-panels");
  layout.append(panels);
  lab.append(layout);
  const axisLabel = node("label", "", "Rotation axis");
  const axis = node("select");
  axis.setAttribute("aria-label", "Rotation axis");
  for (const a of ["x", "y", "z", "custom"]) {
    const o = node(
      "option",
      "",
      a === "custom" ? "Custom unit axis" : a + " axis",
    );
    o.value = a;
    axis.append(o);
  }
  axisLabel.append(axis);
  controls.append(axisLabel);
  const custom = node("div", "explore-axis");
  const components = [0, 1, 2].map((i) => {
    const l = node("label", "", "n" + ["x", "y", "z"][i]);
    const e = node("input");
    e.type = "number";
    e.step = ".1";
    e.value = String(saved.input.axis[i]);
    e.setAttribute("aria-label", "Axis " + ["x", "y", "z"][i]);
    l.append(e);
    custom.append(l);
    return e;
  });
  const applyAxis = () => {
    const v = components.map((e) => Number(e.value)) as Vec3;
    if (!v.every(finite) || Math.abs(Math.hypot(...v) - 1) > 1e-12) {
      notice.textContent =
        "Axis must be finite and unit length. Use Normalize axis explicitly if desired.";
      return;
    }
    saved.input.axis = v;
    changed();
  };
  custom.append(
    btn("Apply axis", applyAxis),
    btn("Normalize axis", () => {
      const v = components.map((e) => Number(e.value));
      const norm = Math.hypot(...v);
      if (!finite(norm) || norm < 1e-12) {
        notice.textContent = "A zero or non-finite axis cannot be normalized.";
        return;
      }
      components.forEach((e, i) => (e.value = String(v[i] / norm)));
      applyAxis();
    }),
  );
  controls.append(custom);
  const timeLabel = node("label", "", "Rotation angle t · radians");
  const time = node("input");
  time.type = "range";
  time.min = "0";
  time.max = String(4 * Math.PI);
  time.step = ".001";
  time.setAttribute("aria-label", "Rotation angle t");
  timeLabel.append(time);
  const timeRead = node("output");
  timeRead.id = "explore-angle";
  controls.append(timeLabel, timeRead);
  const checkpoints = node("div", "explore-checkpoints");
  for (const [label, t] of [
    ["0", 0],
    ["π", Math.PI],
    ["2π", 2 * Math.PI],
    ["4π", 4 * Math.PI],
  ] as const)
    checkpoints.append(
      btn(label, () => {
        pause();
        saved.input.t = t;
        changed(true);
      }),
    );
  controls.append(checkpoints);
  const playback = node("div", "explore-actions");
  const play = btn("Play", () => {
    if (playing) {
      pause();
      return;
    }
    if (saved.input.t >= 4 * Math.PI) saved.input.t = 0;
    playing = true;
    start = performance.now();
    startT = saved.input.t;
    lastRequest = 0;
    play.textContent = "Pause";
    setDockPlayLabel("Pause");
    frame = requestAnimationFrame(tick);
  });
  playback.append(
    play,
    btn("Step +π/12", () => {
      pause();
      saved.input.t = Math.min(4 * Math.PI, saved.input.t + Math.PI / 12);
      changed(true);
    }),
    btn("Reset to identity", () => {
      pause();
      saved.input.t = 0;
      changed(true);
    }),
  );
  controls.append(
    playback,
    node(
      "p",
      "progress-note",
      "Playback: 1 radian per second. t is an angle, not physical time. The 4π endpoint stays visible.",
    ),
  );
  const sliceLabel = node("label", "", "Slice a₀ = c");
  const slice = node("input");
  slice.type = "range";
  slice.min = "-1.2";
  slice.max = "1.2";
  slice.step = ".01";
  slice.setAttribute("aria-label", "Slice coordinate");
  sliceLabel.append(slice);
  const sliceRead = node("output");
  controls.append(sliceLabel, sliceRead);
  const evidence = node("p", "explore-evidence", "Waiting for Rust");
  evidence.setAttribute("role", "status");
  controls.append(evidence);
  const matrix = node("pre", "explore-matrix");
  matrix.setAttribute("aria-label", "SU(2) matrix");
  const equation = node(
    "p",
    "progress-note",
    "U = a₀I + i a·σ = exp(−i t n·σ / 2). α = a₀ + i a₃; β = a₂ + i a₁.",
  );
  controls.append(equation, matrix);
  const plotHost = node("div", "explore-plots");
  panels.append(plotHost);
  const plots = mountPlots(plotHost);
  const pinRow = node("div", "explore-record");
  const title = node("input");
  title.placeholder = "Name this experiment";
  title.value = "A turn through SU(2)";
  title.maxLength = 200;
  title.setAttribute("aria-label", "Experiment title");
  const note = node("textarea");
  note.placeholder = "What did you notice?";
  note.maxLength = 2000;
  note.setAttribute("aria-label", "Experiment observation");
  const attach = node("input");
  attach.type = "checkbox";
  const attachLabel = node("label", "", "Attach field-guide reference");
  attachLabel.prepend(attach);
  const save = btn(
    "Record experiment",
    () => {
      if (playing) {
        pause();
        notice.textContent =
          "Playback paused. Wait for Rust, then record the point.";
        return;
      }
      if (!snapshot || snapshotKey !== inputKey(saved.input)) {
        notice.textContent =
          "Wait for the current Rust calculation before recording.";
        return;
      }
      const entry: Entry = {
        id: crypto.randomUUID(),
        title: title.value.trim() || "SU(2) experiment",
        created: new Date().toISOString(),
        input: clone(saved.input),
        note: note.value,
      };
      const source = options.source?.();
      if (attach.checked && source) entry.source = clone(source);
      if (saved.journal.length >= 100) {
        notice.textContent =
          "Journal is full (100 experiments). Export or remove an entry first.";
        return;
      }
      saved.journal.unshift(entry);
      const durable = persist();
      renderJournal();
      options.sound("record");
      celebrate(durable);
      notice.textContent = durable
        ? "Experiment saved with its exact inputs. Numerical check, not a proof."
        : "Experiment is in memory only. Storage is unavailable; export exploration to keep it.";
    },
    "primary-button",
  );
  pinRow.append(title, note, attachLabel, save);
  lab.append(pinRow);
  const journal = node("section", "explore-journal");
  journal.setAttribute("aria-label", "Experiment journal");
  shell.append(journal);
  const future = node("section", "explore-future");
  future.append(
    node("h2", "", "Beyond the observatory"),
    node(
      "p",
      "",
      "Next destinations under design: a Feynman-diagram workshop, scattering experiment and Standard Model symmetry atlas. These are not working simulations yet.",
    ),
  );
  shell.append(future);
  const destinations = node("div", "explore-destinations");
  for (const [name, chapters, description] of [
    [
      "Feynman workshop",
      "Srednicki 9–10, 45, 72",
      "Build graphs and inspect the rule ledger behind each amplitude.",
    ],
    [
      "Scattering arena",
      "Srednicki 10–11, 48, 59, 81",
      "Vary incoming states, energy and angle; compare verified channels.",
    ],
    [
      "Symmetry atlas",
      "Srednicki 24, 69–70, 87–89",
      "Explore SU(3), SU(2), U(1), representations and symmetry breaking.",
    ],
    [
      "Field & scale gardens",
      "Srednicki 3, 6–9, 28–29",
      "Explore field modes, path histories and renormalization flows.",
    ],
  ]) {
    const card = node("article", "explore-destination");
    card.append(
      node("h3", "", name),
      node("p", "", description),
      node("p", "progress-note", chapters + " · planned destination"),
    );
    destinations.append(card);
  }
  future.append(destinations);
  // The scene is the primary interaction. Detailed instruments are opt-in overlays.
  const scanner = node("dialog", "explore-scanner");
  scanner.setAttribute("aria-label", "Scanner");
  scanner.append(
    btn("Close scan", () => scanner.close(), "overlay-close"),
    node("h2", "", "Scanner · inspect this transformation"),
    layout,
    pinRow,
  );
  shell.append(scanner);
  const journalDialog = node("dialog", "explore-journal-dialog");
  journalDialog.setAttribute("aria-label", "Field journal");
  journalDialog.append(
    btn("Close journal", () => journalDialog.close(), "overlay-close"),
    journal,
  );
  shell.append(journalDialog);
  const atlas = node("dialog", "explore-atlas");
  atlas.setAttribute("aria-label", "World atlas");
  atlas.append(
    btn("Close atlas", () => atlas.close(), "overlay-close"),
    future,
  );
  shell.append(atlas);
  const sceneHost = node("div", "explore-scene-host");
  lab.prepend(sceneHost);
  const sceneStatus = node(
    "output",
    "scene-computation",
    "Connecting to the instrument…",
  );
  sceneStatus.setAttribute("role", "status");
  lab.append(sceneStatus);
  const dock = node("nav", "scene-dock");
  dock.setAttribute("aria-label", "Observatory tools");
  const turnTool = btn("Rotate", () => chooseTool("rotate"), "scene-tool");
  const sliceTool = btn("Slice", () => chooseTool("slice"), "scene-tool");
  function chooseTool(tool: "rotate" | "slice") {
    scene?.setTool(tool);
    turnTool.setAttribute("aria-pressed", String(tool === "rotate"));
    sliceTool.setAttribute("aria-pressed", String(tool === "slice"));
    options.sound();
  }
  dockPlay = btn(
    "Play",
    () => {
      play.click();
      setDockPlayLabel(playing ? "Pause" : "Play");
    },
    "scene-tool",
  );
  const mark = btn(
    "Keep discovery",
    () => {
      if (playing) pause();
      if (save.disabled) {
        keepRequestedKey = inputKey(saved.input);
        mark.setAttribute("aria-busy", "true");
        notice.textContent =
          "Holding this point… it will be kept when its calculation finishes.";
        if (!inFlight) queue(true);
        return;
      }
      save.click();
    },
    "scene-tool",
  );
  dock.append(
    turnTool,
    sliceTool,
    dockPlay,
    btn(
      "Reset",
      () => {
        pause();
        saved.input.t = 0;
        changed(true);
      },
      "scene-tool",
    ),
    mark,
    btn(
      "Scan",
      () => {
        pause();
        openSheet(scanner);
      },
      "scene-tool",
    ),
  );
  lab.append(dock);
  function celebrate(durable: boolean) {
    clearTimeout(feedbackTimer);
    mark.classList.add("discovery-kept");
    mark.dataset.feedback = durable ? "saved" : "memory";
    feedbackTimer = window.setTimeout(() => {
      mark.classList.remove("discovery-kept");
      delete mark.dataset.feedback;
    }, 1100);
  }
  const toolInfo: [HTMLButtonElement, string, string][] = [
    [turnTool, "↻", "R"],
    [sliceTool, "▱", "S"],
    [dockPlay, "▷", "Space"],
    [mark, "✦", "K"],
  ];
  for (const [button, icon, key] of toolInfo) {
    button.setAttribute("aria-label", button.textContent || "");
    const glyph = node("span", "tool-pictogram", icon);
    glyph.setAttribute("aria-hidden", "true");
    button.prepend(glyph);
    const shortcut = node("kbd", "tool-shortcut", key === "Space" ? "␣" : key);
    shortcut.setAttribute("aria-hidden", "true");
    button.append(shortcut);
    button.setAttribute("aria-keyshortcuts", key);
  }
  function setDockPlayLabel(label: string) {
    if (!dockPlay) return;
    const text = Array.from(dockPlay.childNodes).find(
      (n) => n.nodeType === Node.TEXT_NODE,
    );
    if (text) text.textContent = label;
    else dockPlay.append(document.createTextNode(label));
    dockPlay.setAttribute("aria-label", label);
  }
  const tips = node("dialog", "explore-tips");
  tips.setAttribute("aria-label", "Observatory tips");
  const tipHeader = node("div", "tips-heading");
  tipHeader.append(
    node("span", "tips-companion", "✦"),
    node("h2", "", "Try a tiny experiment"),
  );
  tips.append(
    btn("Close tips", () => tips.close(), "overlay-close"),
    tipHeader,
    node(
      "p",
      "",
      "Drag the instrument to turn it. One complete turn reaches −I; a second brings it back to +I. The direction needle can look the same while the group element differs.",
    ),
  );
  const tipActions = node("div", "explore-actions");
  for (const [label, t] of [
    ["One turn · 2π", 2 * Math.PI],
    ["Two turns · 4π", 4 * Math.PI],
  ] as const)
    tipActions.append(
      btn(label, () => {
        tips.close();
        show("lab");
        directChange(t, "angle");
      }),
    );
  tips.append(
    tipActions,
    node(
      "p",
      "progress-note",
      "Prefer exploring freely? Close this and try Rotate, Slice, or a different axis. No unlocks, scores or time pressure.",
    ),
  );
  shell.append(tips);
  nav.append(btn("Tips", () => openSheet(tips)));

  axisCycle = btn(
    "Axis Z ↻",
    () => {
      const v = saved.input.axis;
      axis.value = v[2] === 1 ? "x" : v[0] === 1 ? "y" : "z";
      axis.dispatchEvent(new Event("change"));
    },
    "scene-axis",
  );
  axisCycle.setAttribute("aria-label", "Change rotation axis");
  lab.append(axisCycle);
  nav.append(
    btn("Field journal", () => openSheet(journalDialog)),
    btn("World atlas", () => openSheet(atlas)),
  );
  function directChange(value: number, kind: "angle" | "slice") {
    if (playing) pause();
    if (kind === "angle")
      saved.input.t = Math.max(0, Math.min(4 * Math.PI, value));
    else saved.input.slice = Math.max(-1.2, Math.min(1.2, value));
    sync();
    persist();
    save.disabled = true;
    evidence.dataset.status = "pending";
    evidence.textContent = "Updating · previous views are stale";
    scene?.setPending(true);
    sceneStatus.textContent = "Instrument updating…";
    function sendLatest() {
      if (disposed) return;
      if (inFlight) {
        clearTimeout(interactionTimer);
        interactionTimer = window.setTimeout(sendLatest, 70);
        return;
      }
      queue(true);
    }
    sendLatest();
  }
  scene = mountScene(sceneHost, {
    angle: (value) => directChange(value, "angle"),
    slice: (value) => directChange(value, "slice"),
    sound: () => options.sound(),
  });
  chooseTool("rotate");
  function openSheet(dialog: HTMLDialogElement) {
    dialog.prepend(notice);
    dialog.showModal();
  }
  const dialogEvents = new AbortController();
  sceneHost.addEventListener(
    "keydown",
    (event) => {
      if (
        !(event.target instanceof HTMLCanvasElement) ||
        event.ctrlKey ||
        event.metaKey ||
        event.altKey ||
        event.repeat
      )
        return;
      const key = event.key.toLowerCase();
      const target =
        key === "r"
          ? turnTool
          : key === "s"
            ? sliceTool
            : key === " "
              ? dockPlay
              : key === "k"
                ? mark
                : undefined;
      if (target) {
        event.preventDefault();
        target.click();
      }
    },
    { signal: dialogEvents.signal },
  );
  for (const dialog of [scanner, journalDialog, atlas, tips]) {
    dialog.addEventListener(
      "close",
      () => {
        if (!disposed) shell.insertBefore(notice, world);
      },
      { signal: dialogEvents.signal },
    );
  }
  for (const dialog of [scanner, journalDialog, atlas, tips])
    dialog.addEventListener(
      "click",
      (event) => {
        if (event.target === dialog) {
          const box = dialog.getBoundingClientRect();
          const e = event as MouseEvent;
          if (
            e.clientX < box.left ||
            e.clientX > box.right ||
            e.clientY < box.top ||
            e.clientY > box.bottom
          )
            dialog.close();
        }
      },
      { signal: dialogEvents.signal },
    );

  function persist(): boolean {
    if (blocked) {
      notice.textContent =
        "Exploration changes are in memory only. Export them or resolve saved-data recovery.";
      return false;
    }
    try {
      localStorage.setItem(KEY, JSON.stringify(saved));
      return true;
    } catch {
      notice.textContent =
        "Exploration could not be saved. Export your journal to keep it.";
      return false;
    }
  }
  function pause() {
    if (playing) {
      saved.input.t = Math.min(
        4 * Math.PI,
        startT + (performance.now() - start) / 1000,
      );
      sync(false);
    }
    playing = false;
    cancelAnimationFrame(frame);
    play.textContent = "Play";
    setDockPlayLabel("Play");
    if (!disposed) {
      persist();
      queue(true);
    }
  }
  function tick(now: number) {
    if (!playing || disposed) return;
    const next = Math.min(4 * Math.PI, startT + (now - start) / 1000);
    if (!inFlight && (now - lastRequest >= 100 || next >= 4 * Math.PI)) {
      saved.input.t = next;
      sync(false);
      lastRequest = now;
      queue(true);
    }
    if (next >= 4 * Math.PI) {
      pause();
      return;
    }
    frame = requestAnimationFrame(tick);
  }
  function sync(updateAxis = true) {
    time.value = String(saved.input.t);
    timeRead.textContent = `t = ${saved.input.t.toFixed(5)} rad · ${(saved.input.t / Math.PI).toFixed(3)}π`;
    slice.value = String(saved.input.slice);
    sliceRead.textContent = `c = ${saved.input.slice.toFixed(2)}`;
    if (!updateAxis) return;
    const v = saved.input.axis;
    if (axisCycle)
      axisCycle.textContent =
        v[0] === 1
          ? "Axis X ↻"
          : v[1] === 1
            ? "Axis Y ↻"
            : v[2] === 1
              ? "Axis Z ↻"
              : "Custom axis ↻";
    axis.value =
      v[0] === 1 && v[1] === 0 && v[2] === 0
        ? "x"
        : v[1] === 1 && v[0] === 0 && v[2] === 0
          ? "y"
          : v[2] === 1 && v[0] === 0 && v[1] === 0
            ? "z"
            : "custom";
    custom.hidden = axis.value !== "custom";
    components.forEach((e, i) => (e.value = String(v[i])));
  }
  function changed(immediate = false) {
    sync();
    persist();
    queue(immediate);
  }
  function queue(immediate = false) {
    clearTimeout(timer);
    abort?.abort();
    sequence++;
    evidence.textContent = "Updating · previous views are stale";
    evidence.dataset.status = "pending";
    save.disabled = true;
    scene?.setPending(true);
    sceneStatus.textContent = "Instrument updating…";
    if (immediate) void compute(sequence);
    else timer = window.setTimeout(() => void compute(sequence), 60);
  }
  async function compute(token: number) {
    if (disposed) return;
    inFlight = true;
    abort = new AbortController();
    const input = clone(saved.input);
    try {
      const response = await fetch("/api/su2", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(input),
        signal: abort.signal,
      });
      const result = await response.json();
      if (!response.ok)
        throw Error(result.error?.message || "Rust calculation failed");
      if (disposed || token !== sequence) return;
      snapshot = result as Snapshot;
      snapshotKey = inputKey(input);
      if (snapshotKey !== inputKey(saved.input)) return;
      plots.update(snapshot);
      scene?.update(snapshot);
      scene?.setPending(false);
      sceneStatus.textContent = `t ${(input.t / Math.PI).toFixed(2)}π · Numerically checked`;
      matrix.textContent = matrixText(snapshot);
      evidence.textContent = `Numerically checked · Rust · unitary error ${snapshot.checks.unitarity_error.toExponential(1)} · determinant error ${snapshot.checks.determinant_error.toExponential(1)}`;
      evidence.dataset.status = "checked";
      evidence.dataset.t = String(input.t);
      save.disabled = playing;
      renderComparison();
      if (keepRequestedKey) {
        const wanted = keepRequestedKey;
        keepRequestedKey = undefined;
        mark.removeAttribute("aria-busy");
        if (wanted === snapshotKey && !playing) save.click();
        else
          notice.textContent =
            "The point changed before it was kept. Choose Keep discovery again for the new point.";
      }
    } catch (e) {
      if (disposed || token !== sequence || (e as Error).name === "AbortError")
        return;
      if (keepRequestedKey)
        notice.textContent =
          "Could not keep this discovery: " + (e as Error).message;
      keepRequestedKey = undefined;
      mark.removeAttribute("aria-busy");
      evidence.dataset.status = "error";
      scene?.setPending(true);
      sceneStatus.textContent =
        "Instrument unavailable · retry with Reset or Scan";
      evidence.textContent = "Calculation unavailable: " + (e as Error).message;
      matrix.textContent =
        "No current result. Previous plot, if any, is stale.";
    } finally {
      if (token === sequence) inFlight = false;
    }
  }
  function show(location: "world" | "lab") {
    for (const dialog of [scanner, journalDialog, atlas, tips])
      if (dialog.open) dialog.close();
    shell.dataset.location = location;
    if (location === "world" && playing) pause();
    saved.location = location;
    world.hidden = location !== "world";
    lab.hidden = location !== "lab";
    worldButton.setAttribute("aria-pressed", String(location === "world"));
    labButton.setAttribute("aria-pressed", String(location === "lab"));
    persist();
    if (location === "lab") queue(true);
    options.sound();
  }
  const comparison = node("p", "explore-comparison");
  function renderComparison() {
    if (!compare) return;
    comparison.textContent = `Comparing with “${compare.title}”: Δt = ${(saved.input.t - compare.input.t).toFixed(5)} rad; previous axis (${compare.input.axis.join(", ")}), slice ${compare.input.slice}. Parameter comparison; this is not a group distance.`;
  }
  function renderJournal() {
    journal.replaceChildren(node("h2", "", "Field journal"));
    const tools = node("div", "explore-actions");
    tools.append(
      btn("Export exploration", () =>
        download(JSON.stringify(saved, null, 2), "symmetry-exploration.json"),
      ),
    );
    const importInput = node("input");
    importInput.type = "file";
    importInput.accept = ".json";
    importInput.hidden = true;
    importInput.setAttribute("aria-label", "Import exploration file");
    importInput.onchange = async () => {
      const file = importInput.files?.[0];
      if (!file) return;
      try {
        if (file.size > 2_000_000) throw Error("Import exceeds 2 MB.");
        const text = await file.text();
        if (disposed) return;
        const imported = decode(text);
        if (
          !window.confirm(
            "Replace exploration inputs and journal with this file? Other applications are unchanged.",
          )
        )
          return;
        keepRequestedKey = undefined;
        mark.removeAttribute("aria-busy");
        compare = undefined;
        comparison.textContent = "";
        pause();
        saved = imported;
        blocked = false;
        recovery.remove();
        sync();
        persist();
        renderJournal();
        show(saved.location);
      } catch (e) {
        notice.textContent = "Import rejected: " + (e as Error).message;
      }
    };
    tools.append(
      btn("Import exploration", () => importInput.click()),
      importInput,
    );
    journal.append(tools, comparison);
    if (!saved.journal.length)
      journal.append(
        node(
          "p",
          "progress-note",
          "Your experiments will appear here. Record any point; no unlocks or reaction-speed challenges.",
        ),
      );
    for (const entry of saved.journal) {
      const card = node("article", "explore-journal-entry");
      card.append(
        node("h3", "", entry.title),
        node(
          "p",
          "",
          `t=${entry.input.t.toFixed(5)} rad · n=(${entry.input.axis.map((v) => v.toFixed(3)).join(", ")}) · slice ${entry.input.slice}`,
        ),
        node("p", "", entry.note),
      );
      const actions = node("div", "explore-actions");
      actions.append(
        btn("Reopen", () => {
          pause();
          saved.input = clone(entry.input);
          sync();
          show("lab");
        }),
        btn("Compare", () => {
          compare = entry;
          renderComparison();
        }),
        btn("Export discovery", () =>
          download(
            JSON.stringify(
              {
                format: "mathphysics.discovery.v1",
                convention: "pauli-plus-i-v1",
                ...entry,
              },
              null,
              2,
            ),
            "discovery.json",
          ),
        ),
      );
      if (entry.source && options.openSource)
        actions.append(
          btn("Open source: " + entry.source.name, () =>
            options.openSource?.(entry.source!),
          ),
        );
      actions.append(
        btn("Remove", () => {
          saved.journal = saved.journal.filter((e) => e.id !== entry.id);
          if (compare?.id === entry.id) {
            compare = undefined;
            comparison.textContent = "";
          }
          persist();
          renderJournal();
        }),
      );
      card.append(actions);
      journal.append(card);
    }
  }
  axis.onchange = () => {
    pause();
    if (axis.value === "custom") {
      custom.hidden = false;
      return;
    }
    saved.input.axis =
      axis.value === "x"
        ? [1, 0, 0]
        : axis.value === "y"
          ? [0, 1, 0]
          : [0, 0, 1];
    changed(true);
  };
  time.oninput = () => {
    const selectedAngle = Number(time.value);
    pause();
    saved.input.t = selectedAngle;
    changed();
  };
  slice.oninput = () => {
    saved.input.slice = Number(slice.value);
    changed();
  };
  const visibility = () => {
    if (document.hidden && playing) pause();
  };
  document.addEventListener("visibilitychange", visibility);
  sync();
  renderJournal();
  show(saved.location);
  return {
    exportState() {
      download(JSON.stringify(saved, null, 2), "symmetry-exploration.json");
    },
    dispose() {
      disposed = true;
      playing = false;
      cancelAnimationFrame(frame);
      clearTimeout(timer);
      abort?.abort();
      document.removeEventListener("visibilitychange", visibility);
      clearTimeout(interactionTimer);
      clearTimeout(feedbackTimer);
      scene?.dispose();
      dialogEvents.abort();
      for (const dialog of [scanner, journalDialog, atlas, tips])
        if (dialog.open) dialog.close();
      worldView.dispose();
      plots.dispose();
      persist();
      shell.remove();
    },
  };
}
