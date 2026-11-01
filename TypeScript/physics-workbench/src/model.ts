export type Style = "notebook" | "pixel";
export type Kind = "directory" | "pdf" | "tex" | "markdown" | "text" | "file";
export interface Root {
  id: string;
  label: string;
  path: string;
}
export interface Entry {
  name: string;
  path: string;
  kind: Kind;
  size: number;
  modified_secs: number;
}
export interface Place {
  root: string;
  path: string;
  name: string;
}
export interface FileTab {
  id: string;
  root: string;
  path: string;
  name: string;
  kind: Kind;
  text?: string;
  paired?: { path: string; kind: Kind } | null;
  draft?: string;
  page?: number;
}
export interface Cell {
  id: string;
  kind: "math" | "text" | "trace";
  source: string;
  output?: unknown;
}
export interface Workspace {
  version: 1;
  style: Style;
  sound: boolean;
  volume: number;
  cells: Cell[];
  tabs: FileTab[];
  activeTab: string;
  favorites: Place[];
  recent: Place[];
  visits: Record<string, number>;
  location: Place;
  zoom: number;
}
const key = "mathphysics.workbench.v1";
export function initialWorkspace(): Workspace {
  return {
    version: 1,
    style: "notebook",
    sound: false,
    volume: 0.18,
    activeTab: "notebook",
    tabs: [],
    favorites: [],
    recent: [],
    visits: {},
    location: { root: "mathphysics", path: "", name: "Workspace" },
    zoom: 100,
    cells: [
      {
        id: "intro",
        kind: "text",
        source:
          "## A place to think clearly\nKeep the calculation, its assumptions, and the source together.",
      },
      {
        id: "conventions",
        kind: "math",
        source:
          "g_{\\mu\\nu}=\\operatorname{diag}(-1,+1,+1,+1),\\qquad \\{\\gamma^\\mu,\\gamma^\\nu\\}=-2g^{\\mu\\nu}I",
      },
      { id: "first-trace", kind: "trace", source: "mu, nu, rho, sigma" },
    ],
  };
}
function validPlace(p: unknown): p is Place {
  return (
    !!p &&
    typeof p === "object" &&
    ["root", "path", "name"].every(
      (k) => typeof (p as Record<string, unknown>)[k] === "string",
    )
  );
}
export function restoreWorkspace(): Workspace {
  const base = initialWorkspace();
  try {
    const raw = localStorage.getItem(key);
    if (!raw) return base;
    const s = JSON.parse(raw) as Partial<Workspace>;
    if (
      s.version !== 1 ||
      !Array.isArray(s.cells) ||
      !s.cells.every(
        (c) =>
          c &&
          ["math", "text", "trace"].includes(c.kind) &&
          typeof c.source === "string" &&
          typeof c.id === "string",
      )
    )
      return base;
    return {
      ...base,
      ...s,
      style: s.style === "pixel" ? "pixel" : "notebook",
      sound: s.sound === true,
      volume:
        typeof s.volume === "number" && Number.isFinite(s.volume)
          ? Math.max(0, Math.min(1, s.volume))
          : 0.18,
      cells: s.cells.slice(0, 200),
      tabs: Array.isArray(s.tabs)
        ? s.tabs
            .filter(
              (t) =>
                t &&
                typeof t.id === "string" &&
                typeof t.path === "string" &&
                typeof t.root === "string" &&
                typeof t.name === "string",
            )
            .slice(0, 30)
        : [],
      favorites: Array.isArray(s.favorites)
        ? s.favorites.filter(validPlace).slice(0, 100)
        : [],
      recent: Array.isArray(s.recent)
        ? s.recent.filter(validPlace).slice(0, 7)
        : [],
      visits: s.visits && typeof s.visits === "object" ? s.visits : {},
      location: validPlace(s.location) ? s.location : base.location,
      zoom: [100, 125, 150, 200].includes(s.zoom ?? 0) ? s.zoom! : 100,
    };
  } catch {
    return base;
  }
}
export function saveWorkspace(state: Workspace): boolean {
  try {
    localStorage.setItem(key, JSON.stringify(state));
    return true;
  } catch {
    return false;
  }
}
export function rememberPlace(state: Workspace, place: Place): void {
  const id = place.root + ":" + place.path;
  state.visits[id] = (Number(state.visits[id]) || 0) + 1;
  state.recent = [
    place,
    ...state.recent.filter(
      (p) => p.root !== place.root || p.path !== place.path,
    ),
  ].slice(0, 7);
}
