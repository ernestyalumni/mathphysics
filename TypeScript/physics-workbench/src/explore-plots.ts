/** Rust owns mathematical state; this module owns drawing geometry only. */
export interface PlotSnapshot {
  coordinates: [number, number, number, number];
  alpha: [number, number];
  beta: [number, number];
  stereographic: { status: "finite" | "pole"; point: number[] | null };
  slice: {
    coordinate: number;
    radius: number | null;
    status: "sphere" | "point" | "empty";
  };
  bloch: number[];
}
type V = [number, number, number];
const f = (n: number): string =>
  !Number.isFinite(n)
    ? "invalid"
    : n !== 0 && (Math.abs(n) < 0.0001 || Math.abs(n) > 1e5)
      ? n.toExponential(4)
      : Number(n.toFixed(5)).toString();
const vec = (a: number[]): string => `(${a.map(f).join(", ")})`;
const finite = (a: number[] | null): a is V =>
  !!a && a.length === 3 && a.every(Number.isFinite);
export function mountPlots(host: HTMLElement): {
  update(snapshot: PlotSnapshot): void;
  dispose(): void;
} {
  const d = host.ownerDocument,
    w = d.defaultView!;
  const el = <K extends keyof HTMLElementTagNameMap>(
    tag: K,
    text: string,
    parent: HTMLElement,
    css = "",
  ): HTMLElementTagNameMap[K] => {
    const a = d.createElement(tag);
    a.textContent = text;
    a.style.cssText = css;
    parent.append(a);
    return a;
  };
  const root = el(
    "section",
    "",
    host,
    "display:grid;grid-template-columns:minmax(0,1fr);gap:16px;color:inherit;font:inherit;min-width:0",
  );
  root.className = "explore-plots";
  root.setAttribute("aria-label", "Linked SU(2) coordinate observatory");
  el("h3", "SU(2) · Coordinate observatory", root, "margin:0");
  el(
    "p",
    "One state, linked views · α = a₀ + i a₃ · β = a₂ + i a₁ · |α|² + |β|² = 1",
    root,
    "margin:0",
  );
  const summary = el("output", "Awaiting Rust snapshot", root);
  const grid = el(
    "div",
    "",
    root,
    "display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,310px),1fr));gap:16px",
  );
  const panels = [
    ["01 / α complex plane", "Re α = a₀ · Im α = a₃"],
    ["02 / β complex plane", "Re β = a₂ · Im β = a₁"],
    [
      "03 / S³ stereographic chart → ℝ³",
      "x = (a₁, a₂, a₃)/(1 − a₀) · missing pole I = (1, 0, 0, 0)",
    ],
    [
      "04 / S³ coordinate slice",
      "a₀ = c · a₁² + a₂² + a₃² = 1 − c² · radius √(1 − c²)",
    ],
  ].map(([title, subtitle]) => {
    const card = el(
      "section",
      "",
      grid,
      "padding:16px;border:1px solid currentColor;border-radius:14px;min-width:0;background:color-mix(in srgb,currentColor 3%,transparent)",
    );
    card.className = "plot-card";
    el("h4", title, card, "margin:0 0 8px;font:inherit;font-weight:700");
    el("p", subtitle, card, "margin:0 0 8px;font-size:.9em;line-height:1.5");
    const canvas = el(
      "canvas",
      "",
      card,
      "display:block;width:100%;height:260px;border-radius:8px",
    );
    canvas.className = "plot-canvas";
    canvas.setAttribute("role", "img");
    canvas.setAttribute("aria-label", title);
    const output = el(
      "output",
      "Awaiting Rust snapshot",
      card,
      "display:block;line-height:1.6;overflow-wrap:anywhere;font-variant-numeric:tabular-nums",
    );
    output.className = "plot-readout";
    return { canvas, output };
  });
  const controls = el("div", "", root, "display:flex;flex-wrap:wrap;gap:8px");
  controls.setAttribute("role", "group");
  controls.setAttribute("aria-label", "Shared 3D camera");
  const camera = el("output", "", root);
  el(
    "p",
    "Diamond = selected state · dashed circle = unit reference. Canvas 2D orthographic 3D prototype. Drag either 3D view or use arrows; Home resets. No automatic motion.",
    root,
    "margin:0;font-size:.9em;line-height:1.5",
  );
  const bloch = el(
    "output",
    "Bloch S²: awaiting Rust snapshot",
    root,
    "padding:12px;border-left:3px solid currentColor;overflow-wrap:anywhere",
  );
  let state: PlotSnapshot | undefined,
    dead = false,
    yaw = 0.62,
    pitch = -0.38;
  const cleanup: (() => void)[] = [];
  function on(t: EventTarget, name: string, fn: EventListener): void {
    t.addEventListener(name, fn);
    cleanup.push(() => t.removeEventListener(name, fn));
  }
  function orbit(x: number, y: number): void {
    yaw = (yaw + x) % (2 * Math.PI);
    pitch = Math.max(-1.45, Math.min(1.45, pitch + y));
    draw();
  }
  function reset(): void {
    yaw = 0.62;
    pitch = -0.38;
    draw();
  }
  for (const [label, action] of [
    ["Orbit left", () => orbit(-0.18, 0)],
    ["Orbit right", () => orbit(0.18, 0)],
    ["Orbit up", () => orbit(0, -0.18)],
    ["Orbit down", () => orbit(0, 0.18)],
    ["Reset camera", reset],
  ] as const) {
    const b = el(
      "button",
      label,
      controls,
      "font:inherit;color:inherit;background:transparent;border:1px solid currentColor;border-radius:7px;padding:8px 12px;min-height:44px;cursor:pointer",
    );
    b.type = "button";
    on(b, "click", action);
  }
  for (const { canvas: c } of panels.slice(2)) {
    c.tabIndex = 0;
    c.style.touchAction = "none";
    c.style.cursor = "grab";
    c.setAttribute(
      "aria-label",
      `${c.getAttribute("aria-label")}. Arrow keys orbit; Home resets.`,
    );
    let drag: { id: number; x: number; y: number } | undefined;
    const end = (): void => {
      const old = drag;
      drag = undefined;
      if (old && c.hasPointerCapture(old.id)) c.releasePointerCapture(old.id);
      c.style.cursor = "grab";
    };
    on(c, "pointerdown", ((e: PointerEvent) => {
      if (e.button || drag) return;
      drag = { id: e.pointerId, x: e.clientX, y: e.clientY };
      c.setPointerCapture(e.pointerId);
      c.focus();
      c.style.cursor = "grabbing";
    }) as EventListener);
    on(c, "pointermove", ((e: PointerEvent) => {
      if (!drag || drag.id !== e.pointerId) return;
      orbit((e.clientX - drag.x) * 0.008, (e.clientY - drag.y) * 0.008);
      drag.x = e.clientX;
      drag.y = e.clientY;
    }) as EventListener);
    for (const name of ["pointerup", "pointercancel", "lostpointercapture"])
      on(c, name, end);
    cleanup.push(end);
    on(c, "keydown", ((e: KeyboardEvent) => {
      const step = e.shiftKey ? 0.35 : 0.12;
      const keys: Record<string, [number, number]> = {
        ArrowLeft: [-step, 0],
        ArrowRight: [step, 0],
        ArrowUp: [0, -step],
        ArrowDown: [0, step],
      };
      if (keys[e.key]) {
        e.preventDefault();
        orbit(...keys[e.key]);
      } else if (e.key === "Home") {
        e.preventDefault();
        reset();
      }
    }) as EventListener);
  }
  function diamond(c: CanvasRenderingContext2D, x: number, y: number): void {
    c.save();
    c.setLineDash([]);
    c.lineWidth = 2;
    c.beginPath();
    c.moveTo(x, y - 6);
    c.lineTo(x + 6, y);
    c.lineTo(x, y + 6);
    c.lineTo(x - 6, y);
    c.closePath();
    c.fill();
    c.stroke();
    c.beginPath();
    c.arc(x, y, 11, 0, Math.PI * 2);
    c.stroke();
    c.restore();
  }
  function paint(index: number, values?: number[]): void {
    const { canvas, output } = panels[index];
    const width = Math.max(
        1,
        Math.min(1600, canvas.getBoundingClientRect().width),
      ),
      height = 260,
      ratio = Math.min(2, w.devicePixelRatio || 1);
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    const c = canvas.getContext("2d");
    if (!c) return;
    c.setTransform(ratio, 0, 0, ratio, 0, 0);
    c.strokeStyle = c.fillStyle = w.getComputedStyle(root).color;
    c.font = "12px system-ui,sans-serif";
    const size = Math.max(1, Math.min(width - 80, height - 66) / 2);
    const line = (a: number[], b: number[]): void => {
      c.beginPath();
      c.moveTo(a[0], a[1]);
      c.lineTo(b[0], b[1]);
      c.stroke();
    };
    if (index < 2) {
      const valid = values?.length === 2 && values.every(Number.isFinite);
      const extent = valid ? Math.max(1.2, ...values!.map(Math.abs)) : 1.2;
      const xy = (x: number, y: number): number[] => [
        width / 2 + (x / extent) * size,
        height / 2 - (y / extent) * size,
      ];
      c.globalAlpha = 0.12;
      for (let i = -4; i <= 4; i++) {
        const v = (extent * i) / 4;
        line(xy(v, -extent), xy(v, extent));
        line(xy(-extent, v), xy(extent, v));
      }
      c.globalAlpha = 0.55;
      c.setLineDash([4, 4]);
      c.beginPath();
      c.arc(width / 2, height / 2, size / extent, 0, Math.PI * 2);
      c.stroke();
      c.setLineDash([]);
      c.globalAlpha = 1;
      line(xy(-extent, 0), xy(extent, 0));
      line(xy(0, -extent), xy(0, extent));
      c.fillText("Re", width / 2 + size - 10, height / 2 + 20);
      c.fillText("Im", width / 2 + 8, height / 2 - size);
      c.fillText(f(-extent), width / 2 - size, height / 2 - 8);
      c.fillText(f(extent), width / 2 + size - 20, height / 2 - 8);
      c.fillText("0", width / 2 + 4, height / 2 + 14);
      c.fillText(f(extent), width / 2 + 8, height / 2 - size + 16);
      c.fillText(f(-extent), width / 2 + 8, height / 2 + size);
      if (valid) {
        const p = xy(values![0], values![1]);
        c.setLineDash([3, 4]);
        line(xy(values![0], 0), p);
        line(xy(0, values![1]), p);
        c.setLineDash([]);
        line(xy(0, 0), p);
        diamond(c, p[0], p[1]);
      }
      return;
    }
    const s = state;
    let point: V | null = null,
      radius: number | null = index === 2 ? 1 : null,
      extent = 1.25;
    if (s && index === 2) {
      if (s.stereographic.status === "finite" && finite(s.stereographic.point))
        point = s.stereographic.point;
      if (point) extent = Math.max(1.25, ...point.map(Math.abs));
    } else if (s) {
      const sl = s.slice;
      const valid =
        Number.isFinite(sl.coordinate) &&
        sl.radius !== null &&
        Number.isFinite(sl.radius) &&
        ((sl.status === "sphere" && sl.radius > 0) ||
          (sl.status === "point" && sl.radius === 0));
      if (valid) {
        radius = sl.radius;
        if (
          s.coordinates.every(Number.isFinite) &&
          Math.abs(s.coordinates[0] - sl.coordinate) <= 1e-7
        )
          point = s.coordinates.slice(1) as V;
        extent = Math.max(1.25, radius! * 1.2, ...(point?.map(Math.abs) ?? []));
      }
    }
    const rotate = (a: V): V => {
      const x = a[0] * Math.cos(yaw) + a[2] * Math.sin(yaw),
        z = -a[0] * Math.sin(yaw) + a[2] * Math.cos(yaw);
      return [
        x,
        a[1] * Math.cos(pitch) - z * Math.sin(pitch),
        a[1] * Math.sin(pitch) + z * Math.cos(pitch),
      ];
    };
    const screen = (a: V): number[] => {
      const r = rotate(a);
      return [width / 2 + r[0] * size, height / 2 - r[1] * size];
    };
    if (radius !== null) {
      if (radius === 0) diamond(c, width / 2, height / 2);
      else
        for (let j = 0; j < 11; j++) {
          let previous: V | undefined;
          for (let i = 0; i <= 64; i++) {
            const t = (i * Math.PI) / 32,
              a = j < 6 ? (j * Math.PI) / 6 : ((j - 8) * Math.PI) / 6,
              r = radius / extent;
            const p: V =
              j < 6
                ? [
                    r * Math.cos(t) * Math.cos(a),
                    r * Math.sin(t),
                    r * Math.cos(t) * Math.sin(a),
                  ]
                : [
                    r * Math.cos(a) * Math.cos(t),
                    r * Math.sin(a),
                    r * Math.cos(a) * Math.sin(t),
                  ];
            const back = rotate(p)[2] < 0;
            c.globalAlpha = back ? 0.14 : 0.5;
            c.setLineDash(back || index === 2 ? [3, 4] : []);
            if (previous) line(screen(previous), screen(p));
            previous = p;
          }
        }
    }
    c.globalAlpha = 0.75;
    c.setLineDash([]);
    for (let axis = 0; axis < 3; axis++) {
      const a: V = [0, 0, 0],
        b: V = [0, 0, 0];
      a[axis] = -1;
      b[axis] = 1;
      line(screen(a), screen(b));
      const p = screen(b);
      c.fillText(
        `+${index === 2 ? "x" : "a"}${axis + 1}`,
        Math.max(2, Math.min(width - 35, p[0] + 5)),
        Math.max(14, Math.min(height - 15, p[1] - 5)),
      );
    }
    c.globalAlpha = 1;
    if (point) {
      const p = screen(point.map((n) => n / extent) as V);
      c.setLineDash([3, 4]);
      line([width / 2, height / 2], p);
      c.setLineDash([]);
      diamond(c, p[0], p[1]);
    }
    c.fillText(`Axis extent ±${f(extent)} · orthographic`, 8, height - 8);
    if (s && index === 2)
      output.textContent =
        s.stereographic.status === "pole"
          ? "Pole / off-chart: no finite point supplied by Rust (projection pole or numerical chart limit)."
          : point
            ? `x = ${vec(point)} · fitted extent ±${f(extent)} · unit reference sphere, not S³`
            : "Invalid finite chart point; no selected point drawn.";
    if (s && index === 3)
      output.textContent =
        `c = ${f(s.slice.coordinate)} · ` +
        (s.slice.status === "empty"
          ? "Empty slice."
          : radius === null
            ? "Invalid slice geometry."
            : `${s.slice.status === "point" ? "Degenerate point" : "Sphere"} · radius ${f(radius)} · ${point ? `selected ${vec(point)}` : "selected state off-slice"} · display tolerance |a₀ − c| ≤ 10⁻⁷`);
  }
  function draw(): void {
    if (dead) return;
    camera.textContent = `Shared camera · azimuth ${f((yaw * 180) / Math.PI)}° · elevation ${f((pitch * 180) / Math.PI)}°`;
    if (state) {
      summary.textContent = `SU(2) ≅ S³ · (a₀, a₁, a₂, a₃) = ${vec(state.coordinates)} · U and −U remain distinct`;
      panels[0].output.textContent = `α (Re, Im) = ${vec(state.alpha)}`;
      panels[1].output.textContent = `β (Re, Im) = ${vec(state.beta)}`;
      bloch.textContent = `Bloch S² vector (Rust): ${vec(state.bloch)}. This differs from SU(2) ≅ S³ and from the a₀ = c slice. U and −U have the same Bloch action.`;
    }
    paint(0, state?.alpha);
    paint(1, state?.beta);
    paint(2);
    paint(3);
    for (const p of panels)
      p.canvas.setAttribute("aria-description", p.output.textContent || "");
  }
  const resize = new w.ResizeObserver(draw);
  resize.observe(root);
  const themes = new w.MutationObserver(draw);
  for (let a: HTMLElement | null = host; a; a = a.parentElement)
    themes.observe(a, {
      attributes: true,
      attributeFilter: ["class", "style", "data-theme", "data-style"],
    });
  on(w, "resize", draw);
  on(w.matchMedia("(prefers-color-scheme: dark)"), "change", draw);
  draw();
  return {
    update(s): void {
      if (dead) return;
      state = {
        coordinates: [...s.coordinates],
        alpha: [...s.alpha],
        beta: [...s.beta],
        stereographic: {
          ...s.stereographic,
          point: s.stereographic.point ? [...s.stereographic.point] : null,
        },
        slice: { ...s.slice },
        bloch: [...s.bloch],
      };
      draw();
    },
    dispose(): void {
      if (dead) return;
      dead = true;
      resize.disconnect();
      themes.disconnect();
      cleanup.forEach((fn) => fn());
      root.remove();
      state = undefined;
    },
  };
}
