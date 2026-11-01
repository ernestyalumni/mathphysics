/** Original procedural observatory. Mathematical state is supplied by Rust. */
export interface SceneState {
  t: number;
  axis: [number, number, number];
  coordinates: [number, number, number, number];
  bloch: number[];
  slice: {
    coordinate: number;
    radius: number | null;
    status: "sphere" | "point" | "empty";
  };
}
type Tool = "rotate" | "slice";
export function mountScene(
  host: HTMLElement,
  actions: {
    angle(t: number): void;
    slice(c: number): void;
    sound(): void;
  },
): {
  update(state: SceneState): void;
  setPending(pending: boolean): void;
  setTool(tool: Tool): void;
  dispose(): void;
} {
  const root = document.createElement("div");
  root.className = "observatory-scene";
  root.dataset.grip = "ready";
  root.style.cssText =
    "position:relative;width:100%;height:100%;min-height:0;overflow:hidden;background:#090f21;color:#e1f4ef;font:13px/1.5 monospace";
  const canvas = document.createElement("canvas");
  canvas.tabIndex = 0;
  canvas.setAttribute("role", "slider");
  canvas.style.cssText =
    "display:block;width:100%;height:100%;touch-action:none;outline-offset:-4px;cursor:grab";
  const readout = document.createElement("div");
  readout.style.cssText =
    "position:absolute;left:14px;top:12px;padding:5px 9px;background:#080e20d9;pointer-events:none;max-width:calc(100% - 28px)";
  const status = document.createElement("span");
  status.setAttribute("role", "status");
  const summary = document.createElement("div");
  const values = document.createElement("div");
  // Full numerical readout remains available to assistive technology without a panel.
  values.style.cssText =
    "position:absolute;width:1px;height:1px;padding:0;overflow:hidden;clip-path:inset(50%);white-space:nowrap";
  const legend = document.createElement("div");
  legend.textContent = "Needle: Bloch direction · shape: a₀ hemisphere";
  legend.style.cssText = "font-size:11px;color:#b9d4ce";
  readout.append(status, summary, legend, values);
  const hint = document.createElement("div");
  hint.style.cssText =
    "position:absolute;bottom:12px;left:12px;right:12px;text-align:center;background:#080e20d9;padding:5px;pointer-events:none;font-size:clamp(11px,1.4vw,14px)";
  root.append(canvas, readout, hint);
  host.append(root);
  const ctx = canvas.getContext("2d"),
    events = new AbortController();
  const motion = matchMedia("(prefers-reduced-motion: reduce)");
  let state: SceneState | undefined,
    tool: Tool = "rotate",
    pending = false,
    dead = false,
    inside = false;
  let w = 0,
    h = 0,
    frame = 0,
    angle = 0,
    slice = 0;
  let hovered = false;
  let focused = false;
  let drag: { id: number; x: number; y: number; value: number } | undefined;
  const layout = () => ({
    x: w * 0.5,
    y: h * 0.51,
    r: Math.max(1, Math.min(w * 0.29, h * 0.3)),
  });
  const visible = () => !dead && inside && !document.hidden && w > 0 && h > 0;
  const fmt = (n: number) => (Number.isFinite(n) ? n.toFixed(3) : "—");
  function describe() {
    status.textContent = pending
      ? "Updating"
      : state
        ? "SU(2) · Observatory"
        : "Awaiting Rust state";
    summary.textContent = state
      ? tool === "rotate"
        ? `t ${fmt(state.t)} rad`
        : `Slice ${fmt(state.slice.coordinate)} · ${state.slice.status}`
      : "Waiting for signal";
    values.textContent = state
      ? `t ${fmt(state.t)} radians; axis (${state.axis.map(fmt)}); U coefficients (${state.coordinates.map(fmt)}); Bloch (${state.bloch.map(fmt)}); slice coordinate ${fmt(state.slice.coordinate)}, radius ${state.slice.radius === null ? "none" : fmt(state.slice.radius)}, ${state.slice.status}`
      : "No mathematical state available";
    summary.style.opacity = pending ? ".4" : "1";
    values.style.opacity = pending ? ".4" : "1";
    root.setAttribute("aria-busy", String(pending));
    hint.textContent =
      tool === "rotate"
        ? "Drag the apparatus ↔ • ← → rotate • Home / End checkpoints"
        : "Drag the plane ↔ / ↕ scan • ↑ ↓ move slice";
    canvas.setAttribute("aria-valuemin", tool === "rotate" ? "0" : "-1.2");
    canvas.setAttribute(
      "aria-valuemax",
      String(tool === "rotate" ? 4 * Math.PI : 1.2),
    );
    canvas.setAttribute(
      "aria-valuenow",
      String(tool === "rotate" ? angle : slice),
    );
    canvas.setAttribute(
      "aria-label",
      `SU(2) observatory ${tool === "rotate" ? "rotation ring" : "slice scanner"}. ${hint.textContent}`,
    );
  }
  function draw(now = 0) {
    if (!ctx || !visible()) return;
    const dpr = Math.min(devicePixelRatio || 1, 2);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.imageSmoothingEnabled = false;
    const box = (x: number, y: number, a: number, b: number, c: string) => {
      ctx.fillStyle = c;
      ctx.fillRect(Math.round(x), Math.round(y), Math.ceil(a), Math.ceil(b));
    };
    const line = (
      x: number,
      y: number,
      a: number,
      b: number,
      c: string,
      k = 2,
    ) => {
      ctx.strokeStyle = c;
      ctx.lineWidth = k;
      ctx.beginPath();
      ctx.moveTo(x, y);
      ctx.lineTo(a, b);
      ctx.stroke();
    };
    const oval = (
      x: number,
      y: number,
      a: number,
      b: number,
      c: string,
      k = 2,
    ) => {
      ctx.strokeStyle = c;
      ctx.lineWidth = k;
      ctx.beginPath();
      ctx.ellipse(x, y, Math.max(0.1, a), Math.max(0.1, b), 0, 0, Math.PI * 2);
      ctx.stroke();
    };
    const text = (s: string, x: number, y: number, c = "#c7e5dd") => {
      ctx.fillStyle = c;
      ctx.font = "13px monospace";
      ctx.textAlign = "center";
      ctx.fillText(s, x, y);
    };
    const { x, y, r } = layout();
    box(0, 0, w, h, "#0a1122");
    const sky = ctx.createRadialGradient(x, h * 0.3, 0, x, h * 0.3, w * 0.65);
    sky.addColorStop(0, "#293a58");
    sky.addColorStop(1, "#090f21");
    ctx.fillStyle = sky;
    ctx.fillRect(0, 0, w, h * 0.8);
    for (let i = 0; i < 110; i++) {
      ctx.globalAlpha = motion.matches
        ? 1
        : 0.55 + 0.45 * Math.sin(now / 1700 + i);
      box(
        (((i * 137 + 31) % 997) / 997) * w,
        (((i * 89 + 17) % 541) / 541) * h * 0.78,
        i % 11 ? 1 : 3,
        i % 11 ? 1 : 3,
        i % 3 ? "#8aacca" : "#f8dba3",
      );
    }
    ctx.globalAlpha = 1;
    // Original dusk conservatory: distant landscape, timber, glazed copper and plants.
    for (let i = 0; i < 18; i++) {
      const tx = (i * w) / 17;
      const height = h * (0.1 + ((i * 7) % 11) / 90);
      box(tx, h * 0.76 - height, w / 16, height, i % 2 ? "#324e50" : "#263f49");
      box(tx + 5, h * 0.76 - height - 8, w / 24, 12, "#3f605b");
    }
    for (const px of [w * 0.055, w * 0.9]) {
      box(px - 5, 0, w * 0.055 + 10, h * 0.84, "#332c33");
      box(px, 0, w * 0.045, h * 0.84, "#87573d");
      box(px + 4, 0, 5, h * 0.84, "#c38a57");
      for (let i = 0; i < 9; i++) {
        const yy = (i * h) / 11;
        box(px, yy, w * 0.045, 10, "#bd8051");
        box(px + 6, yy + 3, 3, 3, "#f3ce86");
        box(px + 12, yy + 20, 2, h / 24, "#65402f");
      }
    }
    for (let i = 0; i < 11; i++) {
      const bx = w * (0.1 + i * 0.073),
        by = Math.abs(i - 5) * h * 0.018;
      box(bx, by, w * 0.077, 25, "#553b36");
      box(bx, by, w * 0.077, 5, "#ce995f");
      box(bx + 5, by + 9, w * 0.06, 3, "#98633f");
    }
    for (const side of [-1, 1]) {
      const px = x + side * r * 1.42;
      box(px - 30, y + r * 0.42, 60, 9, "#d4a16e");
      box(px - 25, y + r * 0.42 + 9, 50, 12, "#634737");
      box(px - 14, y + r * 0.24, 28, r * 0.18, "#b86b4a");
      box(px - 17, y + r * 0.23, 34, 6, "#e0a176");
      for (let i = 0; i < 12; i++) {
        const leafX = px + Math.sin(i * 2.3) * 22;
        const leafY = y + r * 0.23 - i * 6;
        line(px, leafY + 10, leafX, leafY, "#678b57", 3);
        box(leafX - 4, leafY - 3, 10, 6, i % 3 ? "#83a766" : "#c4b46f");
      }
      box(px - 8, y + r * 0.57, 16, r * 0.5, "#84583d");
      box(px - 4, y + r * 0.57, 4, r * 0.5, "#c68b53");
    }
    box(0, h * 0.82, w, h * 0.18, "#453c3e");
    const tiles = ["#8c7460", "#a48a6c", "#776957", "#947b68"];
    for (let row = 0; row < 5; row++) {
      for (let col = -1; col < 15; col++) {
        const bx = (col * w) / 13 + ((row % 2) * w) / 26;
        const by = h * 0.82 + row * h * 0.043;
        box(bx, by, w / 13 - 3, h * 0.043 - 3, tiles[(col + row * 3 + 16) % 4]);
        box(bx + 2, by + 2, w / 13 - 8, 2, "#c0a085");
        box(bx + 9, by + 8, 4, 2, "#64564e");
      }
    }
    // Solid wood pedestal and copper feet give the device a physical anchor.
    oval(x, y + r * 1.27, r * 1.14, r * 0.13, "#302e36", 12);
    box(x - r * 0.64, y + r * 0.86, r * 1.28, r * 0.39, "#694736");
    box(x - r * 0.7, y + r * 0.94, r * 1.4, 9, "#d19c60");
    box(x - r * 0.58, y + r * 1.22, r * 1.16, 7, "#b87d4d");
    ctx.save();
    ctx.globalAlpha = pending ? 0.28 : 1;
    ctx.shadowColor = "#68ead8";
    ctx.shadowBlur = pending ? 0 : 16;
    oval(x, y, r + 5, r + 5, "#513d38", 16);
    oval(x, y, r, r, "#b97e50", 10);
    oval(x, y, r - 4, r - 4, "#96f6df", 2);
    ctx.shadowBlur = 0;
    oval(x, y, r * 0.85, r * 0.28, "#3b7c91");
    oval(x, y, r * 0.32, r * 0.85, "#3b7c91");
    for (let i = 0; i < 24; i++) {
      const a = (i * Math.PI) / 12;
      line(
        x + Math.cos(a) * (r - 3),
        y + Math.sin(a) * (r - 3),
        x + Math.cos(a) * (r - 11),
        y + Math.sin(a) * (r - 11),
        "#c0e0d6",
      );
    }
    // Grip furniture is visible at rest; hover/focus adds directional gesture cues.
    const engaged = hovered || focused || !!drag;
    const gripColor = engaged ? "#fff0bd" : "#c99d6b";
    if (engaged) oval(x, y, r + 13, r + 13, gripColor, 2);
    for (const side of [-1, 1]) {
      const gx = x + side * (r + (drag ? 0 : 3));
      box(gx - 9, y - 19, 18, 38, "#493c38");
      box(gx - 6, y - 16, 12, 32, gripColor);
      for (let i = -1; i <= 1; i++) box(gx - 4, y + i * 7, 8, 2, "#79553f");
    }
    if (engaged)
      text(
        drag ? "↔" : tool === "rotate" ? "←  DRAG  →" : "↕  SCAN  ↔",
        x,
        y + r * 0.72,
        gripColor,
      );
    if (state) {
      const [scalar, ...coeff] = state.coordinates,
        negative = scalar < 0,
        color = negative ? "#ffca83" : "#8affdf";
      const b = state.bloch.length >= 3 ? state.bloch : coeff;
      const nx = x + (b[0] - 0.45 * b[2]) * r * 0.63,
        ny = y - (b[1] - 0.3 * b[2]) * r * 0.63;
      line(x, y, nx, ny, "#f0fff6", 4);
      ctx.fillStyle = color;
      ctx.beginPath();
      if (negative) {
        ctx.moveTo(nx, ny - 9);
        ctx.lineTo(nx + 9, ny);
        ctx.lineTo(nx, ny + 9);
        ctx.lineTo(nx - 9, ny);
        ctx.closePath();
      } else ctx.arc(nx, ny, 7, 0, Math.PI * 2);
      ctx.fill();
      const identity =
        coeff.every((c) => Math.abs(c) < 1e-6) &&
        Math.abs(Math.abs(scalar) - 1) < 1e-6;
      text(
        identity
          ? negative
            ? "−I · diamond"
            : "+I · circle"
          : negative
            ? "a₀ < 0 · diamond"
            : "a₀ ≥ 0 · circle",
        x,
        y + r + 45,
        color,
      );
      if (tool === "slice") {
        const s = state.slice,
          sy = y - s.coordinate * r * 0.62;
        ctx.fillStyle = "#ffc77d18";
        ctx.fillRect(x - r, sy - 7, r * 2, 14);
        line(x - r, sy, x + r, sy, "#ffd18f");
        for (const side of [-1, 1]) {
          box(x + side * r - 10, sy - 11, 20, 22, "#694936");
          box(x + side * r - 7, sy - 8, 14, 16, gripColor);
        }
        if (s.status === "sphere" && s.radius !== null) {
          oval(
            x,
            sy,
            Math.max(0, s.radius) * r * 0.72,
            Math.max(0, s.radius) * r * 0.2,
            "#ffe0ab",
            3,
          );
        } else if (s.status === "point") box(x - 4, sy - 4, 8, 8, "#ffe0ab");
        else {
          line(x - 5, sy - 5, x + 5, sy + 5, "#ffc18f");
          line(x + 5, sy - 5, x - 5, sy + 5, "#ffc18f");
        }
        text(
          s.status === "empty"
            ? "Empty slice"
            : s.status === "point"
              ? "Point slice"
              : `Slice radius ${s.radius === null ? "—" : fmt(s.radius)}`,
          x,
          y - r - 20,
          "#ffe0ab",
        );
      }
    } else text("AWAITING SIGNAL", x, y);
    ctx.restore();
    box(
      w - 35,
      18,
      12,
      12,
      pending ? "#efbb71" : state ? "#8debc0" : "#637785",
    );
    if (pending) {
      line(w - 33, 21, w - 25, 21, "#302721");
      line(w - 33, 26, w - 25, 26, "#302721");
    } else oval(w - 29, 24, 3, 3, "#163832");
  }
  function stop() {
    cancelAnimationFrame(frame);
    frame = 0;
  }
  function tick(now: number) {
    frame = 0;
    if (!visible()) return;
    draw(now);
    if (!motion.matches) frame = requestAnimationFrame(tick);
  }
  function refresh() {
    stop();
    if (visible()) {
      draw();
      if (!motion.matches) frame = requestAnimationFrame(tick);
    }
  }
  function resize() {
    const rect = root.getBoundingClientRect();
    w = rect.width;
    h = rect.height;
    const d = Math.min(devicePixelRatio || 1, 2);
    canvas.width = Math.round(w * d);
    canvas.height = Math.round(h * d);
    refresh();
  }
  function release() {
    const held = drag;
    drag = undefined;
    affordance();
    if (held && canvas.hasPointerCapture(held.id))
      canvas.releasePointerCapture(held.id);
  }
  function affordance() {
    root.dataset.grip = drag
      ? "dragging"
      : hovered || focused
        ? "hover"
        : "ready";
    canvas.style.cursor = drag
      ? "grabbing"
      : hovered || focused
        ? "grab"
        : "default";
    canvas.style.outline = focused ? "3px solid #ffe3a0" : "";
    refresh();
  }
  function hit(e: PointerEvent) {
    const rect = canvas.getBoundingClientRect(),
      { x, y, r } = layout();
    return (
      !!state &&
      Math.hypot(e.clientX - rect.left - x, e.clientY - rect.top - y) <=
        r + (e.pointerType === "touch" ? 28 : 14)
    );
  }
  function emit(value: number) {
    if (tool === "rotate") {
      angle = Math.max(0, Math.min(4 * Math.PI, value));
      actions.angle(angle);
    } else {
      slice = Math.max(-1.2, Math.min(1.2, value));
      actions.slice(slice);
    }
    describe();
  }
  canvas.addEventListener(
    "pointerdown",
    (e) => {
      if (e.button !== 0 || drag || !state) return;
      if (!hit(e)) return;
      canvas.focus({ preventScroll: true });
      canvas.setPointerCapture(e.pointerId);
      canvas.style.cursor = "grabbing";
      drag = {
        id: e.pointerId,
        x: e.clientX,
        y: e.clientY,
        value: tool === "rotate" ? angle : slice,
      };
      affordance();
      actions.sound();
      e.preventDefault();
    },
    { signal: events.signal },
  );
  canvas.addEventListener(
    "pointermove",
    (e) => {
      hovered = hit(e);
      if (!drag || drag.id !== e.pointerId) {
        affordance();
        return;
      }
      const dx = e.clientX - drag.x,
        dy = drag.y - e.clientY;
      emit(
        drag.value +
          (tool === "rotate"
            ? dx * 4 * Math.PI
            : (Math.abs(dx) >= Math.abs(dy) ? dx : dy) * 2.4) /
            (2 * layout().r),
      );
    },
    { signal: events.signal },
  );
  for (const type of ["pointerup", "pointercancel", "lostpointercapture"])
    canvas.addEventListener(
      type,
      (e) => {
        if ((e as PointerEvent).pointerId === drag?.id) release();
      },
      { signal: events.signal },
    );
  canvas.addEventListener(
    "keydown",
    (e) => {
      if (
        document.activeElement !== canvas ||
        !state ||
        e.altKey ||
        e.ctrlKey ||
        e.metaKey
      )
        return;
      let value: number;
      if (tool === "rotate") {
        switch (e.key) {
          case "ArrowLeft":
            value = angle - Math.PI / 24;
            break;
          case "ArrowRight":
            value = angle + Math.PI / 24;
            break;
          case "Home":
            value = 0;
            break;
          case "End":
            value = 4 * Math.PI;
            break;
          default:
            return;
        }
      } else {
        switch (e.key) {
          case "ArrowUp":
            value = slice + 0.05;
            break;
          case "ArrowDown":
            value = slice - 0.05;
            break;
          default:
            return;
        }
      }
      e.preventDefault();
      actions.sound();
      emit(value);
    },
    { signal: events.signal },
  );
  canvas.addEventListener(
    "pointerleave",
    () => {
      hovered = false;
      affordance();
    },
    { signal: events.signal },
  );
  canvas.addEventListener(
    "focus",
    () => {
      focused = true;
      affordance();
    },
    { signal: events.signal },
  );
  canvas.addEventListener(
    "blur",
    () => {
      focused = false;
      release();
    },
    { signal: events.signal },
  );
  window.addEventListener("blur", release, { signal: events.signal });
  document.addEventListener(
    "visibilitychange",
    () => {
      if (document.hidden) release();
      refresh();
    },
    { signal: events.signal },
  );
  motion.addEventListener("change", refresh, { signal: events.signal });
  const observer = new ResizeObserver(resize);
  observer.observe(root);
  const intersection = new IntersectionObserver((entries) => {
    inside = entries.some((e) => e.isIntersecting);
    if (!inside) release();
    refresh();
  });
  intersection.observe(root);
  describe();
  resize();
  return {
    update(next) {
      if (dead) return;
      state = {
        t: next.t,
        axis: [...next.axis],
        coordinates: [...next.coordinates],
        bloch: [...next.bloch],
        slice: { ...next.slice },
      };
      if (!drag) {
        angle = next.t;
        slice = next.slice.coordinate;
      }
      pending = false;
      describe();
      refresh();
    },
    setPending(value) {
      if (dead) return;
      pending = value;
      describe();
      refresh();
    },
    setTool(value) {
      if (dead) return;
      release();
      tool = value;
      describe();
      refresh();
    },
    dispose() {
      if (dead) return;
      dead = true;
      stop();
      release();
      events.abort();
      observer.disconnect();
      intersection.disconnect();
      root.remove();
    },
  };
}
