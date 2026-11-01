/** Original hand-built pixel garden. Input belongs to this canvas, never the page. */
export function mountWorld(
  host: HTMLElement,
  onEnter: () => void,
  onSound: () => void,
): { dispose(): void } {
  const root = document.createElement("section");
  root.className = "observatory-world";
  root.innerHTML = `<div class="world-hud"><div><span class="world-eyebrow">FIELD SITE 01 · FREE EXPLORATION</span><h2>The Spin Observatory</h2></div><span class="world-badge">SU(2) · garden entrance</span></div>
    <p class="world-instructions">Wander the garden at your own pace. Focus the world, then use WASD or arrow keys to walk. Near the observatory, press Enter, E or Space to enter. You can also tap the observatory or open the lab directly.</p>
    <canvas class="world-canvas" width="640" height="360" tabindex="0" aria-label="Observatory garden. Use WASD or arrow keys to walk; Enter, E or Space near the observatory opens the SU(2) lab.">Explore the observatory using the Open SU(2) lab button below.</canvas>
    <div class="world-footer"><div class="world-dpad" role="group" aria-label="Walk around the garden"><button type="button" data-direction="up" aria-label="Walk north">↑</button><button type="button" data-direction="left" aria-label="Walk west">←</button><button type="button" data-direction="down" aria-label="Walk south">↓</button><button type="button" data-direction="right" aria-label="Walk east">→</button></div><p class="world-status" role="status" aria-live="polite"></p><button type="button" class="world-enter">Open SU(2) lab <span aria-hidden="true">↗</span></button></div>`;
  host.append(root);
  const canvas = root.querySelector("canvas")!;
  const ctx = canvas.getContext("2d");
  const status = root.querySelector<HTMLElement>(".world-status")!;
  const motion = matchMedia("(prefers-reduced-motion: reduce)");
  const events = new AbortController();
  const keys = new Set<string>();
  const held = new Map<number, string>();
  let x = 312,
    y = 282,
    frame = 0,
    last = 0,
    dead = false,
    near = false,
    targeted = false;
  const directions: Record<string, readonly [number, number]> = {
    up: [0, -1],
    down: [0, 1],
    left: [-1, 0],
    right: [1, 0],
  };
  const mapping: Record<string, string> = {
    w: "up",
    arrowup: "up",
    s: "down",
    arrowdown: "down",
    a: "left",
    arrowleft: "left",
    d: "right",
    arrowright: "right",
  };
  const listen = (target: EventTarget, name: string, handler: EventListener) =>
    target.addEventListener(name, handler, { signal: events.signal });
  const enter = () => {
    if (!dead) {
      onSound();
      onEnter();
    }
  };
  const updateStatus = () => {
    near = Math.hypot(x - 466, y - 196) < 63;
    const message = targeted
      ? "Spin Observatory · click or tap the highlighted doorway to enter."
      : near
        ? "Observatory reached · Enter / E / Space to open the lab."
        : "Follow the lantern path east. Exploration is optional.";
    if (status.textContent !== message) status.textContent = message;
  };
  const blocked = (px: number, py: number) =>
    px < 24 ||
    px > 616 ||
    py < 40 ||
    py > 332 ||
    (px > 80 && px < 239 && py > 115 && py < 223) ||
    (px > 406 && px < 530 && py > 81 && py < 176);
  const move = (dx: number, dy: number, distance: number) => {
    if (!blocked(x + dx * distance, y)) x += dx * distance;
    if (!blocked(x, y + dy * distance)) y += dy * distance;
    updateStatus();
  };
  function draw(_time: number) {
    if (!ctx) return;
    ctx.imageSmoothingEnabled = false;
    const box = (a: number, b: number, w: number, h: number, color: string) => {
      ctx.fillStyle = color;
      ctx.fillRect(Math.round(a), Math.round(b), w, h);
    };
    const text = (
      s: string,
      a: number,
      b: number,
      color = "#fff3cf",
      size = 12,
    ) => {
      ctx.fillStyle = color;
      ctx.font = `${size}px monospace`;
      ctx.fillText(s, a, b);
    };
    box(0, 0, 640, 360, "#244d48");
    box(12, 28, 616, 316, "#527d58");
    box(24, 40, 592, 292, "#789856");
    for (let i = 0; i < 150; i++) {
      const a = 28 + ((i * 97) % 580),
        b = 45 + ((i * 61) % 278);
      box(a, b, 3, 2, i % 3 ? "#65864e" : "#a4b971");
    }
    // Fixed tile motifs: repeatable clusters, never random decoration.
    for (let row = 0; row < 18; row++)
      for (let col = 0; col < 37; col++) {
        const a = 25 + col * 16,
          b = 41 + row * 16;
        const tone = (row * 7 + col * 11) % 5;
        box(a + 3, b + 10, 5, 2, tone < 2 ? "#648447" : "#86a45f");
        box(a + 5, b + 7, 2, 4, "#a1b66d");
        if (tone === 0) box(a + 11, b + 3, 2, 2, "#c0ca81");
      }
    // Layered hedge tiles and warm stone retaining border.
    for (let i = 0; i < 38; i++) {
      const a = 16 + i * 16;
      box(a, 31, 16, 9, "#354e36");
      box(a, 25, 15, 9, "#557548");
      box(a + 2, 25, 10, 2, "#94a968");
      box(a, 338, 15, 6, "#7e7254");
      box(a, 337, 14, 2, "#c5b58a");
    }
    // Stepped sandstone paths and the central round plaza.
    box(286, 185, 52, 147, "#b89d69");
    box(298, 196, 28, 136, "#dbc28c");
    box(312, 185, 182, 38, "#b89d69");
    box(326, 194, 168, 20, "#dbc28c");
    box(270, 245, 84, 43, "#b89d69");
    box(278, 252, 68, 29, "#e4cc96");
    for (let i = 0; i < 9; i++) box(300, 222 + i * 12, 24, 1, "#bba16d");
    // Cut stone joints and a small timber crossing on the garden route.
    for (let i = 0; i < 10; i++) {
      box(331 + i * 16, 195, 1, 18, "#b49a6c");
      box(333 + i * 16, 196, 12, 2, "#f0d9a2");
    }
    box(284, 227, 56, 25, "#584839");
    for (let i = 0; i < 7; i++) {
      box(287 + i * 7, 228, 6, 21, "#b78351");
      box(288 + i * 7, 229, 4, 2, "#e1b078");
      box(289 + i * 7, 245, 2, 2, "#6c513a");
    }
    for (const b of [224, 250]) box(281, b, 62, 3, "#e1b078");
    for (const [a, b] of [
      [312, 298],
      [378, 204],
      [418, 204],
    ]) {
      box(a - 4, b - 3, 6, 2, "#876b47");
      box(a, b - 1, 4, 2, "#876b47");
      box(a - 4, b + 1, 6, 2, "#876b47");
    }
    // Pond: deliberately stepped banks, reeds and lily pads.
    box(91, 133, 136, 73, "#3d6249");
    box(103, 120, 111, 100, "#3d6249");
    box(101, 139, 117, 58, "#346e7a");
    box(113, 130, 91, 79, "#346e7a");
    box(112, 145, 95, 37, "#408c95");
    const shimmer = 0; // Static reflections also honor reduced motion.
    for (let i = 0; i < 6; i++)
      box(119 + i * 13, 147 + (i % 3) * 16 + shimmer, 10, 2, "#81bab1");
    // Broken dome/sky reflections beneath the far bank.
    for (let i = 0; i < 5; i++) {
      box(
        151 - i * 3,
        151 + i * 5,
        28 + i * 6,
        2,
        i % 2 ? "#5a9c9d" : "#9dc5b6",
      );
    }
    box(132, 185, 17, 6, "#8daf63");
    box(138, 181, 6, 4, "#f5cfaa");
    box(190, 137, 12, 6, "#8daf63");
    for (const [a, b] of [
      [94, 151],
      [213, 192],
      [108, 212],
    ]) {
      box(a, b - 13, 2, 18, "#cfb66a");
      box(a + 5, b - 8, 2, 14, "#b4c17c");
    }
    // Evergreen geometry, orchard flowers and boundary fence.
    for (const [a, b] of [
      [48, 72],
      [260, 84],
      [570, 91],
      [576, 278],
      [61, 281],
      [375, 307],
    ]) {
      box(a - 4, b + 5, 9, 22, "#634e3b");
      box(a - 20, b - 3, 42, 19, "#315d48");
      box(a - 15, b - 17, 32, 21, "#3c7050");
      box(a - 9, b - 28, 20, 19, "#4d8058");
      box(a - 4, b - 24, 5, 9, "#86a768");
    }
    for (let i = 0; i < 22; i++) {
      const a = 38 + i * 26;
      box(a, 322, 3, 18, "#b5ad7b");
      box(a, 327, 23, 3, "#c8bd89");
    }
    for (let i = 0; i < 18; i++) {
      const a = 365 + ((i * 13) % 190),
        b = 240 + ((i * 19) % 60);
      box(a, b, 2, 7, "#42664b");
      box(a - 2, b - 2, 6, 4, i % 2 ? "#edcc8a" : "#ddacb0");
    }
    // Ferns, berry shrubs and stepping-stone shadows.
    for (const [a, b] of [
      [46, 197],
      [251, 152],
      [552, 233],
      [228, 292],
      [532, 303],
    ]) {
      box(a - 15, b + 6, 34, 7, "#4d703f");
      box(a - 17, b - 3, 34, 11, "#345d42");
      box(a - 12, b - 11, 25, 14, "#4e7d4e");
      box(a - 9, b - 10, 17, 3, "#91ae65");
      for (let i = 0; i < 3; i++)
        box(a - 8 + i * 8, b - 3 + (i % 2) * 3, 3, 3, "#dda079");
    }
    // Original octagonal dome, telescope slit, masonry and mathematical station.
    box(400, 164, 144, 14, "#476149");
    box(413, 112, 111, 56, "#e0d5ae");
    box(408, 105, 121, 12, "#b9aa83");
    box(413, 83, 111, 24, "#527d89");
    box(424, 69, 89, 18, "#6599a0");
    box(438, 59, 60, 15, "#8ab5b1");
    box(450, 54, 35, 8, "#adcac0");
    box(468, 57, 10, 48, "#254c5b");
    box(478, 71, 39, 8, "#d1ba7f");
    box(507, 66, 12, 15, "#e5d39d");
    // Brick courses and warm window mullions.
    for (let row = 0; row < 4; row++)
      for (let col = 0; col < 7; col++) {
        box(415 + col * 15 + (row % 2) * 5, 117 + row * 12, 9, 1, "#c3b68e");
      }
    box(455, 133, 28, 35, "#36564f");
    box(460, 140, 18, 28, "#243f43");
    box(474, 153, 3, 3, "#f4d487");
    box(422, 126, 18, 17, "#396c79");
    box(497, 126, 18, 17, "#396c79");
    for (const a of [422, 497]) {
      box(a + 2, 128, 14, 2, "#a4c6bb");
      box(a + 8, 126, 2, 17, "#d2ba87");
      box(a, 134, 18, 2, "#d2ba87");
    }
    const active = near || targeted;
    if (active) {
      ctx.strokeStyle = "#ffe8a1";
      ctx.lineWidth = 2;
      ctx.strokeRect(452, 130, 34, 40);
      box(446, 170, 46, 3, "#ffe8a1");
    }
    box(444, 180, 50, 24, "#304b48");
    box(448, 180, 42, 19, active ? "#f8dc92" : "#d6c99e");
    text("SU(2)", 451, 194, "#29483f", 11);
    for (const a of [355, 540]) {
      box(a, 188, 4, 30, "#465947");
      box(a - 4, 181, 12, 12, "#f2d38a");
      box(a - 1, 184, 6, 6, "#fff1b1");
    }
    text("THE SPIN OBSERVATORY", 359, 39, "#ecdfb3", 12);
    box(51, 244, 122, 25, "#385b4a");
    text("REFLECTION POND", 58, 260, "#f0e2b7", 11);
    // Explorer with boots, field jacket and copper hair.
    box(x - 8, y + 5, 18, 5, "#456347");
    box(x - 6, y - 1, 5, 9, "#293f45");
    box(x + 2, y - 1, 5, 9, "#293f45");
    box(x - 7, y - 13, 15, 15, "#e4b668");
    box(x - 10, y - 10, 4, 11, "#e9c99b");
    box(x + 8, y - 10, 4, 11, "#e9c99b");
    box(x - 6, y - 24, 13, 12, "#edcca1");
    box(x - 8, y - 27, 16, 7, "#784d3c");
    box(x - 8, y - 22, 4, 6, "#784d3c");
    box(x + 3, y - 19, 2, 2, "#33464a");
    box(x - 3, y - 10, 6, 8, "#698d88");
    // Field coat lapels, pocket notebook, glasses and sturdy boot highlights.
    box(x - 7, y - 12, 4, 14, "#f4e4ba");
    box(x + 4, y - 12, 4, 14, "#f4e4ba");
    box(x - 6, y - 6, 3, 4, "#729b8e");
    box(x + 8, y - 5, 5, 7, "#785c49");
    box(x + 9, y - 5, 3, 1, "#f4d29a");
    box(x - 5, y - 20, 5, 4, "#364f50");
    box(x + 2, y - 20, 5, 4, "#364f50");
    box(x, y - 19, 2, 1, "#364f50");
    box(x - 4, y - 19, 3, 1, "#d6e7d3");
    box(x + 3, y - 19, 3, 1, "#d6e7d3");
    box(x + 1, y - 14, 3, 1, "#ad6d51");
    box(x - 6, y + 6, 5, 2, "#ac8050");
    box(x + 2, y + 6, 5, 2, "#ac8050");
    if (targeted) {
      box(373, 211, 205, 22, "#283f38");
      text("CLICK / TAP TO ENTER", 386, 226, "#ffe8a1", 12);
    }
    if (near) {
      box(x - 31, y - 49, 96, 16, "#233f43");
      text("ENTER / E", x - 26, y - 37, "#ffe3a2", 11);
    }
    text("GARDEN / 01", 25, 20, "#e8dbb5", 11);
  }
  const tick = (now: number) => {
    if (dead) return;
    if (root.closest("[hidden]") || document.hidden) {
      keys.clear();
      held.clear();
      last = 0;
      frame = requestAnimationFrame(tick);
      return;
    }
    const dt = last ? Math.min((now - last) / 1000, 0.04) : 0;
    last = now;
    let dx = 0,
      dy = 0;
    for (const direction of new Set([...keys, ...held.values()])) {
      const v = directions[direction];
      if (v) {
        dx += v[0];
        dy += v[1];
      }
    }
    const length = Math.hypot(dx, dy);
    if (length) move(dx / length, dy / length, 95 * dt);
    draw(now);
    frame = requestAnimationFrame(tick);
  };
  listen(canvas, "keydown", (event) => {
    const e = event as KeyboardEvent;
    if (document.activeElement !== canvas || e.ctrlKey || e.metaKey || e.altKey)
      return;
    const direction = mapping[e.key.toLowerCase()];
    if (direction) {
      e.preventDefault();
      keys.add(direction);
    } else if (
      (e.key === "Enter" || e.key.toLowerCase() === "e" || e.key === " ") &&
      near
    ) {
      e.preventDefault();
      if (!e.repeat) enter();
    }
  });
  listen(canvas, "keyup", (event) => {
    const e = event as KeyboardEvent;
    const direction = mapping[e.key.toLowerCase()];
    if (direction) {
      e.preventDefault();
      keys.delete(direction);
    }
  });
  listen(canvas, "blur", () => keys.clear());
  listen(window, "blur", () => {
    keys.clear();
    held.clear();
  });
  listen(document, "visibilitychange", () => {
    keys.clear();
    held.clear();
    last = 0;
  });
  const targetAt = (event: Event) => {
    const e = event as PointerEvent,
      bounds = canvas.getBoundingClientRect();
    const px = ((e.clientX - bounds.left) * 640) / bounds.width;
    const py = ((e.clientY - bounds.top) * 360) / bounds.height;
    targeted = px >= 400 && px <= 544 && py >= 54 && py <= 208;
    canvas.style.cursor = targeted ? "pointer" : "default";
    updateStatus();
  };
  listen(canvas, "pointermove", targetAt);
  listen(canvas, "pointerdown", targetAt);
  listen(canvas, "pointerleave", () => {
    targeted = false;
    updateStatus();
  });
  listen(canvas, "pointercancel", () => {
    targeted = false;
    updateStatus();
  });
  listen(canvas, "click", (event) => {
    const e = event as MouseEvent,
      bounds = canvas.getBoundingClientRect();
    canvas.focus({ preventScroll: true });
    const px = ((e.clientX - bounds.left) * 640) / bounds.width,
      py = ((e.clientY - bounds.top) * 360) / bounds.height;
    if (px >= 400 && px <= 544 && py >= 54 && py <= 208) enter();
  });
  listen(root.querySelector(".world-enter")!, "click", enter);
  for (const button of root.querySelectorAll<HTMLButtonElement>(
    "[data-direction]",
  )) {
    const direction = button.dataset.direction!;
    listen(button, "pointerdown", (event) => {
      const e = event as PointerEvent;
      if (e.button !== 0) return;
      button.setPointerCapture(e.pointerId);
      held.set(e.pointerId, direction);
      onSound();
    });
    const release = (event: Event) =>
      held.delete((event as PointerEvent).pointerId);
    listen(button, "pointerup", release);
    listen(button, "pointercancel", release);
    listen(button, "lostpointercapture", release);
    listen(button, "blur", () => held.clear());
    listen(button, "click", (event) => {
      if ((event as MouseEvent).detail === 0) {
        const v = directions[direction]!;
        move(v[0], v[1], 16);
        onSound();
        draw(0);
      }
    });
  }
  listen(motion, "change", () => draw(0));
  updateStatus();
  draw(0);
  if (ctx) frame = requestAnimationFrame(tick);
  return {
    dispose() {
      if (dead) return;
      dead = true;
      cancelAnimationFrame(frame);
      events.abort();
      keys.clear();
      held.clear();
      root.remove();
    },
  };
}
