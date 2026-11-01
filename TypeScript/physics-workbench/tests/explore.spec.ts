import { test, expect } from "@playwright/test";
async function lab(page: import("@playwright/test").Page) {
  await page.goto("/explore.html");
  await page
    .locator(".explore-header")
    .getByRole("button", { name: "Open SU(2) lab", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Scan", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
}
test("standalone Explore loads without notebook UI, code or storage access", async ({
  page,
}) => {
  const errors: string[] = [];
  page.on("pageerror", (e) => errors.push(e.message));
  await page.addInitScript(() => {
    const original = Storage.prototype.getItem;
    Storage.prototype.getItem = function (key) {
      if (key === "mathphysics.workbench.v1")
        throw Error("Notebook storage read");
      return original.call(this, key);
    };
  });
  await page.goto("/explore.html");
  await expect(page.locator(".world-canvas")).toBeVisible();
  await expect(
    page.locator("#tabs,.sidebar,.inspector,.notebook,[data-cell-id]"),
  ).toHaveCount(0);
  await expect(page.locator("body")).not.toContainText(
    /\bMathematica\b|\bnotebook\b|Run cell|Pin to notebook/i,
  );
  expect(errors).toEqual([]);
  await page.reload();
  await expect(page.locator(".world-canvas")).toBeVisible();
});
test("notebook remains separate and its data survives opening the game", async ({
  page,
}) => {
  await page.goto("/");
  await page
    .locator('[data-cell-id="first-trace"] textarea')
    .fill("preserved_mu, preserved_nu");
  const before = await page.evaluate(() =>
    localStorage.getItem("mathphysics.workbench.v1"),
  );
  await page.goto("/explore.html");
  await expect(page.locator(".world-canvas")).toBeVisible();
  expect(
    await page.evaluate(() => localStorage.getItem("mathphysics.workbench.v1")),
  ).toBe(before);
  await page.goto("/");
  await expect(
    page.locator('[data-cell-id="first-trace"] textarea'),
  ).toHaveValue("preserved_mu, preserved_nu");
  await expect(
    page.getByRole("tab", { name: "✦ Explore", exact: true }),
  ).toHaveCount(0);
});
test("Rust SU2 has independent matrix checkpoints, rejection and slice boundaries", async ({
  request,
}) => {
  for (const [t, real, imag] of [
    [0, 1, 0],
    [Math.PI, 0, -1],
    [2 * Math.PI, -1, 0],
    [4 * Math.PI, 1, 0],
  ]) {
    const r = await request.post("/api/su2", {
      data: { axis: [0, 0, 1], t, slice: 0 },
    });
    expect(r.ok()).toBe(true);
    const s = await r.json();
    expect(s.matrix[0][0][0]).toBeCloseTo(real, 12);
    expect(s.matrix[0][0][1]).toBeCloseTo(imag, 12);
    expect(s.checks.unitarity_error).toBeLessThan(1e-12);
  }
  for (const slice of [-1, 1, 1.1]) {
    const r = await request.post("/api/su2", {
      data: { axis: [1, 0, 0], t: 1, slice },
    });
    expect((await r.json()).slice.status).toBe(
      slice === 1.1 ? "empty" : "point",
    );
  }
  const invalid = await request.post("/api/su2", {
    data: { axis: [2, 0, 0], t: 0, slice: 0 },
  });
  expect(invalid.ok()).toBe(false);
});
test("laboratory checkpoints expose sign and current Rust evidence", async ({
  page,
}) => {
  await lab(page);
  await page
    .locator(".explore-checkpoints")
    .getByRole("button", { name: "2π", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await expect(page.locator(".explore-matrix")).toContainText("-1.00000");
  await expect(page.locator("#explore-angle")).toContainText("2.000π");
  await page
    .locator(".explore-checkpoints")
    .getByRole("button", { name: "4π", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await expect(page.locator(".explore-matrix")).toContainText("1.00000");
  await expect(page.locator(".explore-matrix")).not.toContainText("-1.00000");
  await expect(page.locator(".explore-plots canvas").first()).toBeVisible();
});
test("slice empty state and explicit custom axis validation", async ({
  page,
}) => {
  await lab(page);
  await page
    .getByRole("slider", { name: "Slice coordinate", exact: true })
    .fill("1.1");
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await expect(page.locator(".explore-panels")).toContainText(
    /empty|no (real )?slice/i,
  );
  await page
    .getByRole("combobox", { name: "Rotation axis", exact: true })
    .selectOption("custom");
  await page.getByRole("spinbutton", { name: "Axis x", exact: true }).fill("0");
  await page.getByRole("spinbutton", { name: "Axis y", exact: true }).fill("0");
  await page.getByRole("spinbutton", { name: "Axis z", exact: true }).fill("0");
  await page.getByRole("button", { name: "Apply axis", exact: true }).click();
  await expect(
    page.locator(".explore-shell p.progress-note[role=status]"),
  ).toContainText("unit length");
});
test("journal reopens exact inputs after reload and exports a discovery", async ({
  page,
}) => {
  await lab(page);
  await page
    .locator(".explore-checkpoints")
    .getByRole("button", { name: "π", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await page
    .getByRole("textbox", { name: "Experiment title", exact: true })
    .fill("Spinor sign investigation");
  await page
    .getByRole("textbox", { name: "Experiment observation", exact: true })
    .fill("Compare after one and two turns.");
  await page
    .getByRole("button", { name: "Record experiment", exact: true })
    .click();
  await expect(page.locator(".explore-journal-entry")).toHaveCount(1);
  await page.reload();
  await page
    .getByRole("button", { name: "Field journal", exact: true })
    .click();
  await expect(page.locator(".explore-journal-entry")).toContainText(
    "Spinor sign investigation",
  );
  await page
    .locator(".explore-journal-entry")
    .getByRole("button", { name: "Reopen", exact: true })
    .click();
  await expect(page.locator("#explore-angle")).toContainText("1.000π");
  await page
    .getByRole("button", { name: "Field journal", exact: true })
    .click();
  const exported = page.waitForEvent("download");
  await page
    .locator(".explore-journal-entry")
    .getByRole("button", { name: "Export discovery", exact: true })
    .click();
  expect((await exported).suggestedFilename()).toBe("discovery.json");
  await expect(
    page.getByRole("button", { name: "Pin to notebook", exact: true }),
  ).toHaveCount(0);
});
test("playback pauses reproducibly and reset selects identity", async ({
  page,
}) => {
  await lab(page);
  await page
    .locator(".explore-controls")
    .getByRole("button", { name: "Play", exact: true })
    .click();
  await page.waitForTimeout(400);
  await page
    .locator(".explore-controls")
    .getByRole("button", { name: "Pause", exact: true })
    .click();
  const before = await page.locator("#explore-angle").textContent();
  await page.waitForTimeout(250);
  expect(await page.locator("#explore-angle").textContent()).toBe(before);
  await page
    .getByRole("button", { name: "Reset to identity", exact: true })
    .click();
  await expect(page.locator("#explore-angle")).toContainText("0.000π");
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
});
test("unknown exploration data survives safe recovery", async ({ page }) => {
  const raw = JSON.stringify({ version: 999, privateNote: "preserve-me" });
  await page.addInitScript(
    (value) => localStorage.setItem("mathphysics.explore.v1", value),
    raw,
  );
  await page.goto("/explore.html");
  await expect(page.locator(".explore-recovery")).toContainText("preserved");
  expect(
    await page.evaluate(() => localStorage.getItem("mathphysics.explore.v1")),
  ).toBe(raw);
});
test("mobile Explore remains within viewport and direct entry works", async ({
  page,
}) => {
  await page.setViewportSize({ width: 460, height: 900 });
  await page.goto("/explore.html");
  await page
    .locator(".explore-header")
    .getByRole("button", { name: "Open SU(2) lab", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
});

test("slow Rust responses do not starve playback or label stale inputs checked", async ({
  page,
}) => {
  await lab(page);
  let completed = 0;
  await page.route("**/api/su2", async (route) => {
    const response = await route.fetch();
    await new Promise((resolve) => setTimeout(resolve, 220));
    try {
      await route.fulfill({ response });
      completed++;
    } catch {
      /* A manual pause may cancel the last sample. */
    }
  });
  await page
    .locator(".explore-controls")
    .getByRole("button", { name: "Play", exact: true })
    .click();
  await page.waitForTimeout(1050);
  expect(completed).toBeGreaterThanOrEqual(2);
  await page
    .locator(".explore-controls")
    .getByRole("button", { name: "Pause", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  const t = Number(
    await page.locator(".explore-evidence").getAttribute("data-t"),
  );
  expect(t).toBeGreaterThan(0.7);
  await expect(page.locator("#explore-angle")).toContainText(t.toFixed(5));
  await page.reload();
  await expect(page.locator("#explore-angle")).toContainText(t.toFixed(5));
});
test("failed journal persistence reports unsaved data and still permits export", async ({
  page,
}) => {
  await page.addInitScript(() => {
    const original = Storage.prototype.setItem;
    Storage.prototype.setItem = function (key, value) {
      if (key === "mathphysics.explore.v1")
        throw new DOMException("Quota exhausted", "QuotaExceededError");
      return original.call(this, key, value);
    };
  });
  await lab(page);
  await page
    .getByRole("button", { name: "Record experiment", exact: true })
    .click();
  await expect(
    page.locator(".explore-shell p.progress-note[role=status]"),
  ).toContainText("in memory only");
  await expect(page.locator(".explore-journal-entry")).toHaveCount(1);
  await page.getByRole("button", { name: "Close scan", exact: true }).click();
  await page
    .getByRole("button", { name: "Field journal", exact: true })
    .click();
  const event = page.waitForEvent("download");
  await page
    .getByRole("button", { name: "Export exploration", exact: true })
    .click();
  expect((await event).suggestedFilename()).toBe("symmetry-exploration.json");
});

test("Explore export shortcut writes exploration state rather than the notebook", async ({
  page,
}) => {
  await lab(page);
  const event = page.waitForEvent("download");
  await page.keyboard.press("Control+s");
  expect((await event).suggestedFilename()).toBe("symmetry-exploration.json");
});

test("scrubbing during playback selects the requested angle before pausing", async ({
  page,
}) => {
  await lab(page);
  await page
    .locator(".explore-controls")
    .getByRole("button", { name: "Play", exact: true })
    .click();
  await page.waitForTimeout(150);
  await page
    .getByRole("slider", { name: "Rotation angle t", exact: true })
    .fill("2");
  await expect(
    page
      .locator(".explore-controls")
      .getByRole("button", { name: "Play", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await expect(page.locator("#explore-angle")).toContainText("t = 2.00000");
});

async function scene(page: import("@playwright/test").Page) {
  await page.goto("/explore.html");
  await page
    .locator(".explore-header")
    .getByRole("button", { name: "Open SU(2) lab", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
}
test("scene dominates the default view and scanner/journal are opt-in", async ({
  page,
}) => {
  await scene(page);
  await expect(page.locator(".explore-scene-host canvas")).toBeVisible();
  await expect(page.locator(".explore-matrix")).not.toBeVisible();
  await expect(page.locator(".explore-plots").first()).not.toBeVisible();
  await expect(page.locator(".explore-journal")).not.toBeVisible();
  const b = await page.locator(".explore-scene-host").boundingBox();
  expect(b!.height).toBeGreaterThan(500);
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Scan", exact: true })
    .click();
  await expect(page.locator(".explore-matrix")).toBeVisible();
  await page.keyboard.press("Escape");
  await expect(page.locator(".explore-matrix")).not.toBeVisible();
  await expect(
    page
      .locator(".scene-dock")
      .getByRole("button", { name: "Scan", exact: true }),
  ).toBeFocused();
});
test("direct scene rotation and slice controls change real Rust state without opening Scan", async ({
  page,
}) => {
  await scene(page);
  const canvas = page.locator(".explore-scene-host canvas");
  await canvas.focus();
  await page.keyboard.press("End");
  await expect(page.locator(".scene-computation")).toContainText("4.00π");
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Slice", exact: true })
    .click();
  await canvas.focus();
  await page.keyboard.press("ArrowUp");
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  expect(await page.locator(".explore-scene-host").textContent()).toMatch(
    /Slice/,
  );
  await expect(page.locator("dialog[open]")).toHaveCount(0);
  await page
    .getByRole("button", { name: "Change rotation axis", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await expect(
    page.getByRole("button", { name: "Change rotation axis", exact: true }),
  ).toHaveText("Axis X ↻");
});
test("keep a discovery directly from the scene and open journal on demand", async ({
  page,
}) => {
  await scene(page);
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Keep discovery", exact: true })
    .click();
  await expect(page.locator(".explore-journal-entry")).toHaveCount(1);
  await expect(page.locator(".explore-journal")).not.toBeVisible();
  await page
    .getByRole("button", { name: "Field journal", exact: true })
    .click();
  await expect(page.locator(".explore-journal-entry")).toBeVisible();
  await page.keyboard.press("Escape");
  await expect(page.locator(".explore-journal")).not.toBeVisible();
});

test("optional Tips experiment uses the real one-turn and two-turn states", async ({
  page,
}) => {
  await scene(page);
  await page.getByRole("button", { name: "Tips", exact: true }).click();
  await expect(
    page.getByRole("dialog", { name: "Observatory tips", exact: true }),
  ).toBeVisible();
  await page
    .getByRole("button", { name: "One turn · 2π", exact: true })
    .click();
  await expect(page.locator(".scene-computation")).toContainText("2.00π");
  await expect(page.locator(".explore-matrix")).toContainText("-1.00000");
  await expect(page.locator("dialog[open]")).toHaveCount(0);
  await page.getByRole("button", { name: "Tips", exact: true }).click();
  await page
    .getByRole("button", { name: "Two turns · 4π", exact: true })
    .click();
  await expect(page.locator(".scene-computation")).toContainText("4.00π");
  await expect(page.locator(".explore-matrix")).not.toContainText("-1.00000");
});
test("whole apparatus has responsive grab cues and changes actual state", async ({
  page,
}) => {
  await scene(page);
  const canvas = page.locator(".explore-scene-host canvas");
  const b = (await canvas.boundingBox())!;
  await page.mouse.move(b.x + b.width * 0.5, b.y + b.height * 0.51);
  await expect(page.locator(".observatory-scene")).toHaveAttribute(
    "data-grip",
    "hover",
  );
  await page.mouse.down();
  await expect(page.locator(".observatory-scene")).toHaveAttribute(
    "data-grip",
    "dragging",
  );
  await page.mouse.move(b.x + b.width * 0.5 + 40, b.y + b.height * 0.51, {
    steps: 5,
  });
  await page.mouse.up();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  expect(
    Number(await page.locator(".explore-evidence").getAttribute("data-t")),
  ).toBeGreaterThan(0.1);
});
test("keeping a moving point waits for its exact result", async ({ page }) => {
  await scene(page);
  await page.route("**/api/su2", async (route) => {
    const response = await route.fetch();
    await new Promise((r) => setTimeout(r, 200));
    try {
      await route.fulfill({ response });
    } catch {}
  });
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Play", exact: true })
    .click();
  await page.waitForTimeout(150);
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Keep discovery", exact: true })
    .click();
  await expect(page.locator(".explore-journal-entry")).toHaveCount(1);
  const saved = await page.evaluate(() =>
    JSON.parse(localStorage.getItem("mathphysics.explore.v1")!),
  );
  expect(saved.journal[0].input.t).toBeGreaterThan(0);
  expect(saved.journal[0].input.t).toBe(
    Number(await page.locator(".explore-evidence").getAttribute("data-t")),
  );
});
test("failed capture never fabricates a kept discovery", async ({ page }) => {
  await scene(page);
  await page.route("**/api/su2", (route) =>
    route.fulfill({
      status: 503,
      contentType: "application/json",
      body: JSON.stringify({ error: { message: "Signal unavailable" } }),
    }),
  );
  await page
    .getByRole("button", { name: "Change rotation axis", exact: true })
    .click();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "error",
  );
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Keep discovery", exact: true })
    .click();
  await expect(
    page.locator(".explore-shell p.progress-note[role=status]"),
  ).toContainText("Could not keep");
  await expect(page.locator(".explore-journal-entry")).toHaveCount(0);
});
test("contextual shortcuts keep pictograms and only work at the instrument", async ({
  page,
}) => {
  await scene(page);
  const c = page.locator(".explore-scene-host canvas");
  await c.focus();
  await page.keyboard.press("s");
  await expect(
    page
      .locator(".scene-dock")
      .getByRole("button", { name: "Slice", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await page.keyboard.press("r");
  await expect(
    page
      .locator(".scene-dock")
      .getByRole("button", { name: "Rotate", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await page.keyboard.press("Space");
  await expect(
    page
      .locator(".scene-dock")
      .getByRole("button", { name: "Pause", exact: true }),
  ).toBeVisible();
  await page.keyboard.press("Space");
  await expect(
    page
      .locator(".scene-dock")
      .getByRole("button", { name: "Play", exact: true })
      .locator(".tool-pictogram"),
  ).toHaveCount(1);
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Scan", exact: true })
    .click();
  await page
    .getByRole("textbox", { name: "Experiment observation", exact: true })
    .fill("r s k");
  await expect(
    page.getByRole("textbox", { name: "Experiment observation", exact: true }),
  ).toHaveValue("r s k");
});

test("replacing a journey cancels any capture requested in the old one", async ({
  page,
}) => {
  await scene(page);
  let release!: () => void;
  const gate = new Promise<void>((r) => (release = r));
  await page.route("**/api/su2", async (route) => {
    const response = await route.fetch();
    await gate;
    try {
      await route.fulfill({ response });
    } catch {}
  });
  await page
    .getByRole("button", { name: "Change rotation axis", exact: true })
    .click();
  await page
    .locator(".scene-dock")
    .getByRole("button", { name: "Keep discovery", exact: true })
    .click();
  await expect(
    page
      .locator(".scene-dock")
      .getByRole("button", { name: "Keep discovery", exact: true }),
  ).toHaveAttribute("aria-busy", "true");
  const input = await page.evaluate(
    () => JSON.parse(localStorage.getItem("mathphysics.explore.v1")!).input,
  );
  await page
    .getByRole("button", { name: "Field journal", exact: true })
    .click();
  page.once("dialog", (d) => d.accept());
  await page
    .getByLabel("Import exploration file", { exact: true })
    .setInputFiles({
      name: "replacement.json",
      mimeType: "application/json",
      buffer: Buffer.from(
        JSON.stringify({
          version: 1,
          convention: "pauli-plus-i-v1",
          location: "lab",
          input,
          journal: [],
        }),
      ),
    });
  await expect(page.locator("dialog[open]")).toHaveCount(0);
  release();
  await expect(page.locator(".explore-evidence")).toHaveAttribute(
    "data-status",
    "checked",
  );
  await expect(page.locator(".explore-journal-entry")).toHaveCount(0);
  expect(
    await page
      .locator(".scene-dock")
      .getByRole("button", { name: "Keep discovery", exact: true })
      .getAttribute("aria-busy"),
  ).toBe(null);
});
