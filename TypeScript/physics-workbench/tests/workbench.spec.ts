import { test, expect } from "@playwright/test";
test.beforeEach(async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#connection")).toHaveText("● Local workspace");
});
test("styles switch live without replacing the editor and persist its draft", async ({
  page,
}) => {
  const input = page.locator('[data-cell-id="first-trace"] textarea');
  await input.fill("fresh_mu, fresh_nu");
  await input.evaluate((e) => e.setAttribute("data-original", "yes"));
  const size = await input.evaluate((e) => getComputedStyle(e).fontSize);
  await page.locator("#style-pixel").click();
  await expect(page.locator("html")).toHaveAttribute("data-style", "pixel");
  await expect(input).toHaveAttribute("data-original", "yes");
  await expect(input).toHaveValue("fresh_mu, fresh_nu");
  expect(await input.evaluate((e) => getComputedStyle(e).fontSize)).toBe(size);
  await page.reload();
  await expect(input).toHaveValue("fresh_mu, fresh_nu");
  await expect(page.locator("html")).toHaveAttribute("data-style", "pixel");
  await page.locator("#style-notebook").click();
  await expect(input).toHaveValue("fresh_mu, fresh_nu");
});
test("real Rust calculation and invalid inputs retain truthful status", async ({
  page,
}) => {
  const cell = page.locator('[data-cell-id="first-trace"]');
  await cell.locator("textarea").fill("mu, nu");
  await cell.getByRole("button", { name: "Run cell ↵" }).click();
  await expect(
    cell.getByText("Computed · exact algebra", { exact: true }),
  ).toBeVisible();
  await expect(cell.locator(".trace-math .katex")).toBeVisible();
  await cell.locator("textarea").fill("invalid<symbol, nu");
  await cell.getByRole("button", { name: "Run cell ↵" }).click();
  await expect(cell.locator(".status-pill")).toHaveText("invalid input");
  await expect(
    cell.getByText("Computed · exact algebra", { exact: true }),
  ).toHaveCount(0);
});
test("search, open TeX, paired PDF and preserve PDF element during style change", async ({
  page,
}) => {
  await page.locator("#root").selectOption("fixtures");
  await page.locator("#search").fill("sample.tex");
  const file = page.locator(".file-entry").filter({ hasText: "sample.tex" });
  await expect(file).toHaveCount(1);
  await file.click();
  await expect(page.locator(".source-editor")).toHaveValue(/Spinor notes/);
  await expect(page.locator(".latex-preview .katex")).toHaveCount(1);
  await page
    .getByRole("button", { name: "Open compiled PDF", exact: true })
    .click();
  const frame = page.locator(".pdf-page[data-rendered=true]");
  await expect(frame).toBeVisible();
  await expect(frame.locator(".textLayer")).toContainText(
    "Mathphysics PDF reading test",
  );
  expect(
    await frame.locator("canvas").evaluate((c: HTMLCanvasElement) => {
      const d = c.getContext("2d")!.getImageData(0, 0, c.width, c.height).data;
      return d.some((v, i) => i % 4 !== 3 && v < 100);
    }),
  ).toBe(true);
  const response = await page.request.get(
    "/api/pdf?root=fixtures&path=sample.pdf",
    { headers: { Range: "bytes=0-9" } },
  );
  expect(response.status()).toBe(206);
  expect(response.headers()["content-type"]).toContain("application/pdf");
  expect((await response.body()).toString()).toContain("%PDF");
  await frame.evaluate((e) => e.setAttribute("data-original", "pdf"));
  await page.locator("#style-pixel").click();
  await expect(frame).toHaveAttribute("data-original", "pdf");
  await expect(
    page.getByRole("button", { name: "Open TeX source", exact: true }),
  ).toBeVisible();
  await page
    .getByRole("button", { name: "Open TeX source", exact: true })
    .click();
  await expect(page.locator(".source-editor")).toHaveValue(/Spinor notes/);
});
test("favorites survive clearing automatic history", async ({ page }) => {
  await page.locator("#root").selectOption("fixtures");
  await expect(page.locator("#path")).toHaveValue("");
  await page.locator("#favorite").click();
  await expect(page.locator("#favorite")).toHaveText("★");
  await page
    .getByRole("button", { name: "Clear visit history", exact: true })
    .click();
  await expect(page.locator("#favorite")).toHaveText("★");
  await page.reload();
  await expect(page.locator("#favorite")).toHaveText("★");
});
test("local Markdown renders math while HTML handlers and code math stay inert", async ({
  page,
}) => {
  await page.locator("#local-input").setInputFiles({
    name: "local-note.md",
    mimeType: "text/markdown",
    buffer: Buffer.from(
      '# Local\n$E=mc^2$\n\n`$not_math$`\n<img src=x onerror="window.XSS=1"><script>window.XSS=1</script>',
    ),
  });
  await expect(page.locator(".latex-preview .katex")).toHaveCount(1);
  await expect(page.locator(".latex-preview code")).toHaveText("$not_math$");
  expect(await page.evaluate(() => "XSS" in window)).toBe(false);
  expect(await page.locator(".latex-preview [onerror]").count()).toBe(0);
});
test("notebook exports the exact edited input and sound preferences persist", async ({
  page,
}) => {
  await page
    .locator('[data-cell-id="first-trace"] textarea')
    .fill("export_mu, export_nu");
  await page.locator("#sound").check();
  await page.locator("#volume").fill("0.4");
  const download = page.waitForEvent("download");
  await page
    .getByRole("button", { name: "Export notebook", exact: true })
    .click();
  const file = await download;
  expect(file.suggestedFilename()).toBe("research.mathphysics.json");
  const path = await file.path();
  const fs = await import("node:fs/promises");
  const doc = JSON.parse(await fs.readFile(path!, "utf8"));
  expect(
    doc.cells.find((c: { id: string }) => c.id === "first-trace").source,
  ).toBe("export_mu, export_nu");
  await page.reload();
  await expect(page.locator("#sound")).toBeChecked();
  await expect(page.locator("#volume")).toHaveValue("0.4");
});
test("mobile file drawer and keyboard search remain operable without page overflow", async ({
  page,
}) => {
  await page.setViewportSize({ width: 460, height: 900 });
  await page.locator("#mobile-files").click();
  await expect(page.locator("html")).toHaveAttribute("data-sidebar", "open");
  await expect(page.locator("#search")).toBeFocused();
  await page.keyboard.press("Escape");
  await expect(page.locator("html")).toHaveAttribute("data-sidebar", "closed");
  await page.keyboard.press("Control+k");
  await expect(page.locator("#search")).toBeFocused();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
});

test("local PDF renders without uploading its bytes", async ({ page }) => {
  const response = await page.request.get(
    "/api/pdf?root=fixtures&path=sample.pdf",
  );
  await page
    .locator("#local-input")
    .setInputFiles({
      name: "local.pdf",
      mimeType: "application/pdf",
      buffer: await response.body(),
    });
  await expect(
    page.locator(".pdf-page[data-rendered=true] .textLayer"),
  ).toContainText("Mathphysics PDF reading test");
});
