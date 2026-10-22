import { test, expect } from '@playwright/test';

const state = page => page.evaluate(() => window.__simTest.state);
const norm = v => Math.hypot(...v);
const waitPaused = (page, value) => page.waitForFunction(v => window.__simTest.state?.paused === v, value);
const reset = async (page, preset) => {
  await page.getByTestId(`preset-${preset}`).click();
  await page.waitForFunction(p => window.__simTest.state?.preset === p && window.__simTest.state.sim_time < 1, preset);
  await expect(page.getByTestId('connection-status')).toContainText('CONNECTED');
};

test('authoritative browser controls, free coast, cameras and mobile touch', async ({ page, browser }) => {
  const browserErrors = [];
  page.on('pageerror', e => browserErrors.push(e.message));
  await page.request.post('/api/reset', { data: { preset: 'free-tumble' } });
  await page.goto('/');
  await page.waitForFunction(() => window.__simTest?.connected && window.__simTest.state && window.__simTest.webgl);
  await expect(page.locator('#viewport canvas')).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBeTruthy();
  await page.screenshot({ path: '/tmp/mathphysics-rigid-body-desktop.png', fullPage: true });

  await page.getByTestId('pause').click(); await waitPaused(page, true);
  const paused = await state(page);
  await page.waitForTimeout(200);
  expect((await state(page)).sim_time).toBe(paused.sim_time);
  await page.getByTestId('step').click();
  await page.waitForFunction(t => window.__simTest.state.sim_time > t, paused.sim_time);
  expect((await state(page)).sim_time - paused.sim_time).toBeCloseTo(1 / 240, 12);

  await reset(page, 'thruster-sandbox');
  await page.locator('#viewport').focus();
  await page.keyboard.down('w');
  await page.waitForFunction(() => window.__simTest.state.control.force_B[0] > 0);
  await page.waitForTimeout(350);
  await page.keyboard.up('w');
  await page.waitForFunction(() => window.__simTest.state.control.force_B.every(x => x === 0));
  const coast = await state(page);
  expect(coast.v_I[0]).toBeGreaterThan(.04);
  await page.waitForTimeout(250);
  const later = await state(page);
  expect(later.v_I).toEqual(coast.v_I);
  expect(later.r_I[0]).toBeGreaterThan(coast.r_I[0]);
  expect(later.diagnostics.free_segment).toBeTruthy();
  expect(later.diagnostics.energy_relative_drift).toBeLessThan(1e-10);

  await page.keyboard.down('ArrowUp');
  await page.waitForFunction(() => window.__simTest.state.control.torque_B[1] > 0);
  await page.waitForTimeout(200); await page.keyboard.up('ArrowUp');
  await page.waitForFunction(() => window.__simTest.state.control.torque_B.every(x => x === 0));
  expect((await state(page)).omega_B[1]).toBeGreaterThan(.01);
  expect(Math.abs((await state(page)).q[2])).toBeGreaterThan(0);

  // A lost window focus must release a held command before the next step.
  await page.keyboard.down('w');
  await page.waitForFunction(() => window.__simTest.keys.includes('KeyW'));
  await page.evaluate(() => window.dispatchEvent(new Event('blur')));
  await page.waitForFunction(() => !window.__simTest.keys.length && window.__simTest.wrench.force_B.every(x => x === 0));
  await page.keyboard.up('w');
  await page.getByTestId('mass').focus();
  await page.keyboard.down('w');
  expect(await page.evaluate(() => window.__simTest.keys.length)).toBe(0);
  await page.keyboard.up('w');
  await page.getByTestId('mass').fill('120');
  await page.getByTestId('omega-x').fill('0.4');
  await page.getByTestId('apply-initial').click();
  await page.waitForFunction(() => window.__simTest.state.mass === 120 && Math.abs(window.__simTest.state.omega_B[0] - .4) < .005);
  expect((await state(page)).omega_B[0]).toBeCloseTo(.4, 2);
  await page.locator('[data-camera="chase"]').click();
  expect(await page.evaluate(() => window.__simTest.cameraMode)).toBe('chase');
  await page.locator('[data-camera="orbit"]').click();
  await page.getByTestId('speed').selectOption('5');
  await page.waitForFunction(() => window.__simTest.state.speed === 5);
  await reset(page, 'principal-spin');
  expect((await state(page)).omega_B.slice(0, 2)).toEqual([0, 0]);
  await reset(page, 'intermediate-axis-flip');
  expect((await state(page)).omega_B[1]).toBeGreaterThan(.9);
  await reset(page, 'free-tumble');
  expect(norm((await state(page)).omega_B)).toBeGreaterThan(0);
  expect(await page.evaluate(() => window.__simTest.errors)).toEqual([]);
  expect(browserErrors).toEqual([]);
  await page.screenshot({ path: '/tmp/mathphysics-rigid-body-desktop.png', fullPage: true });
  await page.close();

  const mobile = await browser.newContext({ viewport: { width: 390, height: 844 }, hasTouch: true, isMobile: true, deviceScaleFactor: 1 });
  const touch = await mobile.newPage();
  touch.on('pageerror', e => browserErrors.push(e.message));
  await touch.goto('/');
  await touch.waitForFunction(() => window.__simTest?.connected && window.__simTest.webgl);
  expect(await touch.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBeTruthy();
  await reset(touch, 'thruster-sandbox');
  const forward = touch.getByTestId('thrust-forward');
  await forward.scrollIntoViewIfNeeded();
  const box = await forward.boundingBox();
  const cdp = await mobile.newCDPSession(touch);
  await cdp.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x: box.x + box.width / 2, y: box.y + box.height / 2 }] });
  await touch.waitForFunction(() => window.__simTest.state.control.force_B[0] > 0);
  await touch.waitForTimeout(200);
  await cdp.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  await touch.waitForFunction(() => window.__simTest.state.control.force_B[0] === 0);
  expect((await state(touch)).v_I[0]).toBeGreaterThan(.02);
  await touch.screenshot({ path: '/tmp/mathphysics-rigid-body-mobile.png', fullPage: true });
  expect(browserErrors).toEqual([]);
  expect(await touch.evaluate(() => window.__simTest.errors)).toEqual([]);
  await mobile.close();
});
