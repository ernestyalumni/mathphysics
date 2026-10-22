import './style.css';
import { FlightScene } from './scene.js';
import { CONTROL_CODES, HeldControls, isEditing } from './controls.js';

const $ = selector => document.querySelector(selector);
const $$ = selector => [...document.querySelectorAll(selector)];
const held = new HeldControls();
const PRESETS = {
  'free-tumble': ['Free tumble', 'Angular momentum stays fixed. The body does not.'],
  'principal-spin': ['Principal-axis spin', 'One principal axis. A steady, torque-free spin.'],
  'intermediate-axis-flip': ['Intermediate-axis flip', 'A small perturbation. A remarkable change in attitude.'],
  'thruster-sandbox': ['Thruster sandbox', 'Your wrench. Your trajectory. Release the controls to coast.'],
};
let scene;
try { scene = new FlightScene($('#viewport')); }
catch (error) {
  $('#scene-error').hidden = false;
  $('#scene-error').textContent = `The 3D view needs WebGL. ${error.message}`;
}
let snapshot = null, desiredPaused = false, speed = 1, preset = 'free-tumble';
let resetting = false, connected = false, requestCount = 0, lastRequest = null;
let queue = Promise.resolve(), networkPending = false, apiErrors = [];
let graphHistory = [], lastGraphTime = null, noticeTimer;
let lastAdvance = performance.now(), nextRetry = 0;
const zero = () => [0, 0, 0];
const norm = values => Math.hypot(...(values || zero()));
const value = id => Number($(id).value);
const format = (n, digits = 3) => Number.isFinite(n) ? n.toLocaleString('en-US', { minimumFractionDigits: digits, maximumFractionDigits: digits, useGrouping: false }) : '—';
const vectorText = values => (values || zero()).map(x => format(x, 2)).join(' / ');
const driftText = n => Number.isFinite(n) ? Math.abs(n).toExponential(2) : '—';

function notice(message, error = false) {
  $('#notice').textContent = message; $('#notice').classList.toggle('error', error); $('#notice').hidden = false;
  clearTimeout(noticeTimer); noticeTimer = setTimeout(() => { $('#notice').hidden = true; }, error ? 6500 : 2600);
}
function connection(online) {
  connected = online; $('#connection').classList.toggle('online', online);
  $('#connection span').textContent = online ? 'LOCAL SOLVER CONNECTED' : 'SOLVER DISCONNECTED';
}
async function request(path, payload) {
  lastRequest = { path, payload }; requestCount++;
  const response = await fetch(path, {
    method: payload ? 'POST' : 'GET', headers: { 'Content-Type': 'application/json' },
    ...(payload ? { body: JSON.stringify(payload), keepalive: true } : {}),
  });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `Solver returned HTTP ${response.status}`);
  connection(true); return data;
}
function serial(action) {
  const next = queue.then(action, action); queue = next.catch(() => {}); return next;
}
function apiError(error, showNotice = true) {
  connection(false); apiErrors.push(error.message); apiErrors = apiErrors.slice(-20); nextRetry = performance.now() + 1500;
  if (showNotice) notice(`Solver: ${error.message}`, true);
}
function rawWrench() { return held.wrench(value('#force-scale'), value('#torque-scale')); }
function currentWrench() {
  const wrench = rawWrench();
  return desiredPaused || resetting ? { ...wrench, force_B: zero(), torque_B: zero() } : wrench;
}
function advancePayload(overrides = {}) {
  const wrench = currentWrench();
  return { wall_dt: 0, force_B: wrench.force_B, torque_B: wrench.torque_B, speed, paused: desiredPaused, ...overrides };
}
async function explicitAdvance(overrides = {}) {
  const payload = advancePayload(overrides);
  const state = await serial(() => request('/api/advance', payload)); adopt(state); return state;
}
function releaseAll(sendZero = false) {
  held.clear(); updateControlDisplay();
  if (sendZero && connected) explicitAdvance({ wall_dt: 0, force_B: zero(), torque_B: zero() }).catch(error => apiError(error, false));
}
function configureFields(state) {
  $('#mass').value = state.mass ?? 100;
  ['x', 'y', 'z'].forEach((axis, i) => {
    $(`#dimension-${axis}`).value = state.dimensions?.[i] ?? [3, 2, 1][i];
    $(`#omega-${axis}`).value = format(state.omega_B?.[i] ?? 0, 5);
  });
}
function selectPreset(next) {
  preset = PRESETS[next] ? next : 'free-tumble';
  $$('.preset').forEach(button => {
    const selected = button.dataset.preset === preset; button.classList.toggle('selected', selected); button.setAttribute('aria-pressed', String(selected));
  });
  $('#experiment-title').textContent = PRESETS[preset][0]; $('#experiment-description').textContent = PRESETS[preset][1];
}
async function resetExperiment(nextPreset = preset, custom = false) {
  if (resetting) return;
  resetting = true; releaseAll(); selectPreset(nextPreset);
  const payload = {
    preset: nextPreset, mass: value('#mass'), dimensions: ['x', 'y', 'z'].map(axis => value(`#dimension-${axis}`)),
    ...(custom ? { omega_B: ['x', 'y', 'z'].map(axis => value(`#omega-${axis}`)) } : {}),
  };
  try {
    const state = await serial(() => request('/api/reset', payload));
    desiredPaused = Boolean(state.paused); speed = Number(state.speed ?? 1); $('#speed').value = String(speed);
    graphHistory = []; lastGraphTime = null; scene?.clearTrajectory(); adopt(state); configureFields(state); scene?.recenter();
    notice(`${PRESETS[nextPreset][0]} ready. Hold a flight key to apply thrust.`);
  } catch (error) { apiError(error); }
  finally { resetting = false; lastAdvance = performance.now(); updateControlDisplay(); }
}
function setCamera(mode) {
  scene?.setMode(mode);
  $$('[data-camera]').forEach(button => { const active = button.dataset.camera === mode; button.classList.toggle('active', active); button.setAttribute('aria-pressed', String(active)); });
}
function updateControlDisplay() {
  const keys = held.keys;
  $$('.flight-key').forEach(button => { const active = keys.has(button.dataset.code); button.classList.toggle('held', active); button.setAttribute('aria-pressed', String(active)); });
  $('#boost-indicator').classList.toggle('held', keys.has('ShiftLeft') || keys.has('ShiftRight'));
  $('#force-value').value = $('#force-scale').value; $('#torque-value').value = $('#torque-scale').value;
  $('#pause-label').textContent = desiredPaused ? 'Resume' : 'Pause'; $('#pause-icon').textContent = desiredPaused ? '▶' : 'Ⅱ';
  $('#pause').setAttribute('aria-label', desiredPaused ? 'Resume simulation' : 'Pause simulation');
  $('#step').disabled = !desiredPaused || resetting;
}
function drawGraph(omega, time) {
  if (time === lastGraphTime) return;
  if (lastGraphTime !== null && time < lastGraphTime) graphHistory = [];
  lastGraphTime = time; graphHistory.push([...omega]); if (graphHistory.length > 150) graphHistory.shift();
  const limit = Math.max(0.5, ...graphHistory.flat().map(Math.abs));
  ['x', 'y', 'z'].forEach((axis, component) => {
    const path = graphHistory.map((v, i) => `${i ? 'L' : 'M'}${(i / 149 * 240).toFixed(2)},${(33 - v[component] / limit * 27).toFixed(2)}`).join(' ');
    $(`#graph-${axis}`).setAttribute('d', path);
  });
  $('#rate-norm').textContent = format(norm(omega));
}
function adopt(state) {
  if (!state || !Array.isArray(state.q) || !Array.isArray(state.r_I)) throw new Error('Unexpected solver state');
  snapshot = state; scene?.setState(state);
  const time = Math.max(0, state.sim_time || 0), minutes = Math.floor(time / 60), seconds = time % 60;
  $('#sim-time').textContent = `${String(minutes).padStart(2, '0')}:${seconds.toFixed(3).padStart(6, '0')}`;
  $('#position-value').textContent = vectorText(state.r_I); $('#velocity-value').textContent = vectorText(state.v_I); $('#omega-value').textContent = vectorText(state.omega_B);
  const d = state.diagnostics || {};
  $('#energy-value').textContent = format(d.energy, 3); $('#angular-value').textContent = format(norm(d.angular_momentum_I), 3); $('#linear-value').textContent = format(norm(d.linear_momentum_I), 3);
  $('#energy-drift').textContent = driftText(d.energy_relative_drift); $('#angular-drift').textContent = driftText(d.angular_momentum_relative_drift); $('#linear-drift').textContent = driftText(d.linear_momentum_relative_drift);
  $('#quaternion-value').textContent = state.q.map(x => format(x, 3)).join(' / '); $('#qnorm-value').textContent = format(d.qnorm ?? norm(state.q), 10);
  const force = state.control?.force_B || zero(), torque = state.control?.torque_B || zero();
  const firing = norm(force) > 0 || norm(torque) > 0;
  $('#wrench-status').classList.toggle('firing', firing);
  $('#wrench-status>span').textContent = state.paused ? 'PAUSED' : firing ? 'THRUST APPLIED' : 'COASTING';
  $('#wrench-status small').textContent = firing ? `|F| ${format(norm(force), 0)} N · |τ| ${format(norm(torque), 1)} N·m` : 'ZERO EXTERNAL WRENCH';
  $('#reference-label').textContent = firing ? 'EXTERNAL WRENCH ACTIVE' : d.free_segment === false ? 'LAST DRIVEN SEGMENT' : (d.reference || 'coast segment').toUpperCase();
  $('#drift-explanation').textContent = firing ? 'Energy and momentum may change under thrust' : d.free_segment === false ? 'Reference precedes the last applied thrust' : 'Drift scaled by max(reference, 1 SI unit)';
  for (const id of ['#energy-drift', '#angular-drift', '#linear-drift']) $(id).style.color = firing || d.free_segment === false ? '#b6946c' : '';
  drawGraph(state.omega_B || zero(), time); updateControlDisplay();
}
async function togglePause() {
  desiredPaused = !desiredPaused; releaseAll(); updateControlDisplay();
  try { await explicitAdvance({ wall_dt: 0, force_B: zero(), torque_B: zero() }); }
  catch (error) { apiError(error); }
  lastAdvance = performance.now();
}
async function stepOnce() {
  if (!desiredPaused || resetting) return;
  releaseAll();
  try { await explicitAdvance({ wall_dt: 0, step: true, paused: true, force_B: zero(), torque_B: zero() }); }
  catch (error) { apiError(error); }
}

$$('.preset').forEach(button => button.addEventListener('click', () => resetExperiment(button.dataset.preset)));
$('#initial-form').addEventListener('submit', event => { event.preventDefault(); if (event.currentTarget.reportValidity()) resetExperiment(preset, true); });
$('#pause').addEventListener('click', togglePause); $('#step').addEventListener('click', stepOnce);
$('#reset').addEventListener('click', () => resetExperiment(preset));
$('#speed').addEventListener('change', () => { speed = value('#speed'); explicitAdvance({ wall_dt: 0 }).catch(apiError); });
$('#force-scale').addEventListener('input', updateControlDisplay); $('#torque-scale').addEventListener('input', updateControlDisplay);
$$('[data-camera]').forEach(button => button.addEventListener('click', () => setCamera(button.dataset.camera)));
$('#recenter').addEventListener('click', () => scene?.recenter());
$('#show-grid').addEventListener('change', event => { if (scene) scene.reference.visible = event.target.checked; });
$('#show-momentum').addEventListener('change', event => { if (scene) { scene.showMomentum = event.target.checked; if (snapshot) scene.setState(snapshot); } });
$$('.flight-key').forEach(button => {
  button.addEventListener('pointerdown', event => {
    if (event.button !== 0) return;
    event.preventDefault(); button.setPointerCapture(event.pointerId); held.hold(`pointer:${event.pointerId}`, button.dataset.code); updateControlDisplay();
  });
  for (const eventName of ['pointerup', 'pointercancel', 'lostpointercapture']) button.addEventListener(eventName, event => { held.release(`pointer:${event.pointerId}`); updateControlDisplay(); });
  // Enter activates the same momentary command for keyboard-only users.
  button.addEventListener('keydown', event => { if (event.code === 'Enter') { event.preventDefault(); held.hold(`button:${button.dataset.code}`, button.dataset.code); updateControlDisplay(); } });
  button.addEventListener('keyup', event => { if (event.code === 'Enter') { held.release(`button:${button.dataset.code}`); updateControlDisplay(); } });
  button.addEventListener('blur', () => { held.release(`button:${button.dataset.code}`); updateControlDisplay(); });
});
window.addEventListener('keydown', event => {
  if (isEditing(event.target) || event.metaKey || event.ctrlKey || event.altKey) return;
  if (CONTROL_CODES.has(event.code)) { event.preventDefault(); held.hold(`keyboard:${event.code}`, event.code); updateControlDisplay(); }
  else if (!event.repeat && event.code === 'Space') { event.preventDefault(); togglePause(); }
  else if (!event.repeat && event.code === 'Period') { event.preventDefault(); stepOnce(); }
  else if (!event.repeat && event.code === 'KeyC') { event.preventDefault(); setCamera(scene?.mode === 'chase' ? 'orbit' : 'chase'); }
});
window.addEventListener('keyup', event => { held.release(`keyboard:${event.code}`); updateControlDisplay(); });
window.addEventListener('blur', () => releaseAll(true));
document.addEventListener('visibilitychange', () => { if (document.hidden) releaseAll(true); lastAdvance = performance.now(); });
window.addEventListener('pagehide', () => {
  held.clear();
  if (connected) navigator.sendBeacon('/api/advance', new Blob([JSON.stringify(advancePayload({ wall_dt: 0, force_B: zero(), torque_B: zero() }))], { type: 'application/json' }));
});

async function tick() {
  if (networkPending || resetting || document.hidden || performance.now() < nextRetry) return;
  const now = performance.now(); networkPending = true;
  try {
    if (!snapshot || !connected) {
      const state = await serial(() => request('/api/state'));
      if (!snapshot) { desiredPaused = Boolean(state.paused); speed = state.speed || 1; $('#speed').value = String(speed); selectPreset(state.preset); configureFields(state); }
      adopt(state);
    } else {
      const wall_dt = Math.min(0.1, Math.max(0, (now - lastAdvance) / 1000));
      await explicitAdvance({ wall_dt });
    }
  } catch (error) { apiError(error, connected || !snapshot); }
  finally { lastAdvance = now; networkPending = false; }
}
let lastFrame = performance.now();
function frame(now) { const dt = Math.min((now - lastFrame) / 1000, 0.1); lastFrame = now; scene?.render(dt); requestAnimationFrame(frame); }
requestAnimationFrame(frame); setInterval(tick, 1000 / 30); tick(); updateControlDisplay();

// Read-only observation hooks for browser verification. They never run physics
// in JavaScript; the HTTP solver remains authoritative.
window.__simTest = {
  get state() { return snapshot; }, get keys() { return [...held.keys]; }, get requestCount() { return requestCount; },
  get lastRequest() { return lastRequest; }, get errors() { return [...apiErrors]; }, get connected() { return connected; },
  get webgl() { return Boolean(scene?.renderer); }, get cameraMode() { return scene?.mode; },
  get wrench() { return currentWrench(); }, get paused() { return desiredPaused; },
};
