import test from 'node:test';
import assert from 'node:assert/strict';
import { wrenchFromKeys, HeldControls } from '../src/controls.js';

test('no input is exactly zero external wrench; releasing all inputs restores coast', () => {
  const held = new HeldControls();
  assert.deepEqual(held.wrench(25, 10), { force_B: [0, 0, 0], torque_B: [0, 0, 0], boost: false });
  held.hold('keyboard:KeyW', 'KeyW'); held.hold('keyboard:ArrowUp', 'ArrowUp');
  assert.deepEqual(held.wrench(25, 10).force_B, [25, 0, 0]);
  held.clear();
  assert.deepEqual(held.wrench(25, 10).torque_B, [0, 0, 0]);
});
test('each translation and attitude command follows the aerospace body-frame signs', () => {
  const cases = [
    ['KeyW', [25, 0, 0], [0, 0, 0]], ['KeyS', [-25, 0, 0], [0, 0, 0]],
    ['KeyD', [0, 25, 0], [0, 0, 0]], ['KeyA', [0, -25, 0], [0, 0, 0]],
    ['KeyF', [0, 0, 25], [0, 0, 0]], ['KeyR', [0, 0, -25], [0, 0, 0]],
    ['ArrowUp', [0, 0, 0], [0, 10, 0]], ['ArrowDown', [0, 0, 0], [0, -10, 0]],
    ['ArrowRight', [0, 0, 0], [0, 0, 10]], ['ArrowLeft', [0, 0, 0], [0, 0, -10]],
    ['KeyE', [0, 0, 0], [10, 0, 0]], ['KeyQ', [0, 0, 0], [-10, 0, 0]],
  ];
  for (const [key, force, torque] of cases) {
    const actual = wrenchFromKeys(new Set([key]));
    assert.deepEqual(actual.force_B, force, key); assert.deepEqual(actual.torque_B, torque, key);
  }
});
test('opposing commands cancel; boost affects translation only', () => {
  const actual = wrenchFromKeys(new Set(['KeyW', 'KeyS', 'KeyD', 'ArrowUp', 'ShiftLeft']));
  assert.deepEqual(actual.force_B, [0, 100, 0]); assert.deepEqual(actual.torque_B, [0, 10, 0]);
});
test('independent pointer and keyboard holds do not release one another', () => {
  const held = new HeldControls(); held.hold('keyboard:KeyW', 'KeyW'); held.hold('pointer:1', 'KeyW');
  held.release('pointer:1'); assert.equal(held.wrench(25, 10).force_B[0], 25);
  held.release('keyboard:KeyW'); assert.equal(held.wrench(25, 10).force_B[0], 0);
});
test('invalid scales and unknown controls cannot inject a wrench', () => {
  assert.throws(() => wrenchFromKeys(new Set(), Infinity, 10), RangeError);
  assert.throws(() => wrenchFromKeys(new Set(), 25, -1), RangeError);
  const held = new HeldControls(); held.hold('unknown', 'KeyZ'); assert.equal(held.keys.size, 0);
});
