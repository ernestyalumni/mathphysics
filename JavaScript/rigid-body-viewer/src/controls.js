/** Body +X nose, +Y starboard, +Z belly. Holding a key applies wrench;
 * releasing it applies exactly zero wrench. There is no automatic damping. */
export const CONTROL_CODES = new Set([
  'KeyW', 'KeyS', 'KeyA', 'KeyD', 'KeyR', 'KeyF',
  'ArrowUp', 'ArrowDown', 'ArrowRight', 'ArrowLeft', 'KeyE', 'KeyQ',
  'ShiftLeft', 'ShiftRight',
]);
export function wrenchFromKeys(keys, force = 25, torque = 10) {
  if (![force, torque].every(v => Number.isFinite(v) && v >= 0)) throw new RangeError('finite nonnegative wrench scale required');
  const axis = (positive, negative) => Number(keys.has(positive)) - Number(keys.has(negative));
  const boost = keys.has('ShiftLeft') || keys.has('ShiftRight');
  const f = force * (boost ? 4 : 1);
  return {
    force_B: [axis('KeyW', 'KeyS') * f, axis('KeyD', 'KeyA') * f, axis('KeyF', 'KeyR') * f],
    torque_B: [axis('KeyE', 'KeyQ') * torque, axis('ArrowUp', 'ArrowDown') * torque, axis('ArrowRight', 'ArrowLeft') * torque],
    boost,
  };
}
export function isEditing(target) {
  return Boolean(target?.closest?.('input,select,textarea,[contenteditable="true"]'));
}
/** Distinguishes physical keys and pointer buttons, so releasing one source
 * cannot cancel another still-held source. */
export class HeldControls {
  constructor() { this.sources = new Map(); }
  hold(source, code) { if (CONTROL_CODES.has(code)) this.sources.set(source, code); }
  release(source) { this.sources.delete(source); }
  clear() { this.sources.clear(); }
  get keys() { return new Set(this.sources.values()); }
  wrench(force, torque) { return wrenchFromKeys(this.keys, force, torque); }
}
