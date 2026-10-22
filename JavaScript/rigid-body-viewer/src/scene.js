import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const AXIS_COLORS = [0xfa897e, 0x6ae2bc, 0x91abff];
const UNIT_AXES = [new THREE.Vector3(1, 0, 0), new THREE.Vector3(0, 1, 0), new THREE.Vector3(0, 0, 1)];
const vector = (v = [0, 0, 0]) => new THREE.Vector3(...v);
function disposeTree(root) {
  root.traverse(object => {
    object.geometry?.dispose();
    const materials = Array.isArray(object.material) ? object.material : [object.material];
    for (const material of materials) { material?.map?.dispose(); material?.dispose(); }
  });
}
function label(text, color) {
  const canvas = document.createElement('canvas'); canvas.width = 192; canvas.height = 96;
  const ctx = canvas.getContext('2d');
  ctx.font = '44px ui-monospace, monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillStyle = `#${new THREE.Color(color).getHexString()}`; ctx.fillText(text, 96, 48);
  const texture = new THREE.CanvasTexture(canvas);
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, transparent: true, opacity: 0.85, depthTest: false }));
  sprite.scale.set(1.15, 0.575, 1); return sprite;
}
function arrow(direction, color) {
  return new THREE.ArrowHelper(direction, new THREE.Vector3(), 1, color, 0.25, 0.12);
}
export class FlightScene {
  constructor(container) {
    this.container = container;
    this.scene = new THREE.Scene(); this.scene.fog = new THREE.FogExp2(0x080f19, 0.0015);
    this.camera = new THREE.PerspectiveCamera(43, 1, 0.02, 3000);
    this.camera.up.set(0, 0, -1); this.camera.position.set(-8.5, -10.5, -7.5);
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.setClearColor(0x080f19, 0); this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping; this.renderer.toneMappingExposure = 1.25;
    container.append(this.renderer.domElement);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true; this.controls.dampingFactor = 0.075;
    this.controls.minDistance = 4; this.controls.maxDistance = 120; this.controls.enablePan = true;
    this.controls.target.set(0, 0, 0); this.controls.update();
    this.mode = 'orbit'; this.position = new THREE.Vector3(); this.orientation = new THREE.Quaternion();
    this.previousPosition = new THREE.Vector3(); this.chaseDistance = 13;
    this.scene.add(new THREE.HemisphereLight(0xb9dbf9, 0x132139, 2.2));
    const sun = new THREE.DirectionalLight(0xf0f6ff, 3.4); sun.position.set(6, -8, -11); this.scene.add(sun);
    const rim = new THREE.DirectionalLight(0x5fe2c4, 1.8); rim.position.set(-6, 7, 3); this.scene.add(rim);
    this.body = new THREE.Group(); this.scene.add(this.body);
    this.hull = null; this.dimensionKey = ''; this.dimensions = [3, 2, 1];
    this.buildHull(this.dimensions);
    this.reference = new THREE.Group(); this.scene.add(this.reference);
    const grid = new THREE.GridHelper(160, 80, 0x315a6b, 0x233e50);
    grid.rotation.x = Math.PI / 2; grid.material.transparent = true; grid.material.opacity = 0.23;
    this.reference.add(grid);
    for (let i = 0; i < 3; i++) {
      const a = arrow(UNIT_AXES[i], AXIS_COLORS[i]); a.setLength(6, 0.2, 0.1);
      a.line.material.transparent = true; a.line.material.opacity = 0.22;
      a.cone.material.transparent = true; a.cone.material.opacity = 0.25; this.reference.add(a);
    }
    this.velocityArrow = arrow(UNIT_AXES[0], 0x66cee0); this.scene.add(this.velocityArrow);
    this.momentumArrow = arrow(UNIT_AXES[0], 0xeac285); this.scene.add(this.momentumArrow);
    this.forceArrow = arrow(UNIT_AXES[0], 0xffbd75); this.forceArrow.visible = false; this.body.add(this.forceArrow);
    this.torqueArrow = arrow(UNIT_AXES[0], 0xef9d64); this.torqueArrow.visible = false; this.body.add(this.torqueArrow);
    this.plumes = new THREE.Group(); this.body.add(this.plumes);
    this.makePlumes();
    this.trajectoryPoints = []; this.maxTrajectory = 2400;
    this.trajectoryGeometry = new THREE.BufferGeometry();
    this.trajectoryGeometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(this.maxTrajectory * 3), 3));
    this.trajectoryGeometry.setDrawRange(0, 0);
    this.trajectory = new THREE.Line(this.trajectoryGeometry, new THREE.LineBasicMaterial({ color: 0x5aafa8, transparent: true, opacity: 0.55 }));
    this.trajectory.frustumCulled = false; this.scene.add(this.trajectory); this.previousTime = null;
    this.addStars();
    this.resizeObserver = new ResizeObserver(() => this.resize()); this.resizeObserver.observe(container); this.resize();
    container.addEventListener('wheel', event => {
      if (this.mode === 'chase') this.chaseDistance = THREE.MathUtils.clamp(this.chaseDistance * Math.exp(event.deltaY * 0.001), 5, 80);
    }, { passive: true });
  }
  buildHull(dimensions) {
    if (this.hull) { this.body.remove(this.hull); disposeTree(this.hull); }
    this.hull = new THREE.Group(); this.body.add(this.hull); this.dimensions = dimensions;
    const [x, y, z] = dimensions;
    const geometry = new THREE.BoxGeometry(x, y, z);
    const materials = [0x2e505d, 0x1a2e40, 0x2d4356, 0x243c50, 0x172a3c, 0x3a5367].map(color => new THREE.MeshStandardMaterial({ color, metalness: 0.65, roughness: 0.35 }));
    this.hull.add(new THREE.Mesh(geometry, materials));
    const edges = new THREE.LineSegments(new THREE.EdgesGeometry(geometry), new THREE.LineBasicMaterial({ color: 0xa1c9d8, transparent: true, opacity: 0.65 }));
    this.hull.add(edges);
    const panel = (size, position, color, glow = 0) => {
      const mesh = new THREE.Mesh(new THREE.BoxGeometry(...size), new THREE.MeshStandardMaterial({ color, roughness: 0.5, metalness: 0.5, emissive: glow, emissiveIntensity: 0.5 }));
      mesh.position.set(...position); this.hull.add(mesh); return mesh;
    };
    // Surface details stay close to the cuboid whose uniform inertia is simulated.
    panel([x * 0.42, y * 0.72, 0.018], [-x * 0.12, 0, -z / 2 - 0.015], 0x172735);
    panel([x * 0.43, 0.016, z * 0.22], [x * 0.05, -y / 2 - 0.01, 0], 0xc5a477);
    panel([x * 0.43, 0.016, z * 0.22], [x * 0.05, y / 2 + 0.01, 0], 0xc5a477);
    panel([0.018, y * 0.58, z * 0.54], [x / 2 + 0.013, 0, 0], 0x467e7b, 0x3fbb9d);
    for (let i = -2; i <= 2; i++) panel([0.012, y * 0.61, 0.02], [-x * 0.12 + i * x * 0.075, 0, -z / 2 - 0.027], 0x71899a);
    for (const sign of [-1, 1]) panel([x * 0.13, y * 0.08, z * 0.12], [-x * 0.34, sign * y / 2, -z * 0.28], 0x121e2b);
    const max = Math.max(...dimensions); this.bodyScale = max;
    this.controls.minDistance = max * 1.4; this.controls.maxDistance = Math.max(120, max * 60);
    for (let i = 0; i < 3; i++) {
      const length = dimensions[i] / 2 + max * 0.57;
      const a = arrow(UNIT_AXES[i], AXIS_COLORS[i]); a.setLength(length, max * 0.075, max * 0.035); this.hull.add(a);
      const text = label(['X_B', 'Y_B', 'Z_B'][i], AXIS_COLORS[i]); text.position.copy(UNIT_AXES[i]).multiplyScalar(length + max * 0.13); text.scale.multiplyScalar(max / 3); this.hull.add(text);
    }
    this.dimensionKey = dimensions.join(',');
  }
  makePlumes() {
    this.plumeObjects = [];
    for (let i = 0; i < 3; i++) for (const sign of [-1, 1]) {
      const mesh = new THREE.Mesh(new THREE.ConeGeometry(0.18, 1.3, 12, 1, true), new THREE.MeshBasicMaterial({ color: 0x72def1, transparent: true, opacity: 0.65, blending: THREE.AdditiveBlending, depthWrite: false, side: THREE.DoubleSide }));
      const direction = UNIT_AXES[i].clone().multiplyScalar(sign);
      mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, -1, 0), direction); mesh.visible = false;
      this.plumes.add(mesh); this.plumeObjects.push({ mesh, axis: i, sign });
    }
  }
  addStars() {
    let seed = 77123;
    const random = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
    const points = [], colors = [];
    for (let i = 0; i < 1900; i++) {
      const z = random() * 2 - 1, angle = random() * Math.PI * 2, radius = 300 + random() * 900;
      const xy = Math.sqrt(1 - z * z); points.push(radius * xy * Math.cos(angle), radius * xy * Math.sin(angle), radius * z);
      const b = 0.25 + random() * 0.55; colors.push(b * 0.83, b * 0.94, b);
    }
    const geometry = new THREE.BufferGeometry(); geometry.setAttribute('position', new THREE.Float32BufferAttribute(points, 3)); geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    this.stars = new THREE.Points(geometry, new THREE.PointsMaterial({ size: 0.8, sizeAttenuation: false, transparent: true, opacity: 0.8, vertexColors: true, depthWrite: false })); this.scene.add(this.stars);
  }
  resize() {
    const width = this.container.clientWidth, height = this.container.clientHeight;
    if (!width || !height) return;
    this.renderer.setSize(width, height); this.camera.aspect = width / height; this.camera.updateProjectionMatrix();
  }
  setMode(mode) {
    this.mode = mode; this.controls.enabled = mode === 'orbit';
    if (mode === 'orbit') { this.camera.up.set(0, 0, -1); this.controls.target.copy(this.position); this.controls.update(); }
  }
  recenter() {
    const scale = Math.max(1, this.bodyScale / 3);
    this.camera.position.copy(this.position).add(new THREE.Vector3(-8.5, -10.5, -7.5).multiplyScalar(scale));
    this.camera.up.set(0, 0, -1); this.controls.target.copy(this.position); this.controls.update(); this.chaseDistance = 13 * scale;
  }
  setState(state) {
    if (state.dimensions?.join(',') !== this.dimensionKey) this.buildHull(state.dimensions || [3, 2, 1]);
    this.previousPosition.copy(this.position); this.position.fromArray(state.r_I);
    this.orientation.set(state.q[1], state.q[2], state.q[3], state.q[0]);
    this.body.position.copy(this.position); this.body.quaternion.copy(this.orientation);
    if (this.mode === 'orbit') {
      const displacement = this.position.clone().sub(this.previousPosition);
      this.camera.position.add(displacement); this.controls.target.add(displacement);
    }
    if (this.previousTime !== null && state.sim_time < this.previousTime) this.clearTrajectory();
    if (!this.trajectoryPoints.length || this.position.distanceTo(this.trajectoryPoints.at(-1)) > 0.015) {
      this.trajectoryPoints.push(this.position.clone()); if (this.trajectoryPoints.length > this.maxTrajectory) this.trajectoryPoints.shift();
      const attribute = this.trajectoryGeometry.getAttribute('position');
      this.trajectoryPoints.forEach((p, i) => attribute.setXYZ(i, p.x, p.y, p.z));
      attribute.needsUpdate = true; this.trajectoryGeometry.setDrawRange(0, this.trajectoryPoints.length);
    }
    this.previousTime = state.sim_time;
    const velocity = vector(state.v_I), speed = velocity.length(); this.velocityArrow.visible = speed > 1e-6;
    if (speed > 1e-6) { this.velocityArrow.position.copy(this.position); this.velocityArrow.setDirection(velocity.normalize()); this.velocityArrow.setLength(Math.min(8, 1.5 + speed), 0.25, 0.12); }
    const momentum = vector(state.diagnostics?.angular_momentum_I), h = momentum.length(); this.momentumArrow.visible = h > 1e-8 && this.showMomentum !== false;
    if (h > 1e-8) { this.momentumArrow.position.copy(this.position); this.momentumArrow.setDirection(momentum.normalize()); this.momentumArrow.setLength(this.bodyScale * 1.6, 0.23, 0.1); }
    this.setWrench(state.control?.force_B || [0, 0, 0], state.control?.torque_B || [0, 0, 0]);
  }
  clearTrajectory() { this.trajectoryPoints = []; this.trajectoryGeometry.setDrawRange(0, 0); }
  setWrench(force, torque) {
    const f = vector(force), t = vector(torque), fn = f.length(), tn = t.length();
    this.forceArrow.visible = fn > 0; this.torqueArrow.visible = tn > 0;
    if (fn > 0) { this.forceArrow.setDirection(f.normalize()); this.forceArrow.setLength(this.bodyScale * 0.8 + Math.min(fn / 60, 3), 0.3, 0.16); }
    if (tn > 0) { this.torqueArrow.setDirection(t.normalize()); this.torqueArrow.setLength(this.bodyScale * 1.25, 0.24, 0.15); }
    for (const { mesh, axis, sign } of this.plumeObjects) {
      mesh.visible = force[axis] * sign < 0;
      if (mesh.visible) { const length = 0.7 + Math.min(Math.abs(force[axis]) / 100, 1.3); mesh.scale.set(1, length, 1); mesh.position.copy(UNIT_AXES[axis]).multiplyScalar(sign * (this.dimensions[axis] / 2 + 0.4 * length)); }
    }
  }
  render(dt) {
    if (this.mode === 'chase') {
      const offset = new THREE.Vector3(-this.chaseDistance, 0, -this.chaseDistance * 0.32).applyQuaternion(this.orientation);
      const desired = this.position.clone().add(offset); const alpha = 1 - Math.exp(-dt * 5);
      this.camera.position.lerp(desired, alpha);
      this.camera.up.lerp(new THREE.Vector3(0, 0, -1).applyQuaternion(this.orientation), alpha).normalize();
      this.camera.lookAt(this.position.clone().add(new THREE.Vector3(this.bodyScale * 0.3, 0, 0).applyQuaternion(this.orientation)));
    } else this.controls.update();
    this.stars.position.copy(this.camera.position); this.renderer.render(this.scene, this.camera);
  }
}
