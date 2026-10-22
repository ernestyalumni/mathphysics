"""Local authoritative 6-DOF sandbox. HTTP is transport, never the integrator clock."""
from __future__ import annotations
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import math
from pathlib import Path
import threading
from urllib.parse import urlsplit

from rigid_body import RigidBodyModel, RigidBodyState, step, diagnostics, ConvergenceError

FIXED_DT = 1 / 240
PRESETS = {
    'free-tumble': ((.12, .32, .5), (.12, .05, 0.)),
    'principal-spin': ((0., 0., .8), (0., 0., 0.)),
    'intermediate-axis-flip': ((.02, 1., .01), (0., 0., 0.)),
    'thruster-sandbox': ((0., 0., 0.), (0., 0., 0.)),
}
ZERO = (0., 0., 0.)


def number(value, name, lo, hi):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} must be a number')
    value = float(value)
    if not math.isfinite(value) or not lo <= value <= hi:
        raise ValueError(f'{name} must be finite and in [{lo}, {hi}]')
    return value


def vector(value, name, size=3, limit=1e5):
    if not isinstance(value, (list, tuple)) or len(value) != size:
        raise ValueError(f'{name} must contain {size} numbers')
    return tuple(number(v, name, -limit, limit) for v in value)


def boolean(value, name):
    if not isinstance(value, bool):
        raise ValueError(f'{name} must be true or false')
    return value


def relative_error(value, reference):
    if isinstance(value, (list, tuple)):
        return math.dist(value, reference) / max(math.hypot(*reference), 1.)
    return abs(value - reference) / max(abs(reference), 1.)


class Simulation:
    """One local pilot, fixed steps, atomic request updates, bounded catch-up.

    Drift scales use max(initial magnitude, 1 SI unit), so an initially zero
    momentum has a meaningful absolute floor. A new coast segment starts when
    the applied wrench is released; energy added by thrusters is physical.
    """
    def __init__(self):
        self.lock = threading.RLock()
        self.reset({})

    def reset(self, data):
        with self.lock:
            preset = data.get('preset', 'free-tumble')
            if not isinstance(preset, str) or preset not in PRESETS:
                raise ValueError('unknown preset')
            mass = number(data.get('mass', 100.), 'mass', .001, 1e6)
            dimensions = vector(data.get('dimensions', [3., 2., 1.]), 'dimensions', limit=1000)
            if min(dimensions) < .001:
                raise ValueError('dimensions must be at least 0.001 m')
            omega, velocity = PRESETS[preset]
            model = RigidBodyModel.from_box(mass, dimensions)
            state = RigidBodyState(
                vector(data.get('r_I', ZERO), 'r_I', limit=1e9),
                vector(data.get('v_I', velocity), 'v_I', limit=1e6),
                vector(data.get('q', [1., 0., 0., 0.]), 'q', size=4, limit=1.),
                vector(data.get('omega_B', omega), 'omega_B', limit=100.),
            )
            paused = boolean(data.get('paused', False), 'paused')
            speed = number(data.get('speed', 1.), 'speed', .1, 10.)
            self.model, self.state = model, state
            self.dimensions, self.preset = dimensions, preset
            self.paused, self.speed = paused, speed
            self.sim_time, self.accumulator = 0., 0.
            self.force, self.torque = ZERO, ZERO
            self.reference = diagnostics(model, state)
            self.reference_label = 'initial'
            self.free_segment = True
            return self.snapshot()

    def advance(self, data):
        with self.lock:
            elapsed = number(data.get('wall_dt', 0.), 'wall_dt', 0., .25)
            speed = number(data.get('speed', self.speed), 'speed', .1, 10.)
            paused = boolean(data.get('paused', self.paused), 'paused')
            single = boolean(data.get('step', False), 'step')
            force = vector(data.get('force_B', ZERO), 'force_B')
            torque = vector(data.get('torque_B', ZERO), 'torque_B')
            active = any(force) or any(torque)
            accumulator = self.accumulator
            if single:
                steps, accumulator, paused = 1, 0., True
            elif paused:
                steps, accumulator = 0, 0.
            else:
                accumulator += elapsed * speed
                steps = min(int((accumulator + 1e-12) / FIXED_DT), 600)
                accumulator = max(0., accumulator - steps * FIXED_DT)
            state = self.state
            reference, label = self.reference, self.reference_label
            # The transition reference belongs to the state before its first
            # zero-wrench step. While paused, the previous wrench remains history.
            free = self.free_segment
            if steps:
                if not active and not free:
                    reference, label = diagnostics(self.model, state), 'coast segment'
                free = not active
                for _ in range(steps):
                    state = step(self.model, state, FIXED_DT, force, torque)
            # Commit only after every requested step succeeds.
            self.state, self.accumulator = state, accumulator
            self.sim_time += steps * FIXED_DT
            self.paused, self.speed = paused, speed
            self.reference, self.reference_label = reference, label
            self.free_segment = free
            self.force, self.torque = (force, torque) if steps else (ZERO, ZERO)
            return self.snapshot()

    def snapshot(self):
        with self.lock:
            d = diagnostics(self.model, self.state)
            d.update({
                'angular_momentum_I': d['angular_momentum_inertial'],
                'linear_momentum_I': d['linear_momentum_inertial'],
                'qnorm': d['quaternion_norm'],
                'energy_relative_drift': relative_error(d['energy'], self.reference['energy']),
                'angular_momentum_relative_drift': relative_error(d['angular_momentum_inertial'], self.reference['angular_momentum_inertial']),
                'linear_momentum_relative_drift': relative_error(d['linear_momentum_inertial'], self.reference['linear_momentum_inertial']),
                'reference': self.reference_label, 'free_segment': self.free_segment,
            })
            return {'r_I': self.state.position, 'v_I': self.state.velocity,
                    'q': self.state.quaternion, 'omega_B': self.state.omega_body,
                    'sim_time': self.sim_time, 'mass': self.model.mass,
                    'dimensions': self.dimensions, 'paused': self.paused,
                    'speed': self.speed, 'preset': self.preset, 'fixed_dt': FIXED_DT,
                    'diagnostics': d,
                    'control': {'force_B': self.force, 'torque_B': self.torque}}


class Handler(SimpleHTTPRequestHandler):
    server_version = 'RigidBodyLab/0.1'

    def _json(self, status, value):
        content = json.dumps(value, allow_nan=False, separators=(',', ':')).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(content)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(content)

    def do_GET(self):
        route = urlsplit(self.path).path
        if route == '/api/state':
            self._json(200, self.server.simulation.snapshot())
        elif route == '/api/health':
            self._json(200, {'status': 'ok', 'model': 'ideal rigid body', 'fixed_dt': FIXED_DT})
        elif route.startswith('/api/'):
            self._json(404, {'error': 'unknown API route'})
        else:
            super().do_GET()

    def do_POST(self):
        route = urlsplit(self.path).path
        if route not in ('/api/reset', '/api/advance'):
            self._json(404, {'error': 'unknown API route'})
            return
        try:
            origin = self.headers.get('Origin')
            if origin and origin != 'http://' + self.headers.get('Host', ''):
                self._json(403, {'error': 'same-origin requests required'})
                return
            if self.headers.get_content_type() != 'application/json':
                self._json(415, {'error': 'application/json required'})
                return
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 65536:
                self._json(413, {'error': 'JSON body must be 1..65536 bytes'})
                return
            data = json.loads(self.rfile.read(length),
                              parse_constant=lambda _: (_ for _ in ()).throw(ValueError('nonfinite JSON number')))
            if not isinstance(data, dict):
                raise ValueError('JSON object required')
            method = self.server.simulation.reset if route == '/api/reset' else self.server.simulation.advance
            self._json(200, method(data))
        except (ValueError, TypeError, OverflowError) as exc:
            self._json(400, {'error': str(exc)})
        except ConvergenceError as exc:
            self._json(422, {'error': str(exc)})

    def log_message(self, fmt, *args):
        if self.server.verbose:
            super().log_message(fmt, *args)


def make_server(host='127.0.0.1', port=8765, static_dir=None, verbose=False):
    if host not in ('127.0.0.1', 'localhost'):
        raise ValueError('this single-pilot sandbox binds only to loopback')
    directory = Path(static_dir or '/tmp/mathphysics-rigid-body-viewer-dist').resolve()
    if not (directory / 'index.html').is_file():
        raise ValueError(f'frontend build missing: {directory}/index.html')
    server = ThreadingHTTPServer((host, port), partial(Handler, directory=str(directory)))
    server.simulation, server.verbose = Simulation(), verbose
    server.daemon_threads = True
    return server


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8765)
    parser.add_argument('--static-dir', type=Path)
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()
    with make_server(port=args.port, static_dir=args.static_dir, verbose=args.verbose) as server:
        print(f'Rigid Body Lab: http://127.0.0.1:{server.server_port}', flush=True)
        print('Ideal constant-mass rigid body; no gravity/drag. Ctrl-C stops the server.', flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == '__main__':
    main()
