"""Real HTTP integration: clock semantics, controls, coasting and rejected input."""
import json
import math
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from simulator_server import FIXED_DT, Simulation, make_server
from rigid_body import ConvergenceError
from rigid_body import step as physical_step


class SimulationTest(unittest.TestCase):
    def test_clock_pause_step_and_atomic_errors(self):
        sim = Simulation()
        snap = sim.advance({'wall_dt': .1})
        self.assertAlmostEqual(snap['sim_time'], .1)
        paused = sim.advance({'wall_dt': .25, 'paused': True})
        self.assertEqual(paused['sim_time'], snap['sim_time'])
        single = sim.advance({'step': True, 'wall_dt': .25})
        self.assertAlmostEqual(single['sim_time'], snap['sim_time'] + FIXED_DT)
        self.assertTrue(single['paused'])
        for bad in ({'wall_dt': 1.}, {'speed': -1}, {'force_B': [1, 2]},
                    {'paused': 'false'}, {'torque_B': [math.inf, 0, 0]}):
            before = sim.snapshot()
            with self.assertRaises(ValueError): sim.advance(bad)
            self.assertEqual(before, sim.snapshot())
        # Small browser intervals accumulate into fixed physical steps.
        sim.reset({'paused': False})
        for _ in range(10): sim.advance({'wall_dt': .001})
        self.assertAlmostEqual(sim.sim_time, 2 * FIXED_DT)
        self.assertLess(sim.accumulator, FIXED_DT)

    def test_failed_physical_batch_rolls_back_and_pause_drops_backlog(self):
        sim = Simulation()
        sim.advance({'wall_dt': .001})
        before, accumulator = sim.snapshot(), sim.accumulator
        calls = 0
        def fail_after_one(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise ConvergenceError('controlled nonlinear failure')
            return physical_step(*args, **kwargs)
        with patch('simulator_server.step', side_effect=fail_after_one):
            with self.assertRaises(ConvergenceError):
                sim.advance({'wall_dt': .1, 'force_B': [3, 0, 0]})
        self.assertEqual(calls, 2)
        self.assertEqual(before, sim.snapshot())
        self.assertEqual(accumulator, sim.accumulator)
        sim.advance({'wall_dt': .25, 'paused': True})
        resumed = sim.advance({'wall_dt': 0, 'paused': False})
        self.assertEqual(resumed['sim_time'], before['sim_time'])
        self.assertEqual(sim.accumulator, 0.)

    def test_thruster_release_preserves_new_coast_invariants(self):
        sim = Simulation()
        sim.reset({'preset': 'thruster-sandbox'})
        for _ in range(4):
            active = sim.advance({'wall_dt': .25, 'force_B': [25, 0, 0], 'torque_B': [0, 10, 0]})
        self.assertFalse(active['diagnostics']['free_segment'])
        self.assertGreater(math.hypot(*active['v_I']), .2)
        coast = sim.advance({'wall_dt': .25})
        self.assertTrue(coast['diagnostics']['free_segment'])
        self.assertEqual(coast['diagnostics']['reference'], 'coast segment')
        velocity = coast['v_I']
        for _ in range(20): coast = sim.advance({'wall_dt': .25})
        self.assertEqual(coast['v_I'], velocity)
        for key in ('energy_relative_drift', 'angular_momentum_relative_drift', 'linear_momentum_relative_drift'):
            self.assertLess(coast['diagnostics'][key], 1e-10)
        self.assertLess(math.dist(active['diagnostics']['angular_momentum_I'], coast['diagnostics']['angular_momentum_I']), 1e-9)

    def test_reset_presets_and_rejected_model(self):
        sim = Simulation()
        for preset in ('free-tumble', 'principal-spin', 'intermediate-axis-flip', 'thruster-sandbox'):
            self.assertEqual(sim.reset({'preset': preset})['preset'], preset)
        for bad in ({'preset': 'absent'}, {'mass': 0}, {'dimensions': [1, 0, 1]},
                    {'q': [1, 1, 0, 0]}, {'omega_B': [True, 0, 0]}):
            before = sim.snapshot()
            with self.assertRaises(ValueError): sim.reset(bad)
            self.assertEqual(before, sim.snapshot())


class HTTPTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        Path(cls.temp.name, 'index.html').write_text('<!doctype html><title>fixture</title>')
        cls.server = make_server(port=0, static_dir=cls.temp.name)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f'http://127.0.0.1:{cls.server.server_port}'

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()
        cls.temp.cleanup()

    def request(self, path, data=None, **headers):
        payload = json.dumps(data).encode() if data is not None else None
        request = Request(self.base + path, data=payload,
                          headers={'Content-Type': 'application/json', **headers})
        try:
            with urlopen(request, timeout=10) as response:
                return response.status, response.read(), response.headers
        except HTTPError as response:
            return response.code, response.read(), response.headers

    def test_real_http_roundtrip_and_static_assets(self):
        status, body, _ = self.request('/')
        self.assertEqual(status, 200)
        self.assertIn(b'fixture', body)
        status, body, headers = self.request('/api/reset', {'preset': 'principal-spin', 'paused': True})
        self.assertEqual(status, 200)
        self.assertEqual(headers['Cache-Control'], 'no-store')
        initial = json.loads(body)
        status, body, _ = self.request('/api/advance', {'step': True})
        self.assertEqual(status, 200)
        updated = json.loads(body)
        self.assertAlmostEqual(updated['sim_time'], FIXED_DT)
        self.assertNotEqual(initial['q'], updated['q'])
        status, body, _ = self.request('/api/state')
        self.assertEqual(status, 200)
        self.assertEqual(json.loads(body), updated)
        self.assertEqual(self.request('/api/absent')[0], 404)

    def test_bad_requests_do_not_mutate_state(self):
        _, body, _ = self.request('/api/state')
        for payload in ([], {'mass': -2}, {'mass': float('nan')}, {'q': [0, 0, 0, 0]}):
            self.assertEqual(self.request('/api/reset', payload)[0], 400)
        self.assertEqual(self.request('/api/reset', {}, Origin='https://example.invalid')[0], 403)
        self.assertEqual(self.request('/api/reset', {}, **{'Content-Type': 'text/plain'})[0], 415)
        self.assertEqual(self.request('/api/state')[1], body)


if __name__ == '__main__':
    unittest.main()
