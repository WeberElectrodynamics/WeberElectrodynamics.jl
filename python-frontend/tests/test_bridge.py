"""Integration test: JuliaBridge drives the Julia integrator forward.

Requires a working juliacall + WeberElectrodynamics dev environment.
Run with ``PYTHON_JULIAPKG_PROJECT=<repo root>``.
"""

import importlib.util

import numpy as np
import pytest

if importlib.util.find_spec("juliacall") is None:
    pytest.skip("juliacall is not installed", allow_module_level=True)


@pytest.fixture(scope="module")
def bridge():
    from weber_viewer.bridge import JuliaBridge

    return JuliaBridge.from_problem(
        n_particles=2,
        dims=3,
        q_initial=[1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
        p_initial=[0.0, 0.5, 0.0, 0.0, -0.5, 0.0],
        masses=[1.0, 1.0],
        charges=[1.0, -1.0],
        c=100.0,
        dt=0.01,
        tspan=(0.0, 1.0),
    )


def test_snapshot_initial(bridge):
    t, q, p = bridge.snapshot()
    assert t == 0.0
    assert q.shape == (6,)
    assert np.allclose(q[:3], [1.0, 0.0, 0.0])


def test_step_advances_time(bridge):
    t0, _, _ = bridge.snapshot()
    assert bridge.step() is True
    t1, _, _ = bridge.snapshot()
    assert t1 > t0


def test_energy_conserved_to_order_of_dt(bridge):
    # Symplectic integrator: energy oscillates but stays bounded.
    _, q0, p0 = bridge.snapshot()
    E0 = bridge.energy(q0, p0)
    for _ in range(50):
        bridge.step()
    _, q, p = bridge.snapshot()
    E = bridge.energy(q, p)
    # Loose bound: relative drift < 1% over 50 macro-steps.
    assert abs((E - E0) / E0) < 1e-2
