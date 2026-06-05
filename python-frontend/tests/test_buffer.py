"""Unit tests for RollingBuffer — no Julia runtime needed."""

import numpy as np
import pytest

from weber_viewer.buffer import RollingBuffer


def _push_linear(buf: RollingBuffer, k: int) -> None:
    """Push k frames with deterministic content so we can check ordering."""
    masses = np.ones(buf.n)
    for step in range(k):
        # Use step as the position of particle 1 along x; other coords stay zero.
        q = np.zeros(buf.n * buf.dims)
        q[0] = float(step)
        p = np.zeros(buf.n * buf.dims)
        buf.push_step(t=float(step), q=q, p=p, masses=masses, energy=float(step))


def test_push_within_capacity():
    buf = RollingBuffer(capacity=10, n_particles=2, dims=3)
    _push_linear(buf, 5)

    assert buf.count == 5
    assert buf.cursor == 5
    t = buf.linearize_1d(buf.t)
    assert t.size == 5
    assert np.allclose(t, [0, 1, 2, 3, 4])


def test_push_wraps_and_linearize_chronological():
    buf = RollingBuffer(capacity=4, n_particles=2, dims=3)
    _push_linear(buf, 10)   # wraps 2.5×

    assert buf.count == 4
    # After 10 pushes into capacity 4, the last four frames are 6,7,8,9.
    t = buf.linearize_1d(buf.t)
    assert np.allclose(t, [6, 7, 8, 9])


def test_tail_clipping():
    buf = RollingBuffer(capacity=10, n_particles=2, dims=3)
    _push_linear(buf, 8)
    tail = buf.linearize_1d(buf.t, tail=3)
    assert np.allclose(tail, [5, 6, 7])


def test_pair_separation_computed():
    buf = RollingBuffer(capacity=3, n_particles=2, dims=2)
    masses = np.array([1.0, 1.0])
    # Particle 1 at origin, particle 2 at (3, 4) → separation 5.
    q = np.array([0.0, 0.0, 3.0, 4.0])
    p = np.array([0.0, 1.0, 0.0, -1.0])  # rel velocity along y: -2
    buf.push_step(t=0.0, q=q, p=p, masses=masses, energy=0.0)

    assert buf.pair_separation[(1, 2)][0] == pytest.approx(5.0)
    # r̂·Δv = (3·0 + 4·(-2))/5 = -8/5
    assert buf.pair_radial_velocity[(1, 2)][0] == pytest.approx(-1.6)


def test_latest():
    buf = RollingBuffer(capacity=3, n_particles=1, dims=2)
    _push_linear(buf, 5)
    assert buf.latest(buf.t) == pytest.approx(4.0)


def test_empty_latest_raises():
    buf = RollingBuffer(capacity=3, n_particles=1, dims=2)
    with pytest.raises(IndexError):
        buf.latest(buf.t)


def test_linearize_last_axis_shape():
    buf = RollingBuffer(capacity=5, n_particles=2, dims=3)
    _push_linear(buf, 3)
    out = buf.linearize_last_axis(buf.positions, tail=2)
    assert out.shape == (2, 3, 2)
