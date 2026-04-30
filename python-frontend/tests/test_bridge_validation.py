"""Fast validation tests that do not need to start Julia."""

from __future__ import annotations

import pytest

from weber_viewer.bridge import JuliaBridge


BASE_PROBLEM = {
    "n_particles": 2,
    "dims": 2,
    "q_initial": [1.0, 0.0, -1.0, 0.0],
    "p_initial": [0.0, 0.5, 0.0, -0.5],
    "masses": [1.0, 1.0],
    "charges": [1.0, -1.0],
    "c": 100.0,
    "dt": 0.01,
    "tspan": (0.0, 1.0),
}


def _problem(**updates):
    config = BASE_PROBLEM.copy()
    config.update(updates)
    return config


def test_rejects_bad_initial_condition_length_before_julia_starts():
    with pytest.raises(ValueError, match="q_initial"):
        JuliaBridge.from_problem(**_problem(q_initial=[1.0, 0.0]))


def test_rejects_nonpositive_dt_before_julia_starts():
    with pytest.raises(ValueError, match="dt"):
        JuliaBridge.from_problem(**_problem(dt=0.0))


def test_rejects_invalid_tspan_before_julia_starts():
    with pytest.raises(ValueError, match="tspan"):
        JuliaBridge.from_problem(**_problem(tspan=(1.0, 1.0)))


def test_rejects_nonpositive_masses_before_julia_starts():
    with pytest.raises(ValueError, match="masses"):
        JuliaBridge.from_problem(**_problem(masses=[1.0, 0.0]))
