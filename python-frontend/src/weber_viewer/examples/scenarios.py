"""Reusable example configurations for visual and headless demos."""

from __future__ import annotations

from math import cos, pi, sin
from typing import Any

PROBLEM_KEYS = (
    "n_particles",
    "dims",
    "q_initial",
    "p_initial",
    "masses",
    "charges",
    "c",
    "dt",
    "tspan",
)


def problem_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    """Extract the kwargs needed by ``JuliaBridge.from_problem``."""
    return {key: config[key] for key in PROBLEM_KEYS}


def two_body_config() -> dict[str, Any]:
    """Two opposite charges in a near-circular 3D orbit."""
    return {
        "n_particles": 2,
        "dims": 3,
        "masses": [1.0, 1.0],
        "charges": [1.0, -1.0],
        "q_initial": [1.0, 0.0, 0.15, -1.0, 0.0, -0.15],
        "p_initial": [0.0, 0.5, 0.0, 0.0, -0.5, 0.0],
        "c": 100.0,
        "dt": 0.01,
        "tspan": (0.0, 200.0),
        "buffer_size": 2000,
        "tail_length": 400,
        "compute_batch": 2,
        "initial_pair": (1, 2),
        "phase_mode": "pair",
        "figure_size": (1400, 900),
    }


def three_body_polygon_config() -> dict[str, Any]:
    """Compact planar alternating-charge triangle with tangential momenta."""
    n = 3
    radius = 1.35
    speed = 0.32

    q: list[float] = []
    p: list[float] = []
    for i in range(n):
        theta = 2.0 * pi * i / n
        q.extend([radius * cos(theta), radius * sin(theta)])
        p.extend([speed * -sin(theta), speed * cos(theta)])

    return {
        "n_particles": n,
        "dims": 2,
        "masses": [1.0, 1.0, 1.0],
        "charges": [1.0, -1.0, 1.0],
        "q_initial": q,
        "p_initial": p,
        "c": 120.0,
        "dt": 0.005,
        "tspan": (0.0, 120.0),
        "buffer_size": 2500,
        "tail_length": 800,
        "compute_batch": 3,
        "initial_pair": (1, 2),
        "phase_mode": "pair",
        "figure_size": (1400, 900),
    }
