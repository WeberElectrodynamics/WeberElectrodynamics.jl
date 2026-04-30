"""Tests for packaged example entry points and scenario shapes."""

from __future__ import annotations

import importlib

import numpy as np

from weber_viewer.examples import (
    PROBLEM_KEYS,
    problem_kwargs,
    three_body_polygon_config,
    two_body_config,
)


def _assert_valid_config(config):
    n = config["n_particles"]
    dims = config["dims"]
    expected = n * dims
    assert len(config["q_initial"]) == expected
    assert len(config["p_initial"]) == expected
    assert len(config["masses"]) == n
    assert len(config["charges"]) == n
    assert config["dt"] > 0
    assert config["tspan"][1] > config["tspan"][0]
    assert config["tail_length"] <= config["buffer_size"]
    assert config["compute_batch"] >= 1
    assert np.isfinite(np.asarray(config["q_initial"], dtype=float)).all()
    assert np.isfinite(np.asarray(config["p_initial"], dtype=float)).all()


def test_example_configs_have_valid_shapes():
    _assert_valid_config(two_body_config())
    _assert_valid_config(three_body_polygon_config())


def test_problem_kwargs_extracts_bridge_surface():
    config = two_body_config()
    kwargs = problem_kwargs(config)
    assert tuple(kwargs) == PROBLEM_KEYS
    assert "buffer_size" not in kwargs
    assert kwargs["n_particles"] == 2


def test_packaged_example_entrypoints_are_importable():
    for module_name in (
        "weber_viewer.examples.two_body_streaming",
        "weber_viewer.examples.three_body_polygon",
        "weber_viewer.examples.smoke_probe",
    ):
        module = importlib.import_module(module_name)
        assert callable(module.main)
