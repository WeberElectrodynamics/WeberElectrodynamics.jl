"""In-process Julia bridge backed by :mod:`juliacall`.

Wraps the :func:`init` / :func:`step!` streaming interface exposed by
``WeberElectrodynamics.jl`` so the Python UI can drive the Julia integrator one
macro-step at a time. No physics is implemented here — everything computational
delegates to the compiled Julia backend.
"""

from __future__ import annotations

import atexit
import os
from typing import Sequence

import numpy as np

from weber_viewer._validation import (
    as_flat_float_array,
    as_int,
    as_positive_float,
    validate_tspan,
)

__all__ = ["JuliaBridge"]


def _jl():
    # Lazy import so importing this module doesn't boot the Julia runtime.
    if os.environ.get("WEBER_VIEWER_ENABLE_JULIACALL_ATEXIT") == "1":
        from juliacall import Main as Main_

        return Main_

    original_register = atexit.register

    def register_without_julia_exit(func, *args, **kwargs):
        # JuliaCall 0.9.31 can bus-error during jl_atexit_hook on Python 3.13
        # after a successful run. Let the process tear down Julia instead.
        if (
            getattr(func, "__module__", "") == "juliacall"
            and getattr(func, "__name__", "") == "at_jl_exit"
        ):
            return func
        return original_register(func, *args, **kwargs)

    atexit.register = register_without_julia_exit
    try:
        from juliacall import Main as Main_
    finally:
        atexit.register = original_register

    return Main_


class JuliaBridge:
    """Thin handle over a Julia ``HamiltonianIntegrator``.

    One instance owns one integrator. Construct via :meth:`from_problem`, then
    drive the simulation with :meth:`step` and :meth:`snapshot`.
    """

    def __init__(self):
        jl = _jl()
        # WeberElectrodynamics re-exports step!, init, and SymmetricProjectionIntegrator,
        # so a single `using` brings everything we need into Main.
        jl.seval("using WeberElectrodynamics")
        self._jl = jl
        self._step_fn = jl.seval("step!")
        self._init_fn = jl.seval("init")
        self._compute_energy = jl.seval("WeberElectrodynamics.compute_total_energy")
        self.prob = None
        self.integrator = None
        self.n: int = 0
        self.dims: int = 0
        self.masses: np.ndarray = np.empty(0)
        self.charges: np.ndarray = np.empty(0)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_problem(
        cls,
        *,
        n_particles: int,
        dims: int,
        q_initial: Sequence[float],
        p_initial: Sequence[float],
        masses: Sequence[float],
        charges: Sequence[float],
        c: float,
        dt: float,
        tspan: tuple[float, float],
        convergence_tolerance: float = 1e-13,
        maximum_iterations: int = 100,
    ) -> "JuliaBridge":
        """Build a Julia ``HamiltonianProblem`` and initialize its integrator."""
        n_particles = as_int("n_particles", n_particles, min_value=1)
        dims = as_int("dims", dims, min_value=1, max_value=3)
        dt = as_positive_float("dt", dt)
        c = as_positive_float("c", c)
        tspan = validate_tspan(tspan)
        convergence_tolerance = as_positive_float(
            "convergence_tolerance", convergence_tolerance,
        )
        maximum_iterations = as_int(
            "maximum_iterations", maximum_iterations, min_value=1,
        )
        expected = n_particles * dims
        q0 = as_flat_float_array("q_initial", q_initial, length=expected)
        p0 = as_flat_float_array("p_initial", p_initial, length=expected)
        m = as_flat_float_array("masses", masses, length=n_particles)
        qch = as_flat_float_array("charges", charges, length=n_particles)
        if np.any(m <= 0.0):
            raise ValueError("masses must contain only positive values")

        self = cls()
        jl = self._jl

        system = jl.HamiltonianSystem(n_particles, dims)
        prob = jl.HamiltonianProblem(
            system,
            tspan,
            q0,
            p0,
            masses=m,
            charges=qch,
            c=c,
            dt=dt,
            convergence_tolerance=convergence_tolerance,
            maximum_iterations=maximum_iterations,
        )
        self.prob = prob
        self.integrator = self._init_fn(prob, jl.SymmetricProjectionIntegrator())

        self.n = n_particles
        self.dims = dims
        self.masses = m.copy()
        self.charges = qch.copy()
        self._tspan = (float(tspan[0]), float(tspan[1]))
        return self

    # ------------------------------------------------------------------
    # Runtime
    # ------------------------------------------------------------------

    def step(self) -> bool:
        """Advance one macro-step. Returns ``True`` while more steps remain."""
        return bool(self._step_fn(self.integrator))

    def snapshot(self) -> tuple[float, np.ndarray, np.ndarray]:
        """Return ``(t, q, p)`` for the current integrator state (copies)."""
        t = float(self.integrator.t)
        q = np.asarray(self.integrator.q, dtype=np.float64).copy()
        p = np.asarray(self.integrator.p, dtype=np.float64).copy()
        return t, q, p

    def energy(self, q: np.ndarray | None = None, p: np.ndarray | None = None) -> float:
        """Total Weber energy at the given state (defaults to current)."""
        if (q is None) != (p is None):
            raise ValueError("q and p must be provided together")
        if q is None:
            _, q, p = self.snapshot()
        else:
            expected = self.n * self.dims
            q = as_flat_float_array("q", q, length=expected)
            p = as_flat_float_array("p", p, length=expected)
        return float(self._compute_energy(q, p, self.prob))

    def recycle(self) -> None:
        """Rebuild the integrator from scratch to continue streaming past ``t_end``."""
        jl = self._jl
        self.integrator = self._init_fn(self.prob, jl.SymmetricProjectionIntegrator())
