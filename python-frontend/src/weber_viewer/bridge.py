"""In-process Julia bridge backed by :mod:`juliacall`.

Wraps the :func:`init` / :func:`step!` streaming interface exposed by
``WeberElectrodynamics.jl`` so the Python UI can drive the Julia integrator one
macro-step at a time. No physics is implemented here — everything computational
delegates to the compiled Julia backend.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

__all__ = ["JuliaBridge"]


def _jl():
    # Lazy import so importing this module doesn't boot the Julia runtime.
    from juliacall import Main as Main_

    return Main_


class JuliaBridge:
    """Thin handle over a Julia ``WeberIntegrator``.

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
        """Build a Julia ``WeberProblem`` and initialize its integrator."""
        self = cls()
        jl = self._jl

        expected = n_particles * dims
        q0 = np.asarray(q_initial, dtype=np.float64)
        p0 = np.asarray(p_initial, dtype=np.float64)
        m = np.asarray(masses, dtype=np.float64)
        qch = np.asarray(charges, dtype=np.float64)
        if q0.size != expected or p0.size != expected:
            raise ValueError(
                f"q_initial/p_initial must have length n_particles*dims = {expected}; "
                f"got q0={q0.size}, p0={p0.size}"
            )
        if m.size != n_particles or qch.size != n_particles:
            raise ValueError(
                f"masses and charges must have length n_particles = {n_particles}"
            )

        system = jl.WeberSystem(n_particles, dims)
        prob = jl.WeberProblem(
            system,
            (float(tspan[0]), float(tspan[1])),
            q0,
            p0,
            masses=m,
            charges=qch,
            c=float(c),
            dt=float(dt),
            convergence_tolerance=float(convergence_tolerance),
            maximum_iterations=int(maximum_iterations),
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
        if q is None or p is None:
            _, q, p = self.snapshot()
        return float(self._compute_energy(q, p, self.prob))

    def recycle(self) -> None:
        """Rebuild the integrator from scratch to continue streaming past ``t_end``."""
        jl = self._jl
        self.integrator = self._init_fn(self.prob, jl.SymmetricProjectionIntegrator())
