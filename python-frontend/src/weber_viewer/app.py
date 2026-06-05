"""MainWindow + public ``animate_weber`` entry point."""

from __future__ import annotations

import sys
from typing import Sequence

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from weber_viewer._validation import as_int, validate_tspan
from weber_viewer.bridge import JuliaBridge
from weber_viewer.buffer import RollingBuffer
from weber_viewer.state import AnimationState, PhaseSelection
from weber_viewer.widgets import (
    ControlsWidget,
    HUDWidget,
    PhaseSpaceWidget,
    TrajectoryWidget,
)

__all__ = ["animate_weber", "MainWindow"]


class MainWindow(QMainWindow):
    """Dashboard: trajectory panel + phase-space sidebar + HUD + controls."""

    def __init__(
        self,
        state: AnimationState,
        *,
        figure_size: tuple[int, int] = (1400, 900),
    ):
        super().__init__()
        self.state = state
        self.setWindowTitle("Weber Electrodynamics Viewer")
        self.resize(*figure_size)

        # Reparent so Qt's object tree keeps the state alive for the window's
        # lifetime — avoids PyQt sip deleting it between ticks.
        state.setParent(self)

        self.trajectory = TrajectoryWidget(state)
        self.phase_space = PhaseSpaceWidget(state)
        self.hud = HUDWidget(state)
        self.controls = ControlsWidget(state)

        # Top row: trajectory | phase-space | HUD (splitter so the user can rearrange).
        top = QSplitter(Qt.Orientation.Horizontal)
        top.addWidget(self.trajectory)
        top.addWidget(self.phase_space)
        top.addWidget(self.hud)
        top.setStretchFactor(0, 3)
        top.setStretchFactor(1, 2)
        top.setStretchFactor(2, 1)

        central = QWidget(self)
        outer = QVBoxLayout(central)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.addWidget(top, 1)
        outer.addWidget(self.controls)
        self.setCentralWidget(central)


def animate_weber(
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
    buffer_size: int = 2000,
    tail_length: int = 200,
    compute_batch: int = 1,
    initial_pair: tuple[int, int] = (1, 2),
    phase_mode: str = "pair",
    initial_particle: int = 1,
    initial_component: int = 1,
    figure_size: tuple[int, int] = (1400, 900),
    autoplay: bool = True,
    exec_app: bool = True,
) -> MainWindow:
    """Launch the streaming animation viewer against a fresh Julia integrator.

    Mirrors the Julia ``animate_weber(prob::HamiltonianProblem; ...)`` kwarg surface.
    Returns the :class:`MainWindow` (caller keeps a reference so Qt doesn't
    garbage-collect it mid-event-loop).
    """
    n_particles = as_int("n_particles", n_particles, min_value=2)
    dims = as_int("dims", dims, min_value=2, max_value=3)
    buffer_size = as_int("buffer_size", buffer_size, min_value=1)
    tail_length = as_int("tail_length", tail_length, min_value=1)
    compute_batch = as_int("compute_batch", compute_batch, min_value=1)
    initial_particle = as_int(
        "initial_particle", initial_particle, min_value=1, max_value=n_particles,
    )
    initial_component = as_int(
        "initial_component", initial_component, min_value=1, max_value=dims,
    )
    tspan = validate_tspan(tspan)

    if dims < 2:
        raise ValueError(f"Animation viewer requires dims >= 2 (got {dims})")
    if tail_length > buffer_size:
        raise ValueError("tail_length must be <= buffer_size")
    if phase_mode not in ("pair", "particle"):
        raise ValueError("phase_mode must be 'pair' or 'particle'")
    if len(initial_pair) != 2:
        raise ValueError("initial_pair must contain exactly two particle indices")
    pair_i = as_int("initial_pair[0]", initial_pair[0], min_value=1, max_value=n_particles)
    pair_j = as_int("initial_pair[1]", initial_pair[1], min_value=1, max_value=n_particles)
    if pair_i == pair_j:
        raise ValueError("initial_pair indices must be distinct")
    initial_pair = tuple(sorted((pair_i, pair_j)))

    app = QApplication.instance() or QApplication(sys.argv)
    # Keep a module-level reference so the QApplication survives across
    # ``animate_weber`` calls (PyQt deletes it when the last Python ref drops).
    globals()["_qapp"] = app

    bridge = JuliaBridge.from_problem(
        n_particles=n_particles,
        dims=dims,
        q_initial=q_initial,
        p_initial=p_initial,
        masses=masses,
        charges=charges,
        c=c,
        dt=dt,
        tspan=tspan,
    )
    buffer = RollingBuffer(buffer_size, n_particles, dims)
    phase = PhaseSelection(
        mode=phase_mode,
        pair=initial_pair,
        particle=initial_particle,
        component=initial_component,
    )
    state = AnimationState(
        bridge, buffer,
        tail_length=tail_length,
        compute_batch=compute_batch,
        phase=phase,
    )

    window = MainWindow(state, figure_size=figure_size)
    window.show()
    state.start()
    if autoplay:
        state.set_playing(True)

    if exec_app:
        app.exec()
    return window
