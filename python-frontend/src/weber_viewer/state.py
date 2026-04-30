"""Qt-aware animation state, analogous to ``AnimationState`` in the Makie extension."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from PyQt6.QtCore import QObject, QTimer, pyqtSignal

from weber_viewer._validation import as_int
from weber_viewer.bridge import JuliaBridge
from weber_viewer.buffer import RollingBuffer

__all__ = ["AnimationState", "PhaseSelection"]


@dataclass
class PhaseSelection:
    """What the phase-space panel is currently displaying."""

    mode: str               # "pair" or "particle"
    pair: tuple[int, int]   # 1-indexed
    particle: int           # 1-indexed
    component: int          # 1-indexed (x=1, y=2, z=3)

    def label(self) -> str:
        if self.mode == "pair":
            return f"Pair ({self.pair[0]},{self.pair[1]})"
        return f"Particle {self.particle}"


class AnimationState(QObject):
    """Owns the bridge + rolling buffer + UI-observable state.

    Drives a 60 Hz ``QTimer``; on each tick advances the integrator
    ``compute_batch`` times and emits :attr:`buffer_updated` for widgets
    to refresh.
    """

    buffer_updated = pyqtSignal()
    phase_changed = pyqtSignal()
    play_state_changed = pyqtSignal(bool)

    FRAME_INTERVAL_MS = 1000 // 60

    def __init__(
        self,
        bridge: JuliaBridge,
        buffer: RollingBuffer,
        *,
        tail_length: int,
        compute_batch: int,
        phase: PhaseSelection,
        parent: QObject | None = None,
    ):
        super().__init__(parent)
        tail_length = as_int("tail_length", tail_length, min_value=1)
        compute_batch = as_int("compute_batch", compute_batch, min_value=1)
        if tail_length > buffer.capacity:
            raise ValueError("tail_length must be <= buffer capacity")

        self.bridge = bridge
        self.buffer = buffer
        self.masses = bridge.masses

        self._tail_length = int(tail_length)
        self._compute_batch = int(compute_batch)
        self._is_playing = False
        self._phase = phase
        self._total_steps = 0

        # Seed with the initial state so widgets have something to draw.
        t0, q0, p0 = bridge.snapshot()
        self.E0 = bridge.energy(q0, p0)
        self.buffer.push_step(t0, q0, p0, self.masses, self.E0)

        self._timer = QTimer(self)
        self._timer.setInterval(self.FRAME_INTERVAL_MS)
        self._timer.timeout.connect(self._tick)

    # ------------------------------------------------------------------
    # Observable-style properties
    # ------------------------------------------------------------------

    @property
    def is_playing(self) -> bool:
        return self._is_playing

    def set_playing(self, playing: bool) -> None:
        if playing == self._is_playing:
            return
        self._is_playing = playing
        self.play_state_changed.emit(playing)

    def toggle_playing(self) -> None:
        self.set_playing(not self._is_playing)

    @property
    def tail_length(self) -> int:
        return self._tail_length

    def set_tail_length(self, n: int) -> None:
        self._tail_length = max(1, min(int(n), self.buffer.capacity))
        self.buffer_updated.emit()

    @property
    def compute_batch(self) -> int:
        return self._compute_batch

    def set_compute_batch(self, n: int) -> None:
        self._compute_batch = max(1, int(n))

    @property
    def phase(self) -> PhaseSelection:
        return self._phase

    def set_phase(self, sel: PhaseSelection) -> None:
        self._phase = sel
        self.phase_changed.emit()
        self.buffer_updated.emit()

    @property
    def total_steps(self) -> int:
        return self._total_steps

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self) -> None:
        """Start the 60 Hz render timer (does not auto-play)."""
        self._timer.start()

    def stop(self) -> None:
        self._timer.stop()

    def reset(self) -> None:
        """Rebuild the integrator and re-seed the buffer from t=0."""
        self.bridge.recycle()
        self.buffer.reset()
        t0, q0, p0 = self.bridge.snapshot()
        self.buffer.push_step(t0, q0, p0, self.masses, self.E0)
        self._total_steps = 0
        self.buffer_updated.emit()

    # ------------------------------------------------------------------
    # Tick loop
    # ------------------------------------------------------------------

    def _tick(self) -> None:
        if not self._is_playing:
            return
        bridge = self.bridge
        buf = self.buffer
        masses = self.masses
        for _ in range(self._compute_batch):
            more = bridge.step()
            if not more:
                # Loop the integrator so streaming continues indefinitely, matching
                # _recycle_integrator! in the Makie extension.
                bridge.recycle()
            t, q, p = bridge.snapshot()
            E = bridge.energy(q, p)
            buf.push_step(t, q, p, masses, E)
            self._total_steps += 1
        self.buffer_updated.emit()

    # ------------------------------------------------------------------
    # Computed read helpers (keep widgets out of buffer internals)
    # ------------------------------------------------------------------

    def energy_error_text(self) -> str:
        if self.buffer.count < 2:
            return "ΔE/E₀ (max) = --"
        E = self.buffer.linearize_1d(self.buffer.total_energy)
        dE = np.abs(E - self.E0)
        if abs(self.E0) < 1e-30:
            return f"|ΔE| (max) = {dE.max():.3e}"
        return f"ΔE/E₀ (max) = {dE.max() / abs(self.E0) * 100:.3e} %"

    def time_text(self) -> str:
        if self.buffer.count == 0:
            return "t = --"
        t = float(self.buffer.latest(self.buffer.t))
        return f"t = {t:.4f}"

    def step_text(self) -> str:
        return f"step = {self._total_steps}"
