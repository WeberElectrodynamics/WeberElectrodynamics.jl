"""Phase-space sidebar: pair (r, ṙ) or particle (qₖ, pₖ) portrait."""

from __future__ import annotations

import pyqtgraph as pg
from PyQt6.QtWidgets import QWidget

from weber_viewer.colors import particle_color_255
from weber_viewer.state import AnimationState

__all__ = ["PhaseSpaceWidget"]


class PhaseSpaceWidget(pg.PlotWidget):
    def __init__(self, state: AnimationState, parent: QWidget | None = None):
        super().__init__(parent=parent)
        self.state = state
        self.setBackground("w")
        self.showGrid(x=True, y=True, alpha=0.2)

        self._curve = self.plot([], [], pen=pg.mkPen("#333333", width=2))
        self._marker = pg.ScatterPlotItem(size=10, brush=pg.mkBrush("#000000"))
        self.addItem(self._marker)

        state.buffer_updated.connect(self.refresh)
        state.phase_changed.connect(self._on_phase_changed)
        self._on_phase_changed()
        self.refresh()

    def _on_phase_changed(self) -> None:
        phase = self.state.phase
        if phase.mode == "pair":
            self.setLabel("bottom", "r")
            self.setLabel("left", "dr/dt")
            # Marker color tracks the first particle of the pair for consistency.
            color = particle_color_255(phase.pair[0])
        else:
            self.setLabel("bottom", f"q{phase.component}")
            self.setLabel("left", f"p{phase.component}")
            color = particle_color_255(phase.particle)
        self._curve.setPen(pg.mkPen(color=color, width=2))
        self._marker.setBrush(pg.mkBrush(color))

    def refresh(self) -> None:
        buf = self.state.buffer
        if buf.count == 0:
            return
        tail = self.state.tail_length
        phase = self.state.phase
        if phase.mode == "pair":
            key = tuple(sorted(phase.pair))
            if key not in buf.pair_separation:
                return
            x = buf.linearize_1d(buf.pair_separation[key], tail)
            y = buf.linearize_1d(buf.pair_radial_velocity[key], tail)
        else:
            i = phase.particle - 1
            d = phase.component - 1
            if not (0 <= i < buf.n and 0 <= d < buf.dims):
                return
            x = buf.linearize_last_axis(buf.particle_q[i, d : d + 1, :], tail)[0]
            y = buf.linearize_last_axis(buf.particle_p[i, d : d + 1, :], tail)[0]

        if x.size == 0:
            return
        self._curve.setData(x, y)
        self._marker.setData([x[-1]], [y[-1]])
