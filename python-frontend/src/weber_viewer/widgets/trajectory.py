"""Trajectory panel: 3D GL scene for ``dims == 3``, 2D PlotWidget for ``dims == 2``."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtWidgets import QWidget

from weber_viewer.colors import particle_color_255
from weber_viewer.state import AnimationState

__all__ = ["TrajectoryWidget"]


def TrajectoryWidget(state: AnimationState, parent: QWidget | None = None) -> QWidget:
    """Factory returning the right widget for the problem's dimensionality."""
    if state.bridge.dims == 3:
        return _Trajectory3D(state, parent)
    return _Trajectory2D(state, parent)


class _Trajectory2D(pg.PlotWidget):
    def __init__(self, state: AnimationState, parent: QWidget | None = None):
        super().__init__(parent=parent)
        self.state = state
        self.setBackground("w")
        self.setLabel("bottom", "x")
        self.setLabel("left", "y")
        self.setAspectLocked(True)
        self.showGrid(x=True, y=True, alpha=0.2)
        self.addLegend()

        n = state.bridge.n
        self._trails: list[pg.PlotDataItem] = []
        self._markers: list[pg.ScatterPlotItem] = []
        for i in range(1, n + 1):
            color = particle_color_255(i)
            pen = pg.mkPen(color=color, width=2)
            trail = self.plot([], [], pen=pen, name=f"P{i}")
            marker = pg.ScatterPlotItem(
                size=12, pen=pg.mkPen("k", width=1), brush=pg.mkBrush(color),
            )
            self.addItem(marker)
            self._trails.append(trail)
            self._markers.append(marker)

        state.buffer_updated.connect(self.refresh)
        self.refresh()

    def refresh(self) -> None:
        buf = self.state.buffer
        if buf.count == 0:
            return
        tail = self.state.tail_length
        pos = buf.linearize_last_axis(buf.positions, tail)   # (n, 2, k)
        for i, (trail, marker) in enumerate(zip(self._trails, self._markers)):
            x = pos[i, 0, :]
            y = pos[i, 1, :]
            trail.setData(x, y)
            marker.setData([x[-1]], [y[-1]])


class _Trajectory3D(QWidget):
    def __init__(self, state: AnimationState, parent: QWidget | None = None):
        super().__init__(parent=parent)
        from pyqtgraph.opengl import (   # imported lazily so 2D-only envs don't need OpenGL
            GLAxisItem,
            GLGridItem,
            GLLinePlotItem,
            GLScatterPlotItem,
            GLViewWidget,
        )
        from PyQt6.QtWidgets import QVBoxLayout

        self.state = state
        self._view = GLViewWidget(parent=self)
        self._view.setCameraPosition(distance=10)
        self._view.addItem(GLAxisItem())
        grid = GLGridItem()
        grid.setSize(10, 10)
        self._view.addItem(grid)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._view)

        n = state.bridge.n
        self._trails: list[GLLinePlotItem] = []
        self._markers: list[GLScatterPlotItem] = []
        for i in range(1, n + 1):
            r, g, b = particle_color_255(i)
            rgba = (r / 255, g / 255, b / 255, 1.0)
            trail = GLLinePlotItem(
                pos=np.zeros((1, 3)), color=rgba, width=2.0, antialias=True,
            )
            marker = GLScatterPlotItem(
                pos=np.zeros((1, 3)), color=rgba, size=10,
            )
            self._view.addItem(trail)
            self._view.addItem(marker)
            self._trails.append(trail)
            self._markers.append(marker)

        state.buffer_updated.connect(self.refresh)
        self.refresh()

    def refresh(self) -> None:
        buf = self.state.buffer
        if buf.count == 0:
            return
        tail = self.state.tail_length
        pos = buf.linearize_last_axis(buf.positions, tail)   # (n, 3, k)
        for i, (trail, marker) in enumerate(zip(self._trails, self._markers)):
            pts = pos[i].T                                    # (k, 3)
            trail.setData(pos=pts)
            marker.setData(pos=pts[-1:])
