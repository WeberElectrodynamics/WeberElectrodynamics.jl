"""Right-column HUD: energy error, time, step counter."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QLabel, QVBoxLayout, QWidget

from weber_viewer.state import AnimationState

__all__ = ["HUDWidget"]


class HUDWidget(QWidget):
    def __init__(self, state: AnimationState, parent: QWidget | None = None):
        super().__init__(parent=parent)
        self.state = state

        self._energy = _mono_label()
        self._time = _mono_label()
        self._step = _mono_label()

        layout = QVBoxLayout(self)
        layout.addWidget(self._energy)
        layout.addWidget(self._time)
        layout.addWidget(self._step)
        layout.addStretch(1)

        state.buffer_updated.connect(self.refresh)
        self.refresh()

    def refresh(self) -> None:
        self._energy.setText(self.state.energy_error_text())
        self._time.setText(self.state.time_text())
        self._step.setText(self.state.step_text())


def _mono_label() -> QLabel:
    lbl = QLabel("--")
    lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    lbl.setStyleSheet("font-family: Menlo, Consolas, monospace; font-size: 12pt;")
    return lbl
