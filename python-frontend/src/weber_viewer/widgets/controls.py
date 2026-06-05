"""Bottom control strip: play/pause, reset, sliders, phase mode dropdown."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from weber_viewer.state import AnimationState, PhaseSelection

__all__ = ["ControlsWidget"]


def _log_linear_range(max_val: int) -> list[int]:
    """1,2,...,9, 10,20,...,90, 100,200,...,1000 — matches the Makie ext."""
    values: list[int] = []
    decade = 1
    while decade <= max_val:
        for k in range(1, 10):
            v = k * decade
            if v <= max_val:
                values.append(v)
        decade *= 10
    return values


class ControlsWidget(QWidget):
    SPEED_VALUES = _log_linear_range(1000)

    def __init__(self, state: AnimationState, parent: QWidget | None = None):
        super().__init__(parent=parent)
        self.state = state

        self._play_btn = QPushButton("▶ Play")
        self._reset_btn = QPushButton("⟲ Reset")
        self._play_btn.clicked.connect(state.toggle_playing)
        self._reset_btn.clicked.connect(state.reset)
        state.play_state_changed.connect(self._on_play_state)

        self._tail_slider = QSlider(Qt.Orientation.Horizontal)
        self._tail_slider.setRange(1, state.buffer.capacity)
        self._tail_slider.setSingleStep(10)
        self._tail_slider.setValue(state.tail_length)
        self._tail_slider.valueChanged.connect(state.set_tail_length)
        self._tail_label = QLabel(f"tail: {state.tail_length}")
        self._tail_slider.valueChanged.connect(lambda v: self._tail_label.setText(f"tail: {v}"))

        speeds = self.SPEED_VALUES
        self._speed_slider = QSlider(Qt.Orientation.Horizontal)
        self._speed_slider.setRange(0, len(speeds) - 1)
        try:
            initial_idx = speeds.index(state.compute_batch)
        except ValueError:
            initial_idx = 0
        self._speed_slider.setValue(initial_idx)
        self._speed_slider.valueChanged.connect(self._on_speed_changed)
        self._speed_label = QLabel(f"speed: {state.compute_batch}×")

        self._phase_combo = QComboBox()
        self._populate_phase_options()
        self._phase_combo.currentIndexChanged.connect(self._on_phase_selected)

        row = QHBoxLayout()
        row.addWidget(self._play_btn)
        row.addWidget(self._reset_btn)
        row.addSpacing(12)
        row.addWidget(self._tail_label)
        row.addWidget(self._tail_slider, 2)
        row.addSpacing(12)
        row.addWidget(self._speed_label)
        row.addWidget(self._speed_slider, 2)
        row.addSpacing(12)
        row.addWidget(QLabel("phase:"))
        row.addWidget(self._phase_combo, 1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addLayout(row)

    # ------------------------------------------------------------------
    # Populate / react
    # ------------------------------------------------------------------

    def _populate_phase_options(self) -> None:
        n = self.state.bridge.n
        dims = self.state.bridge.dims
        self._phase_entries: list[PhaseSelection] = []
        current = self.state.phase

        # Pair entries
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                sel = PhaseSelection(mode="pair", pair=(i, j),
                                     particle=current.particle, component=current.component)
                self._phase_entries.append(sel)
                self._phase_combo.addItem(sel.label())

        # Particle / component entries
        comp_labels = ["x", "y", "z"]
        for p in range(1, n + 1):
            for d in range(1, dims + 1):
                sel = PhaseSelection(
                    mode="particle", pair=current.pair, particle=p, component=d,
                )
                self._phase_entries.append(sel)
                self._phase_combo.addItem(f"Particle {p} ({comp_labels[d - 1]})")

        # Pick the entry matching current state
        for idx, sel in enumerate(self._phase_entries):
            if sel.mode == current.mode and (
                (sel.mode == "pair" and sel.pair == current.pair)
                or (sel.mode == "particle" and sel.particle == current.particle
                    and sel.component == current.component)
            ):
                self._phase_combo.setCurrentIndex(idx)
                break

    def _on_phase_selected(self, idx: int) -> None:
        if 0 <= idx < len(self._phase_entries):
            self.state.set_phase(self._phase_entries[idx])

    def _on_speed_changed(self, idx: int) -> None:
        speed = self.SPEED_VALUES[idx]
        self.state.set_compute_batch(speed)
        self._speed_label.setText(f"speed: {speed}×")

    def _on_play_state(self, playing: bool) -> None:
        self._play_btn.setText("⏸ Pause" if playing else "▶ Play")
