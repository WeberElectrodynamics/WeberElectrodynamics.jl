"""Headless smoke test: build MainWindow, tick, verify buffer growth.

Requires ``QT_QPA_PLATFORM=offscreen`` for CI. Requires juliacall.
"""

import os

import pytest

pytest.importorskip("juliacall")
pytest.importorskip("PyQt6")


@pytest.fixture(scope="module", autouse=True)
def _qt_offscreen():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_mainwindow_builds_and_ticks(qapp_args=None):
    # qapp_args preserved for pytest-qt compatibility; we use our own app here
    # because most envs don't ship pytest-qt.
    from PyQt6.QtCore import QCoreApplication
    from PyQt6.QtWidgets import QApplication

    from weber_viewer.bridge import JuliaBridge
    from weber_viewer.buffer import RollingBuffer
    from weber_viewer.state import AnimationState, PhaseSelection
    from weber_viewer.app import MainWindow

    app = QApplication.instance() or QApplication([])

    bridge = JuliaBridge.from_problem(
        n_particles=2,
        dims=3,
        q_initial=[1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
        p_initial=[0.0, 0.5, 0.0, 0.0, -0.5, 0.0],
        masses=[1.0, 1.0],
        charges=[1.0, -1.0],
        c=100.0,
        dt=0.01,
        tspan=(0.0, 1.0),
    )
    buffer = RollingBuffer(100, 2, 3)
    phase = PhaseSelection(mode="pair", pair=(1, 2), particle=1, component=1)
    state = AnimationState(
        bridge, buffer, tail_length=50, compute_batch=1, phase=phase,
    )
    win = MainWindow(state, figure_size=(800, 600))

    # Manually tick 20 frames synchronously — no event loop needed.
    state.set_playing(True)
    for _ in range(20):
        state._tick()
        QCoreApplication.processEvents()

    assert state.buffer.count >= 20
    assert state.total_steps == 20
    win.close()
