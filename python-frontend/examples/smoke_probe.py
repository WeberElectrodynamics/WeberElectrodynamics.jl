"""Short headless probe: bootstrap bridge, tick 20 frames, print diagnostics.

Intended for CI / sanity-checking the Python↔Julia bridge without opening a
window. Run with ``QT_QPA_PLATFORM=offscreen`` via the wrapper script.
"""

from __future__ import annotations

from PyQt6.QtCore import QCoreApplication

from weber_viewer import animate_weber


def main() -> None:
    window = animate_weber(
        n_particles=2,
        dims=3,
        masses=[1.0, 1.0],
        charges=[1.0, -1.0],
        q_initial=[1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
        p_initial=[0.0, 0.5, 0.0, 0.0, -0.5, 0.0],
        c=100.0,
        dt=0.01,
        tspan=(0.0, 10.0),
        compute_batch=1,
        autoplay=True,
        exec_app=False,
    )
    state = window.state
    for _ in range(20):
        state._tick()
        QCoreApplication.processEvents()

    print(f"OK: {state.total_steps} steps, buffer count={state.buffer.count}")
    print(f"    energy  : {state.energy_error_text()}")
    print(f"    time    : {state.time_text()}")
    window.close()


if __name__ == "__main__":
    main()
