"""Streaming animation of a circular 3D two-body orbit.

Mirrors ``examples/two_body_reference.ipynb`` from the Julia side. Run:

    PYTHON_JULIAPKG_PROJECT=$(pwd) python python-frontend/examples/two_body_streaming.py
"""

from __future__ import annotations

from weber_viewer import animate_weber


def main() -> None:
    # Two equal-mass, opposite-charge particles on a circular orbit in the x-y plane,
    # offset slightly out of plane so the 3D view is visually meaningful.
    animate_weber(
        n_particles=2,
        dims=3,
        masses=[1.0, 1.0],
        charges=[1.0, -1.0],
        q_initial=[1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
        p_initial=[0.0, 0.5, 0.0, 0.0, -0.5, 0.0],
        c=100.0,
        dt=0.01,
        tspan=(0.0, 200.0),
        buffer_size=2000,
        tail_length=400,
        compute_batch=2,
    )


if __name__ == "__main__":
    main()
