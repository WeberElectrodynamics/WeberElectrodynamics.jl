"""Headless bridge probe for CI and local sanity checks."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from weber_viewer.examples.scenarios import problem_kwargs, two_body_config


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run a headless Weber bridge smoke probe.")
    parser.add_argument("--steps", type=int, default=25, help="number of macro-steps to run")
    args = parser.parse_args(argv)
    if args.steps < 1:
        raise SystemExit("--steps must be positive")

    from weber_viewer.bridge import JuliaBridge

    bridge = JuliaBridge.from_problem(**problem_kwargs(two_body_config()))
    _, q0, p0 = bridge.snapshot()
    e0 = bridge.energy(q0, p0)

    steps_run = 0
    for _ in range(args.steps):
        if not bridge.step():
            break
        steps_run += 1

    t, q, p = bridge.snapshot()
    e = bridge.energy(q, p)
    rel = abs((e - e0) / e0) if abs(e0) > 1e-30 else abs(e - e0)
    print(f"OK: {steps_run} steps")
    print(f"    time     : {t:.6f}")
    print(f"    energy   : {e:.12e}")
    print(f"    rel drift: {rel:.3e}")


if __name__ == "__main__":
    main()
