"""Streaming animation of a circular 3D two-body orbit."""

from __future__ import annotations

from weber_viewer.examples.scenarios import two_body_config


def main() -> None:
    from weber_viewer import animate_weber

    animate_weber(**two_body_config())


if __name__ == "__main__":
    main()
