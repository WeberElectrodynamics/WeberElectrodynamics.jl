"""Streaming animation of a compact 2D three-body charge polygon."""

from __future__ import annotations

from weber_viewer.examples.scenarios import three_body_polygon_config


def main() -> None:
    from weber_viewer import animate_weber

    animate_weber(**three_body_polygon_config())


if __name__ == "__main__":
    main()
