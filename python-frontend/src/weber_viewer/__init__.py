"""Python animation frontend for the WeberElectrodynamics.jl integrator."""

from __future__ import annotations

__all__ = ["animate_weber"]
__version__ = "0.1.0"


def __getattr__(name):
    # Lazy-import so importing ``weber_viewer`` (e.g. for ``RollingBuffer``) does
    # not pull in PyQt6 / juliacall unless the top-level entry point is used.
    if name == "animate_weber":
        from weber_viewer.app import animate_weber as _fn
        return _fn
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
