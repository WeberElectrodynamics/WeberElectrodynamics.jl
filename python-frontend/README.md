# weber-viewer

Python animation frontend for the [WeberElectrodynamics.jl](../) integrator.

Mirrors the Makie extension's `animate_weber(prob)` streaming mode as a PyQt6
desktop application. The Julia package runs the integrator in-process via
`juliacall`; Python owns the UI.

## Requirements

[uv](https://docs.astral.sh/uv/) and `julia` on `PATH`. uv is the only supported
Python workflow here — it manages the interpreter, the virtualenv and the locked
dependency set. There is nothing to `pip install`.

## Run

All commands run from `python-frontend/`:

```bash
uv run weber-viewer-two-body           # 3D two-body streaming viewer
uv run weber-viewer-three-body         # 2D three-body polygon viewer
uv run weber-viewer-smoke --steps 50   # headless Julia bridge probe
```

The first invocation creates `.venv/` from `uv.lock` and bootstraps a scratch
Julia project at `.julia-env/` that `dev`s the local package, so juliacall can
add PythonCall/OpenSSL_jll without touching the repository's own `Project.toml`.
Both directories are gitignored and regenerate on demand.

Tests:

```bash
uv run pytest
```

### Interpreter pin

`.python-version` pins CPython 3.14. juliacall constrains `OpenSSL_jll` to
whatever OpenSSL the interpreter links against, and uv's CPython 3.13 build
links OpenSSL 3.0, for which no installable `OpenSSL_jll` remains — resolution
fails. The 3.14 build links OpenSSL 3.6 and resolves cleanly.

Override the scratch Julia project location with `WEBER_VIEWER_JULIA_ENV`, or
bypass detection entirely with the standard `PYTHON_JULIAPKG_PROJECT` /
`PYTHON_JULIAPKG_EXE` variables.

## Use as a library

```bash
uv run python
```



```python
from weber_viewer import animate_weber

animate_weber(
    n_particles=2,
    dims=3,
    masses=[1.0, 1.0],
    charges=[1.0, -1.0],
    q_initial=[1.0, 0.0, 0.15, -1.0, 0.0, -0.15],
    p_initial=[0.0, 0.5, 0.0, 0.0, -0.5, 0.0],
    c=100.0,
    dt=0.01,
    tspan=(0.0, 100.0),
)
```

## Architecture

```
PyQt6 MainWindow
  ├── TrajectoryWidget  (3D GLViewWidget or 2D PlotWidget)
  ├── PhaseSpaceWidget  (2D PlotWidget, pair or particle mode)
  ├── HUDWidget         (energy / time / step labels)
  └── ControlsWidget    (play/pause, reset, sliders, dropdown)

QTimer @ 60 Hz → AnimationState._tick()
                  ├── JuliaBridge.step()      ── juliacall ──► HamiltonianIntegrator
                  └── RollingBuffer.push_step()
                       └── emits buffer_updated → widgets refresh
```
