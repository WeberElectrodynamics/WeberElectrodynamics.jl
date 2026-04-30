# weber-viewer

Python animation frontend for the [WeberElectrodynamics.jl](../) integrator.

Mirrors the Makie extension's `animate_weber(prob)` streaming mode as a PyQt6
desktop application. The Julia package runs the integrator in-process via
`juliacall`; Python owns the UI.

## Install

```bash
# 1. Python package (editable)
pip install -e python-frontend/

# 2. Point juliacall at this repo's Julia project so the local WeberElectrodynamics
#    package is used instead of a registered version
export PYTHON_JULIAPKG_PROJECT="$(pwd)"
export PYTHON_JULIAPKG_OFFLINE=yes    # optional: skip package resolution
```

## Run

One-shot wrapper (bootstraps a venv on first run, then just runs):

```bash
./python-frontend/run-two-body.sh                 # 3D two-body viewer
./python-frontend/run-three-body.sh               # 2D three-body polygon viewer
./python-frontend/run-smoke.sh --steps 50         # headless Julia bridge probe

./python-frontend/run-example.sh                 # two-body viewer (default)
./python-frontend/run-example.sh --list          # list bundled examples
./python-frontend/run-example.sh two-body        # 3D two-body viewer
./python-frontend/run-example.sh three-body      # 2D three-body polygon viewer
./python-frontend/run-example.sh smoke           # headless Julia bridge probe
./python-frontend/run-example.sh smoke --steps 50
./python-frontend/run-example.sh path/to/other.py
```

The wrapper keeps JuliaCall's scratch dependencies in `.venv-viewer/julia-env`
and `dev`s the local package there, so the repository `Project.toml` is not
modified by example runs.

Or manually (after `pip install -e python-frontend/`):

```bash
export PYTHON_JULIAPKG_PROJECT="$(pwd)"
export PYTHON_JULIAPKG_EXE="$(which julia)"
python python-frontend/examples/two_body_streaming.py
python python-frontend/examples/three_body_polygon.py
python python-frontend/examples/smoke_probe.py --steps 25
```

Installed console scripts are also available:

```bash
weber-viewer-two-body
weber-viewer-three-body
weber-viewer-smoke --steps 25
```

Or from Python:

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
