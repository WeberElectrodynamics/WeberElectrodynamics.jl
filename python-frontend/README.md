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
./python-frontend/run-example.sh                             # two_body_streaming.py
./python-frontend/run-example.sh examples/other_script.py    # any other example
```

Or manually (after `pip install -e python-frontend/`):

```bash
export PYTHON_JULIAPKG_PROJECT="$(pwd)"
export PYTHON_JULIAPKG_EXE="$(which julia)"
python python-frontend/examples/two_body_streaming.py
```

Or from Python:

```python
from weber_viewer import animate_weber

animate_weber(
    n_particles=2,
    dims=3,
    masses=[1.0, 1.0],
    charges=[1.0, -1.0],
    q_initial=[1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
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
                  ├── JuliaBridge.step()      ── juliacall ──► WeberIntegrator
                  └── RollingBuffer.push_step()
                       └── emits buffer_updated → widgets refresh
```
