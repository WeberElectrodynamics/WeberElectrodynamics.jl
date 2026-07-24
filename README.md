# WeberElectrodynamics.jl

[![CI](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/actions/workflows/CI.yml)
[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl/stable)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl/dev)
[![Coverage](https://codecov.io/gh/WeberElectrodynamics/WeberElectrodynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WeberElectrodynamics/WeberElectrodynamics.jl)
[![Julia ≥1.9](https://img.shields.io/badge/julia-%E2%89%A51.9-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19239678-blue)](https://doi.org/10.5281/zenodo.19239678)

A Julia package for symplectic numerical integration of n-body Weber electrodynamics.

## Paper

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19337293-blue)](https://doi.org/10.5281/zenodo.19337293) [![Paper License: CC BY 4.0](https://img.shields.io/badge/Paper%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**[Computational Weber electrodynamics: symplectic n-body integration](https://doi.org/10.5281/zenodo.19337293)**

## Theory

- [Weber Electrodynamics](theory/WeberElectrodynamics.md) — potential, force, Lagrangian and Hamiltonian formulation, equations of motion
- [Semi-Explicit Symplectic Integrator](theory/SemiExplicitIntegrator.md) — extended phase space, Strang splitting, symmetric projection algorithm

## Extensions

- **WeberElectrodynamicsPlotsExt** — static plotting of trajectories, energy, forces, momentum and phase space
- **WeberElectrodynamicsMakieExt** — real-time animated dashboard with any Makie backend (GLMakie, CairoMakie, WGLMakie)
- **WeberElectrodynamicsJLD2Ext** — optional `HamiltonianSolution` archives via JLD2

## Regularization capability matrix

Regularization is opt-in through `RegularizedIntegrator`. The lifted backend
regularizes binary close encounters by dimension; chain encounters still use
adaptive Cartesian sub-stepping.

| Situation | Backend | Status |
| --- | --- | --- |
| 1D binary close encounter | `:lifted_pair` | Lifted square-root chart |
| 2D binary close encounter | `:lifted_pair` | Lifted Levi-Civita chart |
| 3D binary close encounter | `:lifted_pair` | Lifted KS chart with constraint projection diagnostics |
| Any dimension binary fallback | `:adaptive_cartesian` | Cartesian close-encounter substeps |
| Multi-particle close cluster | chain mode | Adaptive Cartesian substeps over the active component |

The regularization machinery handles the Coulomb/Kepler singular part. Weber's
velocity-dependent correction is evaluated through the existing equations of
motion and is not analytically regularized.

## Examples

The [examples/](examples/) directory contains:

- [`two_body_reference.ipynb`](examples/two_body_reference.ipynb) — the canonical annotated two-body tutorial with every diagnostic plot.
- [`api_showcase.ipynb`](examples/api_showcase.ipynb) — an API tour covering custom Hamiltonians, term introspection, accessors, regularization, callbacks, and plotting.

Additional **exploratory** studies live under [`_research/`](_research/). The
folder is a research sandbox — its contents are not definitive and should not
be cited as package behaviour. See [`_research/README.md`](_research/README.md).

## Installation

```julia
using Pkg
Pkg.add("WeberElectrodynamics")
```

## Quick start

The canonical two-body reference problem — symmetric whole-number initial
conditions producing a precessing ellipse ($e = 3/4$) over five Kepler
periods.

```julia
using WeberElectrodynamics, Plots

system = HamiltonianSystem(2, 2)
prob = HamiltonianProblem(system, (0.0, 27.14),
    [-1.0, 0.0,  1.0, 0.0],        # q(0): particles at (±1, 0)
    [ 0.0, -0.25, 0.0, 0.25];      # p(0): |p| = 1/4 tangential
    masses = [1.0, 1.0], charges = [1.0, -1.0], c = 4.0, dt = 0.001,
)
sol = solve(prob)

plot_trajectories(compute_trajectory_data(sol, 2, 2; stride = 10))
energy = compute_energy_timeseries(sol; stride = 10)
plot_energy(energy)
plot_energy_errors(energy)
```

![Two-body trajectories](examples/figures/two_body_trajectories.png)
![Energy conservation](examples/figures/two_body_energy.png)
![Energy error diagnostics](examples/figures/two_body_energy_errors.png)

See [`examples/two_body_reference.ipynb`](examples/two_body_reference.ipynb)
for the full annotated tutorial with every diagnostic plot.

## License

[MIT](LICENSE)
