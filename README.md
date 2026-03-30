# WeberElectrodynamics.jl

[![CI](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/actions/workflows/CI.yml)
[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl/stable)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl/dev)
[![Coverage](https://codecov.io/gh/WeberElectrodynamics/WeberElectrodynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WeberElectrodynamics/WeberElectrodynamics.jl)
[![Julia ≥1.9](https://img.shields.io/badge/julia-%E2%89%A51.9-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19239678.svg)](https://doi.org/10.5281/zenodo.19239678)

A Julia package for symplectic numerical integration of n-body Weber electrodynamics.

## Paper

**Computational Weber electrodynamics: symplectic n-body integration** — [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19337293.svg)](https://doi.org/10.5281/zenodo.19337293)

## Theory

- [Weber Electrodynamics](research/theory/WeberElectrodynamics.md) — potential, force, Lagrangian and Hamiltonian formulation, equations of motion
- [Semi-Explicit Symplectic Integrator](research/theory/SemiExplicitIntegrator.md) — extended phase space, Strang splitting, symmetric projection algorithm

## Extensions

- **WeberElectrodynamicsPlotsExt** — static plotting of trajectories, energy, forces, momentum and phase space
- **WeberElectrodynamicsMakieExt** — real-time animated dashboard with any Makie backend (GLMakie, CairoMakie, WGLMakie)

## Examples

See the [examples/](examples/) directory for Jupyter notebooks demonstrating two-body, three-body and four-body simulations with regularization, critical radius dynamics and Zollner electrogravity.

## Installation

```julia
using Pkg
Pkg.add("WeberElectrodynamics")
```

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
