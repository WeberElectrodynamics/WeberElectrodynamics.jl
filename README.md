# WeberElectrodynamics.jl

[![CI](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/actions/workflows/CI.yml)
[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl/stable)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://WeberElectrodynamics.github.io/WeberElectrodynamics.jl/dev)
[![Coverage](https://codecov.io/gh/WeberElectrodynamics/WeberElectrodynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WeberElectrodynamics/WeberElectrodynamics.jl)
[![Julia ≥1.9](https://img.shields.io/badge/julia-%E2%89%A51.9-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19239678.svg)](https://doi.org/10.5281/zenodo.19239678)

A Julia package for n-body simulation of charged particles interacting via Weber's electrodynamics. Implements a symplectic Strang-splitting integrator with symmetric projection for Weber's non-separable Hamiltonian, with Levi-Civita/KS regularization for close encounters. The integrator follows [Jayawardana & Ohsawa (2021)](https://arxiv.org/abs/2111.10915).

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

MIT
