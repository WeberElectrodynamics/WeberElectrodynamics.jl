# WeberElectrodynamics.jl

A Julia package for n-body simulation of charged particles interacting via Weber's electrodynamics. Implements a symplectic Strang-splitting integrator with symmetric projection for Weber's non-separable Hamiltonian, with Levi-Civita/KS regularization for close encounters. The integrator follows [Jayawardana & Ohsawa (2021)](https://arxiv.org/abs/2111.10915).

## Theory

- [Weber Electrodynamics](docs/theory/WeberElectrodynamics.md) — potential, force, Lagrangian and Hamiltonian formulation, equations of motion
- [Semi-Explicit Symplectic Integrator](docs/theory/SemiExplicitIntegrator.md) — extended phase space, Strang splitting, symmetric projection algorithm

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
