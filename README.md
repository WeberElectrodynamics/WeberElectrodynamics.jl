# WeberElectrodynamics.jl

A Julia package for n-body simulation of charged particles interacting via [Weber's electrodynamics](https://arxiv.org/abs/2111.10915). Implements a symplectic Strang-splitting integrator with symmetric projection for Weber's non-separable Hamiltonian, with Levi-Civita/KS regularization for close encounters.

## Theory

- [Weber Electrodynamics](docs/theory/WeberElectrodynamics.md) — potential, force, Lagrangian and Hamiltonian formulation, equations of motion
- [Semi-Explicit Symplectic Integrator](docs/theory/SemiExplicitIntegrator.md) — extended phase space, Strang splitting, symmetric projection algorithm

## Extensions

- **WeberElectrodynamicsPlotsExt** — static plotting of trajectories, energy, forces, momentum and phase space (loaded automatically when [Plots.jl](https://github.com/JuliaPlots/Plots.jl) is imported)
- **WeberElectrodynamicsMakieExt** — real-time animated dashboard with any Makie backend (GLMakie, CairoMakie, WGLMakie) (loaded automatically when [Makie.jl](https://github.com/MakieOrg/Makie.jl) is imported)

## Examples

See the [examples/](examples/) directory for Jupyter notebooks demonstrating two-body, three-body and four-body simulations with regularization, critical radius dynamics and Zollner electrogravity.

## Installation

```julia
using Pkg
Pkg.add("WeberElectrodynamics")
```

## License

MIT
