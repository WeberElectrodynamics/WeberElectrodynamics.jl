# WeberElectrodynamics.jl

A Julia package for n-body simulation of charged particles interacting via
Weber's electrodynamics, with optional Zöllner electrogravitational extension.

## Features

- **Symplectic integrator** — Strang-splitting with symmetric projection
  (Jayawardana & Ohsawa 2021) for Weber's velocity-dependent Hamiltonian.
- **Regularization** — Levi-Civita/KS for close encounters in 2D; adaptive
  Cartesian substeps in 3D. Chain regularization for multi-particle encounters.
- **Zöllner extension** — unlike-sign pairs receive κ = 1 + a, producing an
  emergent gravitational correction.
- **Statistics** — energy, force, momentum, and trajectory analysis.
- **Visualization** — Plots.jl extension for static figures; Makie extension
  for interactive real-time and replay animation.

## Installation

```julia
using Pkg
Pkg.add("WeberElectrodynamics")
```

## Quick navigation

- [Quick Start](@ref) — minimal working example
- [API Reference](@ref "System") — complete type and function documentation
- [Theory](@ref) — links to derivations and design documents
