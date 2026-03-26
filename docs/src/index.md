# WeberElectrodynamics.jl

A Julia package for n-body simulation of charged particles interacting via Weber's electrodynamics.

## Features

- **Symplectic integrator** — Extended phase space semi-explicit symplectic integrator for non-separable Hamiltonians with Strang-splitting and symmetric projection ([Jayawardana & Ohsawa 2021](https://arxiv.org/abs/2111.10915)).
- **Regularization** — Levi-Civita/KS for close encounters in 2D; adaptive
  Cartesian substeps in 3D. Chain regularization for multi-particle encounters.
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

