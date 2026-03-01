# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Julia package for n-body Weber electrodynamics simulation with Zöllner electrogravitational extension. Implements a symplectic Strang-splitting symmetric-projection integrator with Levi-Civita/KS regularization for close encounters.

## Commands

```bash
# Run full test suite
julia -e 'using Pkg; Pkg.test()'

# Run a single test file (e.g., test_physics.jl)
julia -e 'using Test; using WeberElectrodynamics; using WeberElectrodynamics: SymmetricProjectionIntegrator; using LinearAlgebra; using Symbolics; @testset "single" begin include("test/test_utils.jl"); include("test/test_physics.jl") end'

# Format all Julia files
make format   # requires JuliaFormatter
```

## Architecture

**Pipeline**: Symbolic Hamiltonian → compiled equations of motion → symplectic integration → statistics/plotting

### Source files (`src/`)

- `WeberElectrodynamics.jl` — Module definition, exports, plot extension stubs
- `weber_system.jl` — `WeberSystem`: uses Symbolics.jl to build the Weber Hamiltonian symbolically, then compiles `dq_dt`, `dp_dt`, and `hamiltonian` functions via `build_function`
- `types.jl` — All core structs: `WeberProblem`, `WeberSolution`, `WeberIntegrator`, `SymmetricProjectionIntegrator`, `RegularizationOptions`, `ZollnerOptions`, buffer/diagnostics types
- `regularization.jl` — Internal helpers: pair distance detection, adjacency graph (BFS), Levi-Civita 2D projection, KS quaternion helpers
- `solve.jl` — Main integrator: Strang splitting flow, symmetric projection via fixed-point iteration on Lagrange multipliers, regularization dispatch, collision bounce, CommonSolve interface (`init`/`step!`/`solve!`/`solve`)
- `statistics/` — `energy.jl`, `forces.jl`, `momentum.jl`, `trajectories.jl` — post-solution analysis producing typed data structs

### Extension (`ext/`)

- `WeberElectrodynamicsPlotsExt.jl` — Plots.jl weak dependency; provides `plot_trajectories`, `plot_energy`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum`, and Zöllner-specific plot functions. Stubs declared in main module.

### Tests (`test/`)

- `test_utils.jl` — Problem builders (`make_weber_problem()`, `make_coulomb_like_problem()`) and reference energy functions; **must be included before other test files**
- Test files cover: types, system generation, solve interface, statistics, integration, physics validation, regularization, Zöllner

## Critical Conventions

### Params vector layout

```
params = [m₁, ..., mₙ, q₁, ..., qₙ, c, κ₁₂, κ₁₃, ..., κ_{N-1,N}]
length = 2N + 1 + N*(N-1)/2
```

Any code calling `sys.dq_dt_compiled(out, q, p, params)` or `sys.dp_dt_compiled(out, q, p, params)` directly **must** include the κ (kappa) entries. When Zöllner is disabled, all κ values are 1.0.

Pair index: `_pair_index(i, j, n) = (i-1)*(2n-i)÷2 + (j-i)` (1-based, i < j)

### Regularization backends

Only two valid values for `regularization_backend`:
- `:adaptive_cartesian` — KS-style, works for 2D and 3D
- `:lifted_pair` — Levi-Civita, **2D only** (auto-falls back to `:adaptive_cartesian` for 3D)

Neither backend regularizes Weber's velocity-dependent force — only the Coulomb/Kepler singularity.

### Collision bounce

- Enabled via `regularization_collision_bounce_radius` kwarg on `WeberProblem` (default 0.0 = off)
- Only valid for ℓ=0 (head-on) collisions
- Works best with the **unregularized** integrator (symplectic error stays bounded)

### Zöllner extension

- `ZollnerOptions(enabled, a)` — mismatch parameter `a`
- κ_ij = 1+a for unlike-sign charge pairs, 1.0 for like-sign
- Stored in `WeberProblem.kappas` and appended to the params vector automatically

### Immutable options pattern

`RegularizationOptions`, `ZollnerOptions` are immutable structs created once per problem. Pass configuration through `WeberProblem` keyword arguments rather than mutating options.

## EnergyStatistics fields

`en.statistics` has: `local_error_max`, `local_error_min`, `local_error_avg`, `global_error_ratio_max/min/avg`, `global_error_percent_max/min/avg`. There is **no** `local_error_percent_max`.
