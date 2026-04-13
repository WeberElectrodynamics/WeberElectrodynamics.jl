# Internals

Reference for contributors, developers, and AI agents working with WeberElectrodynamics.jl internals.

## Architecture

**Pipeline**: Symbolic Hamiltonian → compiled equations of motion → symplectic integration → statistics/plotting/animation

### Source files (`src/`)

- `WeberElectrodynamics.jl` — Module definition, exports, extension stubs (`plot_*`, `animate_weber`)
- `weber_system.jl` — `WeberSystem`: uses Symbolics.jl to build the Weber Hamiltonian symbolically, then compiles `dq_dt`, `dp_dt`, and `hamiltonian` functions via `build_function`
- `types.jl` — All core structs: `WeberProblem`, `WeberSolution`, `WeberIntegrator`, `SymmetricProjectionIntegrator`, `RegularizationOptions`, `ZollnerOptions`, buffer/diagnostics types
- `regularization.jl` — Internal helpers: pair distance detection, adjacency graph (BFS), Levi-Civita 2D projection, KS quaternion helpers
- `solve.jl` — Main integrator: Strang splitting flow, symmetric projection via fixed-point iteration on Lagrange multipliers, regularization dispatch, collision bounce, CommonSolve interface (`init`/`step!`/`solve!`/`solve`)
- `statistics/` — `energy.jl`, `forces.jl`, `momentum.jl`, `trajectories.jl` — post-solution analysis producing typed data structs

### Extensions (`ext/`)

- `WeberElectrodynamicsPlotsExt.jl` — Plots.jl weak dependency; provides `plot_trajectories`, `plot_energy`, `plot_pair_energy`, `plot_energy_errors`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum_errors`, and Zöllner-specific plot functions.
- `WeberElectrodynamicsMakieExt.jl` — Makie weak dependency (any backend: GLMakie, CairoMakie, WGLMakie); provides `animate_weber` for real-time streaming or solution replay with rolling trajectory/energy/momentum/phase-space dashboard.

### Tests (`test/`)

- `test_utils.jl` — Problem builders (`make_weber_problem()`, `make_coulomb_like_problem()`) and reference energy functions; **must be included before other test files**
- `runtests.jl` — Entry point, includes all test files in order
- Test files: `test_types.jl`, `test_weber_system.jl`, `test_solve.jl`, `test_statistics.jl`, `test_integration.jl`, `test_physics.jl`, `test_regularization.jl`, `test_zollner.jl`

### Examples (`examples/`)

Jupyter notebooks run via IJulia from the default Julia environment (where this package is `dev`'d). They use `Plots` for static figures and optionally `GLMakie` for animation.

### Docs (`docs/`)

Documenter.jl scaffold: `make.jl`, `Project.toml`, `src/` (page sources), `build/` (generated output).

### Research (`research/`)

- `theory/` — Weber electrodynamics, semi-explicit integrator, regularization, critical radius, initial conditions, Zöllner theory
- `exploratory/` — Collision bounce lessons learned, three-body bound states
- `sub_critical_weber_research/` — Sub-critical Weber exploration and literature searches

---

## Internal Conventions

### Params vector layout

```
params = [m₁, ..., mₙ, q₁, ..., qₙ, c, κ₁₂, κ₁₃, ..., κ_{N-1,N}]
length = 2N + 1 + N*(N-1)/2
```

Any code calling `sys.dq_dt_compiled(out, q, p, params)` or `sys.dp_dt_compiled(out, q, p, params)` directly **must** include the κ (kappa) entries. When Zöllner is disabled, all κ values are 1.0.

Pair index: `_pair_index(i, j, n) = (i-1)*(2n-i)÷2 + (j-i)` (1-based, i < j)

### Regularization backends

Only two valid values for `RegularizationOptions.backend`:
- `:adaptive_cartesian` — KS-style, works for 2D and 3D
- `:lifted_pair` — Levi-Civita, **2D only** (auto-falls back to `:adaptive_cartesian` for 3D)

Neither backend regularizes Weber's velocity-dependent force — only the Coulomb/Kepler singularity.

### Collision bounce

- Enabled via `RegularizationOptions(collision_bounce_radius = <r>)` passed to `WeberProblem` (default 0.0 = off)
- Only valid for ℓ=0 (head-on) collisions
- Works best with the **unregularized** integrator (symplectic error stays bounded)

### Zöllner extension

- `ZollnerOptions(enabled, a)` — mismatch parameter `a`
- κ_ij = 1+a for unlike-sign charge pairs, 1.0 for like-sign
- Stored in `WeberProblem.kappas` and appended to the params vector automatically

### Makie animation extension

- Weak dependency is `Makie` (not `GLMakie`) — any backend triggers the extension
- `animate_weber(prob)` for live streaming, `animate_weber(sol)` for replay
- Compat: `Makie = "0.21, 0.22, 0.23, 0.24"`

### Immutable options pattern

`RegularizationOptions`, `ZollnerOptions` are immutable structs created once per problem. Pass configuration through `WeberProblem` keyword arguments rather than mutating options.

---

## EnergyStatistics fields

`en.statistics` has: `local_error_max`, `local_error_min`, `local_error_avg`, `global_error_ratio_max/min/avg`, `global_error_percent_max/min/avg`. There is **no** `local_error_percent_max` — use `local_error_max` for local error magnitude.
