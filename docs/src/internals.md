# Internals

Reference for contributors, developers, and AI agents working with WeberElectrodynamics.jl internals.

## Architecture

**Pipeline**: Symbolic Hamiltonian → compiled equations of motion → symplectic integration → statistics/plotting/animation

### Source files (`src/`)

- `WeberElectrodynamics.jl` — Module definition, exports, extension stubs (`plot_*`, `animate_weber`)
- `hamiltonian_system.jl` — `HamiltonianSystem`: uses Symbolics.jl to build the Weber Hamiltonian symbolically, then compiles `dq_dt`, `dp_dt`, and `hamiltonian` functions via `build_function`
- `types.jl` — All core structs: `HamiltonianProblem`, `HamiltonianSolution`, `HamiltonianIntegrator`, `SymmetricProjectionIntegrator`, `RegularizationOptions`, `ZollnerOptions`, buffer/diagnostics types
- `regularization.jl` — Internal helpers: pair distance detection, adjacency graph (BFS), Levi-Civita 2D projection, KS quaternion helpers
- `solve.jl` — Main integrator: Strang splitting flow, symmetric projection via fixed-point iteration on Lagrange multipliers, regularization dispatch, collision bounce, CommonSolve interface (`init`/`step!`/`solve!`/`solve`)
- `statistics/` — `energy.jl`, `forces.jl`, `momentum.jl`, `trajectories.jl` — post-solution analysis producing typed data structs

### Extensions (`ext/`)

- `WeberElectrodynamicsPlotsExt.jl` — Plots.jl weak dependency; provides `plot_trajectories`, `plot_energy`, `plot_pair_energy`, `plot_energy_errors`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum_errors`, and Zöllner-specific plot functions.
- `WeberElectrodynamicsMakieExt.jl` — Makie weak dependency (any backend: GLMakie, CairoMakie, WGLMakie); provides `animate_weber` for real-time streaming or solution replay with rolling trajectory/energy/momentum/phase-space dashboard.

### Tests (`test/`)

- `test_utils.jl` — Problem builders (`make_weber_problem()`, `make_coulomb_like_problem()`) and reference energy functions; **must be included before other test files**
- `runtests.jl` — Entry point, includes all test files in order
- Test files: `test_types.jl`, `test_hamiltonian_system.jl`, `test_solve.jl`, `test_statistics.jl`, `test_integration.jl`, `test_physics.jl`, `test_regularization.jl`, `test_zollner.jl`

### Examples (`examples/`)

Jupyter notebooks run via IJulia from the default Julia environment (where this package is `dev`'d). They use `Plots` for static figures and optionally `GLMakie` for animation.

### Docs (`docs/`)

Documenter.jl scaffold: `make.jl`, `Project.toml`, `src/` (page sources), `build/` (generated output).

### Research sandbox (`_research/`)

> **Research sandbox — contents are exploratory, not definitive.**
> Do not cite material in this folder as package behaviour or authoritative theory.
> Reports, notes, and scripts here may be incomplete, superseded, or contradicted
> by code and theory elsewhere in the repo. See [`_research/README.md`](https://github.com/WeberElectrodynamics/WeberElectrodynamics/blob/main/_research/README.md).

- `FourBodyTwoPlusTwoMinus/` — 4-body (2+/2−) multi-agent investigation
- `homology/` — Rabinowitz-Floer homology study
- `exploratory/` — foundational lessons learned (collision bounce, three-body, …)
- `hypergeometric_analysis/` — special-function analyses
- `investigations/` — geometric/topological questions
- `notebooks/` — Jupyter reference and scratch notebooks
- `sub_critical_weber_research/` — sub-critical regime exploration + literature taxonomy
- `topology/` — phase-space topology notebooks

---

## Internal Conventions

### Params vector layout

```
params = [m₁, ..., mₙ, q₁, ..., qₙ, c, κ₁₂, κ₁₃, ..., κ_{N-1,N}]
length = 2N + 1 + N*(N-1)/2
```

Any code calling `sys.dq_dt_compiled(out, q, p, t, params)` or `sys.dp_dt_compiled(out, q, p, t, params)` directly **must** include the κ (kappa) entries. When Zöllner is disabled, all κ values are 1.0.

Pair index: `_pair_index(i, j, n) = (i-1)*(2n-i)÷2 + (j-i)` (1-based, i < j)

### Regularization backends

Only two valid values for `RegularizationOptions.backend`:
- `:adaptive_cartesian` — KS-style, works for 2D and 3D
- `:lifted_pair` — Levi-Civita, **2D only** (auto-falls back to `:adaptive_cartesian` for 3D)

Neither backend regularizes Weber's velocity-dependent force — only the Coulomb/Kepler singularity.

### Collision bounce

- Implemented as the `CollisionBounce(radius)` callback; pass it to `solve` (or `init`) via the `callbacks` kwarg
- `RegularizedIntegrator` also accepts a `collision_bounce_radius` kwarg and synthesises a matching callback automatically when no `CollisionBounce` is supplied
- Only valid for ℓ=0 (head-on) collisions
- Works best with the **unregularized** integrator (symplectic error stays bounded)

### Zöllner extension

- `ZollnerOptions(enabled, a)` — mismatch parameter `a`
- κ_ij = 1+a for unlike-sign charge pairs, 1.0 for like-sign
- Computed at `HamiltonianProblem` construction and packed into `params`; read via the `kappas(prob)` accessor

### Makie animation extension

- Weak dependency is `Makie` (not `GLMakie`) — any backend triggers the extension
- `animate_weber(prob)` for live streaming, `animate_weber(sol)` for replay
- Compat: `Makie = "0.21, 0.22, 0.23, 0.24"`

### Immutable options pattern

`RegularizationOptions` and `ZollnerOptions` are immutable structs — create them once, never mutate. They sit at different layers:

- `ZollnerOptions` is problem-level: pass it via the `zollner = …` kwarg of `HamiltonianProblem`. It is used at construction to compute per-pair κ values and is not retained on the problem.
- `RegularizationOptions` is algorithm-level: it lives inside `RegularizedIntegrator`. Use the `RegularizedIntegrator(base_alg; kwargs...)` constructor (its kwargs mirror `RegularizationOptions`) rather than building the options struct by hand.

---

## EnergyStatistics fields

`en.statistics` has: `local_error_max`, `local_error_min`, `local_error_avg`, `global_error_ratio_max/min/avg`, `global_error_percent_max/min/avg`. There is **no** `local_error_percent_max` — use `local_error_max` for local error magnitude.
