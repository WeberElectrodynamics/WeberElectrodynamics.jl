# Internals

Reference for contributors, developers, and AI agents working with WeberElectrodynamics.jl internals.

## Architecture

**Pipeline**: Symbolic Hamiltonian → compiled equations of motion → symplectic integration → statistics/plotting/animation

### Source files (`src/`)

- `WeberElectrodynamics.jl` — Module definition, exports, extension stubs (`plot_*`, `animate_weber`, archive helpers)
- `initial_conditions.jl` — COM-frame two-body, polygon, and rigid-rotation initial-condition helpers
- `hamiltonian_system.jl` — `HamiltonianSystem`: uses Symbolics.jl to build the Weber Hamiltonian symbolically, then compiles `dq_dt`, `dp_dt`, and `hamiltonian` functions via `build_function`
- `types.jl` — All core structs: `HamiltonianProblem`, `HamiltonianSolution`, `HamiltonianIntegrator`, `SymmetricProjectionIntegrator`, `RegularizationOptions`, `ZollnerOptions`, buffer/diagnostics types
- `regularization.jl` — Internal helpers: pair distance detection, adjacency graph (BFS), Levi-Civita 2D projection, KS quaternion helpers
- `solve.jl` — Main integrator: Strang splitting flow, symmetric projection via fixed-point iteration on Lagrange multipliers, regularization dispatch, collision bounce, CommonSolve interface (`init`/`step!`/`solve!`/`solve`)
- `statistics/` — `energy.jl`, `forces.jl`, `momentum.jl`, `trajectories.jl` — post-solution analysis producing typed data structs

### Extensions (`ext/`)

- `WeberElectrodynamicsPlotsExt.jl` — Plots.jl weak dependency; provides `plot_trajectories`, `plot_energy`, `plot_pair_energy`, `plot_energy_errors`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum_errors`, and Zöllner-specific plot functions.
- `WeberElectrodynamicsMakieExt.jl` — Makie weak dependency (any backend: GLMakie, CairoMakie, WGLMakie); provides `animate_weber` for real-time streaming or solution replay with rolling trajectory/energy/momentum/phase-space dashboard.
- `WeberElectrodynamicsJLD2Ext.jl` — JLD2 weak dependency; provides `save_solution` and `load_solution`.

### Tests (`test/`)

- `test_utils.jl` — Problem builders (`make_weber_problem()`, `make_coulomb_like_problem()`) and reference energy functions; **must be included before other test files**
- `runtests.jl` — Entry point, includes all test files in order
- Test files: `test_types.jl`, `test_hamiltonian_system.jl`, `test_initial_conditions.jl`, `test_solve.jl`, `test_statistics.jl`, `test_integration.jl`, `test_physics.jl`, `test_regularization.jl`, `test_zollner.jl`

### Examples (`examples/`)

Jupyter notebooks run via IJulia from the default Julia environment (where this package is `dev`'d). They use `Plots` for static figures and optionally `GLMakie` for animation.

### Docs (`docs/`)

Documenter.jl scaffold: `make.jl`, `Project.toml`, `src/` (page sources), `build/` (generated output).

### Research sandbox (`_research/`)

> **Research sandbox — contents are exploratory, not definitive.**
> Do not cite material in this folder as package behaviour or authoritative theory.
> Reports, notes, and scripts here may be incomplete, superseded, or contradicted
> by code and theory elsewhere in the repo. See [`_research/README.md`](https://github.com/WeberElectrodynamics/WeberElectrodynamics/blob/main/_research/README.md).

- `Topology/FourBodyTwoPlusTwoMinus/` — 4-body (2+/2−) multi-agent investigation (reports only)
- `Topology/Homology/` — Rabinowitz-Floer homology study (reports only)
- `Investigations/` — geometric/topological questions, collision bounce, three-body, hypergeometric structure, sub-critical exploration
- `Notebooks/` — Jupyter reference and scratch notebooks
- `LiteratureSearches/` — bibliographic TOML snapshots
- `NextSteps.md` — forward-looking notes on promotion candidates and theory/code gaps

---

## Internal Conventions

### Params and κ vector layout

```
params = [m₁, ..., mₙ, q₁, ..., qₙ, c]                     # length 2N + 1
kappas = [κ₁₂, κ₁₃, ..., κ_{N-1,N}]                         # length N*(N-1)/2
```

Compiled EOM signature (all three accept the same 6 positional args):

```
sys.dq_dt_compiled(out, q, p, t, params, kappas)
sys.dp_dt_compiled(out, q, p, t, params, kappas)
sys.hamiltonian_compiled(q, p, t, params, kappas)
```

Any code calling these directly **must** pass both `params` and `kappas`. When Zöllner is disabled, `kappas = ones(N*(N-1)÷2)`. For single-particle systems (no pairs), `kappas = Float64[]`.

On `HamiltonianProblem`, use `kappas(prob)` to obtain the whole vector, or `kappa(prob, i, j)` for a single pair. `kappa` canonicalizes pair order and is preferred over manual `_pair_index` arithmetic in downstream code.

Pair index (internal): `_pair_index(i, j, n) = (i-1)*(2n-i)÷2 + (j-i)` (1-based, i < j)

### Regularization backends

Only two valid values for `RegularizationOptions.backend`:
- `:adaptive_cartesian` — Cartesian close-encounter substeps, works for all dimensions
- `:lifted_pair` — lifted square-root (1D), Levi-Civita (2D), or KS (3D) binary pair stepping

Multi-particle close clusters use chain mode, which is adaptive Cartesian over
the active component rather than analytic chain-coordinate regularization.
Neither backend analytically regularizes Weber's velocity-dependent force —
only the Coulomb/Kepler singularity.

### Collision bounce

- Implemented as the `CollisionBounce(radius)` callback; pass it to `solve` (or `init`) via the `callbacks` kwarg
- `RegularizedIntegrator` also accepts a `collision_bounce_radius` kwarg and synthesises a matching callback automatically when no `CollisionBounce` is supplied
- Under `RegularizedIntegrator`, the bounce radius is checked after regularized substeps as well as at macro-step boundaries
- Geometric and charge-sign agnostic: it reflects any pair inside the radius, including two-body unlike-charge pass-through when explicitly enabled
- Only valid for `L = 0` head-on, C0-continuable collisions/pass-through events; it is not generally energy-preserving for `N > 2`

### Zöllner extension

- `ZollnerOptions(enabled, a)` — mismatch parameter `a`
- κ_ij = 1+a for opposite nonzero charge signs, 1.0 for like-sign or neutral pairs
- Computed at `HamiltonianProblem` construction and stored separately from `params`; read via the `kappas(prob)` accessor

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
