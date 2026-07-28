# Internals

Reference for contributors, developers, and AI agents working with WeberElectrodynamics.jl internals.

## Architecture

**Pipeline**: Hamiltonian (analytic or symbolic) → compiled equations of motion → symplectic integration → statistics/plotting/animation

The default Weber system takes the **analytic** path: its exact canonical
Hamiltonian requires solving a coupled linear system over particle pairs at
every evaluation to recover physical velocities from canonical momenta, which
has no practical closed symbolic form for general `n`. Custom Hamiltonians take
the **symbolic** path through Symbolics.jl. `has_symbolic_hamiltonian(sys)`
distinguishes them; the three symbolic fields are `nothing` on the analytic
path. Both fill identical compiled slots, so nothing downstream branches.

### Source files (`src/`)

- `WeberElectrodynamics.jl` — Module definition, exports, extension stubs (`plot_*`, `animate_weber`, archive helpers)
- `initial_conditions.jl` — COM-frame two-body, polygon, and rigid-rotation initial-condition helpers
- `hamiltonian_system.jl` — `HamiltonianSystem`: routes the default Weber system to the analytic evaluators in `hamiltonian/builders/weber.jl`, and any custom symbolic `H` through Symbolics.jl `build_function`
- `hamiltonian/builders/weber.jl` — the exact canonical Weber core: `WeberWorkspace`, the pair radial-velocity solve (`_weber_state!`), `weber_dq_dt!`, `weber_dp_dt!`, `weber_hamiltonian`, and the public `physical_velocities`. **Every** Weber path obtains velocities here; nothing in the package may reconstruct a velocity as `p/m`
- `types.jl` — All core structs: `HamiltonianProblem`, `HamiltonianSolution`, `HamiltonianIntegrator`, `SymmetricProjectionIntegrator`, `RegularizationOptions`, buffer/diagnostics types
- `regularization.jl` — Internal helpers: pair distance detection, adjacency graph (BFS), Levi-Civita 2D projection, KS quaternion helpers
- `solve.jl` — Main integrator: Strang splitting flow, symmetric projection via fixed-point iteration on Lagrange multipliers, regularization dispatch, collision bounce, CommonSolve interface (`init`/`step!`/`solve!`/`solve`)
- `statistics/` — `energy.jl`, `forces.jl`, `momentum.jl`, `trajectories.jl` — post-solution analysis producing typed data structs

### Extensions (`ext/`)

- `WeberElectrodynamicsPlotsExt.jl` — Plots.jl weak dependency; provides `plot_trajectories`, `plot_energy`, `plot_pair_energy`, `plot_energy_errors`, `plot_pair_forces`, `plot_phase_space`, and `plot_momentum_errors`.
- `WeberElectrodynamicsMakieExt.jl` — Makie weak dependency (any backend: GLMakie, CairoMakie, WGLMakie); provides `animate_weber` for real-time streaming or solution replay with rolling trajectory/energy/momentum/phase-space dashboard.
- `WeberElectrodynamicsJLD2Ext.jl` — JLD2 weak dependency; provides `save_solution` and `load_solution`.

### Tests (`test/`)

- `test_utils.jl` — Problem builders (`make_weber_problem()`, `make_coulomb_like_problem()`) and reference energy/velocity functions built from the independent two-body closed-form inverse; **must be included before other test files**
- `runtests.jl` — Entry point, includes all test files in order
- `test_canonical_weber.jl` — Independent physics checks of the canonical formulation: canonical momentum vs `∂L/∂v`, Legendre transform, both canonical equations against finite differences of `H`, recovery of the mechanical Weber force, Coulomb and zero-radial limits, translation/rotation invariance, and critical-radius behaviour. Expectations are derived from the Lagrangian, never from the routine under test
- Test files: `test_types.jl`, `test_hamiltonian_system.jl`, `test_initial_conditions.jl`, `test_solve.jl`, `test_statistics.jl`, `test_integration.jl`, `test_physics.jl`, `test_regularization.jl`

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

### Parameter layout

```
params = [m₁, ..., mₙ, q₁, ..., qₙ, c]                     # length 2N + 1
```

Compiled function signatures:

```
sys.dq_dt_compiled(out, q, p, t, params)
sys.dp_dt_compiled(out, q, p, t, params)
sys.hamiltonian_compiled(q, p, t, params)
```

### Canonical momentum

`p_i = ∂L/∂v_i = m_i v_i − Σ_{j≠i} (q_i q_j/c²) ṙ_ij (r_i − r_j)/r_ij²`.

`p/m` is **not** velocity. Use `physical_velocities(prob, q, p)` (or the
in-place `_weber_state!` on a `WeberWorkspace` in hot paths). Evaluating at a
like-charge pair's critical radius `ρ = q₁q₂/(μc²)` throws
`WeberCriticalRadiusError` — the canonical mass matrix is singular there.

### Regularization backends

Only two valid values for `RegularizationOptions.backend`:
- `:adaptive_cartesian` — Cartesian close-encounter substeps, works for all dimensions
- `:lifted_pair` — lifted square-root (1D), Levi-Civita (2D), or KS (3D) binary pair stepping

Multi-particle close clusters use chain mode, which is adaptive Cartesian over
the active component rather than analytic chain-coordinate regularization.
Neither backend analytically regularizes Weber's velocity-dependent force —
only the Coulomb/Kepler singularity. Under the exact canonical Hamiltonian the
Weber term modifies the pair's effective radial inertia `μ − q₁q₂/(rc²)`, which
the Levi-Civita/KS charts do not absorb. Measured against an independent
fine-step Cartesian reference, the lifted backend therefore no longer beats
plain Cartesian on state error in mild encounters; it does still improve energy
conservation. All three paths retain the integrator's second order.

`_extract_pair_relative_state!` in `regularization.jl` builds
`μ(p_i/m_i − p_j/m_j)`, which is the **canonical** relative momentum and is
correct as written — it is not a velocity.

### Collision bounce

- Implemented as the `CollisionBounce(radius)` callback; pass it to `solve` (or `init`) via the `callbacks` kwarg
- `RegularizedIntegrator` also accepts a `collision_bounce_radius` kwarg and synthesises a matching callback automatically when no `CollisionBounce` is supplied
- Under `RegularizedIntegrator`, the bounce radius is checked after regularized substeps as well as at macro-step boundaries
- Geometric and charge-sign agnostic: it reflects any pair inside the radius, including two-body unlike-charge pass-through when explicitly enabled
- Only valid for `L = 0` head-on, C0-continuable collisions/pass-through events; it is not generally energy-preserving for `N > 2`

### Makie animation extension

- Weak dependency is `Makie` (not `GLMakie`) — any backend triggers the extension
- `animate_weber(prob)` for live streaming, `animate_weber(sol)` for replay
- Compat: `Makie = "0.21, 0.22, 0.23, 0.24"`

### Immutable options pattern

`RegularizationOptions` is immutable — create it once, never mutate it. It is
algorithm-level and lives inside `RegularizedIntegrator`. Use the
`RegularizedIntegrator(base_alg; kwargs...)` constructor (its kwargs mirror
`RegularizationOptions`) rather than building the options struct by hand.

---

## EnergyStatistics fields

`en.statistics` has: `local_error_max`, `local_error_min`, `local_error_avg`, `global_error_ratio_max/min/avg`, `global_error_percent_max/min/avg`. There is **no** `local_error_percent_max` — use `local_error_max` for local error magnitude.
