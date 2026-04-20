# WeberElectrodynamics Architectural Refactor — Phased Plan

## Context

The package today hardcodes the Weber Hamiltonian in `WeberSystem(n, dims)`, bundles Zöllner coupling factors (κ) into the `params` vector, and embeds regularization + collision bounce into the core `step!` function. This blocks extensibility for future work: time-dependent Hamiltonians, Lagrange-multiplier constraints, user-provided custom Hamiltonians, promoting physical parameters (mass, charge) to generalized coordinates, and cleanly layered feature extensions.

This plan refactors the package into four explicit layers — **System → Problem → Algorithm → Callbacks** — with a stable integrator interface, symbolic Hamiltonian builder, term composition for physics extensions, and algorithm-wrapper composition for step-substituting features (regularization). Backwards compatibility is deliberately broken once; no deprecation shims. The refactor prepares for but does not implement Lagrange multipliers or non-Weber physics; those are follow-ups that will land cheaply once the architecture is in place.

Pre-refactor research (see below) confirmed that: ModelingToolkit has no usable `HamiltonianSystem` in current versions; SciML callbacks cannot substitute a step (forcing regularization into algorithm composition, not callbacks); and Symbolics.jl has a known compile-time cliff around 100+ variables that warrants a pre-refactor benchmark.

## Confirmed design decisions

1. **Symbolic Hamiltonian builder**: `HamiltonianSystem(H, q, p; params, t)` with first-class `terms::Vector{NamedTerm}` preserving the decomposition (not a black box).
2. **Zöllner as additive term**: `H = weber_term(...) + zollner_term(...; a=...)`. `ZollnerOptions` removed. κ values leave the core params layout and become generic term-owned params.
3. **Regularization as algorithm wrapper** (NOT callback): `RegularizedIntegrator(SymmetricProjection(...); options...)`. SciML callbacks cannot substitute steps; this is a technical constraint.
4. **Callback system for observer hooks**: `CallbackSet(CollisionBounce(r), EnergyObserver(stride), ...)`. Pre-step-only hooks like collision bounce fit cleanly.
5. **Stable integrator interface**: `step!(cache, prob, alg, t, dt)` dispatches on `alg` type. New algorithms (ImplicitMidpoint, GaussLegendre, …) plug in via new `step!` methods without editing existing code.
6. **`t` added to EOM contract**: `dq_dt_compiled(out, q, p, t, params)`. Current H is autonomous; `t` is a free symbolic variable with zero derivative until someone writes a time-dependent term.
7. **Independent of SciMLBase**: keep the `DynamicalODEProblem` *shape* (f1/f2, in-place) without the type hierarchy. `CommonSolve` remains the dispatch protocol.
8. **`WeberProblem` → `HamiltonianProblem`**: no `zollner=` / `regularization=` kwargs (those move to algorithm / callbacks).

## Resolved defaults (can be revised)

- **Kinetic term convention**: `weber_term(...)` includes kinetic + potential (matches the sketch `H = weber_term(q,p,params) + zollner_term(...)`). A separate `kinetic_term(q, p; masses)` builder exists for future non-Weber custom Hamiltonians. Alternative (split kinetic out of `weber_term`) is viable but changes the user-facing sketch.
- **Regression fixtures committed**: `test/regression/fixtures/*.jld2` checked into git so CI can validate equivalence without a local regeneration step.
- **Extension metadata access**: add accessor functions `n_particles(sys)`, `dims(sys)`, `masses(prob)`, etc. Backed by a `metadata::NamedTuple` on `HamiltonianSystem` to keep extensions stable across future Hamiltonian shapes.

## Target architecture

```
User code
    │
    │   H = weber_term(q,p,params) + zollner_term(q,p,params; a=0.1)
    │   sys = HamiltonianSystem(H, q, p; params, t)
    │   prob = HamiltonianProblem(sys, q0, p0, (0.,T); param_values, dt)
    │   alg = RegularizedIntegrator(SymmetricProjection(relaxation=0.25);
    │                               r_on_factor=0.15, backend=:lifted_pair)
    │   cbs = CallbackSet(CollisionBounce(0.02), EnergyObserver(10))
    │   sol = solve(prob, alg; callbacks=cbs)
    │
    ▼
┌────────────────────────────────────────────────────────────────────┐
│ Layer 1: System (src/hamiltonian/)                                 │
│   HamiltonianSystem{…} { terms::Vector{NamedTerm}, compiled EOMs, │
│                          metadata::NamedTuple }                   │
│   Builders: weber_term, zollner_term, kinetic_term, coulomb_term  │
└────────────────────────────────────────────────────────────────────┘
    │  compiled dq_dt(out,q,p,t,params), dp_dt(out,q,p,t,params)
    ▼
┌────────────────────────────────────────────────────────────────────┐
│ Layer 2: Problem (src/problem.jl)                                  │
│   HamiltonianProblem { system, tspan, q0, p0, params, dt, tol, .} │
│   No regularization, no Zöllner options — those are elsewhere.    │
└────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌────────────────────────────────────────────────────────────────────┐
│ Layer 3: Algorithm (src/integrators/)                              │
│   abstract HamiltonianAlgorithm                                   │
│   step!(cache, prob, alg, t, dt) :: Bool                          │
│   SymmetricProjection(relaxation) — base algorithm                │
│   RegularizedIntegrator(base_alg; options) — algorithm wrapper    │
│     owns RegularizationCache; calls into base_alg.projection_kernel│
└────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌────────────────────────────────────────────────────────────────────┐
│ Layer 4: Callbacks (src/callbacks.jl)                              │
│   CollisionBounce (pre-step), EnergyObserver (post-step), …       │
│   Shaped like DiffEq DiscreteCallback but independent impl.       │
└────────────────────────────────────────────────────────────────────┘
```

## Phased implementation

Each phase is a single branch, reviewable in isolation, merged to `main` only when all acceptance criteria pass. CI stays green on each merge. Phases 0–3 are sequential; 4–7 parallelizable once 3 lands.

### Phase 0 — Symbolics compile benchmark spike (GATE)

**Branch**: `bench/symbolics-compile-spike` → lands on `main` standalone.

**Why first**: if generic `HamiltonianSystem` build at n=8–12 hits the known Symbolics compile cliff, we must adopt `cse=true` / `ShardedForm` / `MultithreadedForm` in the core compile pipeline from day one rather than discovering it mid-refactor.

**Deliverables**:
- New [benchmarks/bench_symbolics.jl](benchmarks/bench_symbolics.jl) measuring: symbolic build time, `Symbolics.derivative` gradient cost, `build_function` compile time, first-call latency, expression byte size, eval throughput. Sweep n ∈ {3, 5, 8, 12}, dims ∈ {2, 3}.
- Four build variants: current Weber symbolic path, current + `cse=true`, a *candidate* generic build, generic + `cse=true` and/or sharded form if available.
- Results written to [benchmarks/results_phase0.md](benchmarks/results_phase0.md) (committed).

**GATE criteria**:
- Generic build at n=8, dims=3 compiles in ≤ 2× current Weber time, OR
- `cse=true` / sharded form closes the gap to ≤ 2× (adopt as core default), OR
- Document a deferred/lazy per-term compilation fallback.

**Size**: M. **Risk**: low — benchmark only, no API changes.

### Phase 1 — `HamiltonianSystem` + term structure + `t` signature

**Branch**: `refactor/hamiltonian-system`

**This is the big break.** Old public API (`WeberSystem`, `WeberProblem`, `ZollnerOptions`) is deleted. Most tests will be red until Phases 2–5 land.

**Files created**:
- [src/hamiltonian/system.jl](src/hamiltonian/system.jl) — `HamiltonianSystem{…}` with `terms::Vector{NamedTerm}`, aggregate compiled EOMs, `metadata::NamedTuple` for extensions.
- [src/hamiltonian/terms.jl](src/hamiltonian/terms.jl) — `NamedTerm{S}` with per-term `H_compiled`, `dq`, `dp`, optional `pair_decomposition(i, j, q, p, params)` for statistics.
- [src/hamiltonian/builders/weber.jl](src/hamiltonian/builders/weber.jl) — `weber_term(q, p, params; masses, charges, c, n_particles, dims)`. Reuses the existing body of `_build_weber_hamiltonian` (weber_system.jl:79-131) and `_pair_index` (weber_system.jl:75). Attaches a `pair_decomposition` closure implementing the decomposition currently hardcoded in `compute_pair_weber_components` (energy.jl:110-161).
- [src/hamiltonian/builders/zollner.jl](src/hamiltonian/builders/zollner.jl) — `zollner_term(q, p, params; charges, a, c)` or `zollner_term(...; kappas)`. Migrates `_compute_zollner_kappas` (types.jl:78-97).
- [src/hamiltonian/builders/kinetic.jl](src/hamiltonian/builders/kinetic.jl) — `kinetic_term` for future non-Weber H. Not used by `weber_term` itself (which keeps its current kinetic+potential shape).
- [src/hamiltonian/builders/coulomb.jl](src/hamiltonian/builders/coulomb.jl) — optional; pure Coulomb for tests.
- [src/problem.jl](src/problem.jl) — `HamiltonianProblem{…}` with only system-agnostic state: `system`, `tspan`, `q_initial`, `p_initial`, `params::Vector{Float64}`, `dt`, `convergence_tolerance`, `maximum_iterations`. Drops `masses`, `charges`, `c`, `kappas`, `regularization`, `zollner` as fields.

**Files modified**:
- [src/WeberElectrodynamics.jl](src/WeberElectrodynamics.jl) — include graph + exports (remove `WeberSystem`, `WeberProblem`, `ZollnerOptions`; add `HamiltonianSystem`, `HamiltonianProblem`, `weber_term`, `zollner_term`, `kinetic_term`, `get_term`, `has_term`, `term_names`, `n_particles`, `dims`).
- [src/types.jl](src/types.jl) — delete `WeberProblem` (466-533), `ZollnerOptions` (61-71), `_compute_zollner_kappas` moved. `SymmetricProjectionBuffers` (621-693) and `RegularizationBuffers` (271-424) remain until Phase 2/3.
- [src/solve.jl](src/solve.jl) — signature-only change: every call site of `dq_dt_compiled` / `dp_dt_compiled` gains `t` argument (~20 sites: lines 31, 37, 43, 103, 143, 168, 193, 214, 357-360, 635-636, 695-696, and the same in regularization paths).
- [src/regularization.jl](src/regularization.jl) — signature-only `t` propagation.
- [src/statistics/*.jl](src/statistics/) — signature-only `t` propagation. (Decomposition rewrite happens in Phase 4.)

**Files deleted**:
- [src/weber_system.jl](src/weber_system.jl) — logic migrated into `src/hamiltonian/builders/weber.jl`.

**Acceptance**:
- Package loads. `HamiltonianSystem`, `HamiltonianProblem`, `weber_term`, `zollner_term`, `get_term` exist.
- Single new minimal smoke test `test/test_phase1_smoke.jl`: builds a 2-body Weber problem via the new API, solves with the still-embedded Strang path, matches the Phase-0 captured trajectory within 1e-12.
- Old tests are **expected broken** — they are rewritten in Phase 5.

**Size**: XL. **Risks**: compile cliff (mitigated by Phase 0); any of ~20 `t` call sites missed; extension field access breakage (partially handled by `metadata` accessors).

### Phase 2 — Stable integrator interface + algorithm dispatch

**Branch**: `refactor/integrator-interface`

**Files created**:
- [src/integrators/api.jl](src/integrators/api.jl) — `abstract type HamiltonianAlgorithm end`, `step!(cache, prob, alg, t, dt) :: Bool`, `allocate_cache(alg, prob)`, `init(prob, alg)`. Traits: `supports_projection(::Alg)`, `projection_kernel(cache, alg)` returns a callable kernel `(Z, Ẑ, μ, dt) → residual` used by wrapping algorithms.
- [src/integrators/symmetric_projection.jl](src/integrators/symmetric_projection.jl) — `SymmetricProjection(; relaxation=0.25)`. Reuses existing `_projected_cartesian_step!` body (solve.jl:92-233), `compute_constraint_residual!` (solve.jl:52-90), `strang_splitting_flow!` (solve.jl:4-50) unchanged modulo `t` threading. Owns `SymmetricProjectionBuffers` (moved out of types.jl:621-693).
- [src/integrators/cache.jl](src/integrators/cache.jl) — cache allocation helpers.

**Files modified**:
- [src/solve.jl](src/solve.jl) — `CommonSolve.step!(integrator)` (1118) simplified to dispatch on `alg` via `step!(cache, prob, alg, t, dt)`. Collision bounce inlining (1140-1148) and regularization branching (1150-1154) removed from core `step!` — they re-enter via Phase 3 wrappers/callbacks.
- [src/types.jl](src/types.jl) — remove `SymmetricProjectionIntegrator` type (31-38), `SymmetricProjectionBuffers` (621-693). Rename `WeberIntegrator` → `HamiltonianIntegrator` (718-731) with generic `cache::Any`, `alg`, `callbacks` (empty set by default).

**Acceptance**:
- `solve(prob, SymmetricProjection())` reproduces Phase-0 captured *unregularized* trajectory to 1e-12.
- `step!(cache, prob, alg, t, dt)` is the sole hot-path method. No regularization or bounce branching in core `solve!` loop.
- New algorithm can be added without touching existing code (demonstrate with a stub `ExplicitEuler` that the tests reject for non-separable H — proof of dispatch).

**Size**: L. **Risks**: `projection_kernel` trait surface may prove insufficient for Phase 3's needs; expect to refine.

### Phase 3 — Algorithm wrappers, callbacks, regression validation

**Branch**: `refactor/wrappers-callbacks-regression`

**This phase establishes numerical equivalence against `main`.** Without it, the refactor cannot be verified.

**Files created**:
- [src/integrators/regularized.jl](src/integrators/regularized.jl) — `RegularizedIntegrator{BaseAlg}(base_alg; r_on_factor, r_off_factor, backend, max_substeps, constraint_tolerance, g_floor, chain_enabled, warn_on_fallback)`. Implements `step!` by reimplementing `_step_regularized_dispatch!` logic (solve.jl:964-997) using `projection_kernel(base_cache, base_alg)` rather than hardcoding `_projected_cartesian_step!`. Owns `RegularizedCache` wrapping base cache + `RegularizationBuffers` (migrated from types.jl:271-424).
- [src/callbacks.jl](src/callbacks.jl) — `abstract type HamiltonianCallback end`, `DiscreteCallback(condition, affect!; save_positions)`, `CallbackSet(cbs...)`. Pre-step and post-step dispatch points inside `solve!`.
- [src/callbacks/collision_bounce.jl](src/callbacks/collision_bounce.jl) — `CollisionBounce(radius)` pre-step callback. Migrates `_apply_collision_bounces!` (solve.jl:1088-1105) and `_reflect_pair!` (solve.jl:558-581) verbatim.
- [src/callbacks/observers.jl](src/callbacks/observers.jl) — `EnergyObserver(stride)` (post-step), `RegularizationDiagnosticsObserver` (replaces the embedded `RegularizationDiagnostics` struct updates scattered through solve.jl).
- [test/regression/capture.jl](test/regression/capture.jl) — run once against `main` BEFORE Phase 1 merges (via a throwaway checkout). Produces 4 JLD2 fixtures in `test/regression/fixtures/`:
  - `twobody_ellipse.jld2` — 2-body unregularized, full (t, q, p).
  - `threebody_chain.jld2` — 3-body unregularized.
  - `close_approach_lifted.jld2` — 2-body close approach with `backend=:lifted_pair`.
  - `zollner_offmatch.jld2` — 3-body Zöllner with a ≠ 0.
- [test/regression/validate.jl](test/regression/validate.jl) — runs after every phase from Phase 1 onward. Builds the equivalent new-API problem, runs, asserts `maximum(abs, q_new .- q_old) < 1e-12` componentwise.

**Files modified**:
- [src/solve.jl](src/solve.jl) — add pre-step / post-step callback dispatch around `step!`. All regularization branches deleted; they live in `RegularizedIntegrator` now.
- [src/regularization.jl](src/regularization.jl) — `_step_regularized_pair_lifted_2d!`, `_step_regularized_pair_adaptive!`, `_step_regularized_chain!` moved into `src/integrators/regularized.jl` as methods of `RegularizedIntegrator`'s `step!`. Pure helpers (`_lc_lift!`, `_lc_project!`, `_ks_*`) remain but are relocated to `src/integrators/regularization_helpers.jl`.

**Acceptance**:
- All 4 regression fixtures pass at 1e-12.
- `solve(prob, RegularizedIntegrator(SymmetricProjection(); ...))` with callback set mixing `CollisionBounce` + `EnergyObserver` produces identical regularization diagnostics as current `main`.
- Projection-kernel trait proves sufficient; if not, document what had to widen.

**Size**: XL. **Risks**: projection-kernel trait too narrow → may need to expose `compute_constraint_residual!` directly. `RegularizationBuffers` migration into wrapper cache is plumbing-heavy.

### Phase 4 — Statistics re-platform (parallel with 5–7)

**Branch**: `refactor/statistics-term-queries`

**Files modified**:
- [src/statistics/energy.jl](src/statistics/energy.jl:110-161) — `compute_pair_weber_components` rewritten to call `term.pair_decomposition(i, j, q, p, params)` on the `:weber` and `:zollner` named terms. Caller at line 309 updated.
- [src/statistics/forces.jl](src/statistics/forces.jl:124) — same treatment.
- [src/statistics/momentum.jl](src/statistics/momentum.jl), [src/statistics/trajectories.jl](src/statistics/trajectories.jl) — switch to `n_particles(sys)`, `dims(sys)`, `masses(prob)` accessors.
- [ext/WeberElectrodynamicsPlotsExt.jl](ext/WeberElectrodynamicsPlotsExt.jl:962-970) — replace direct `sol1.prob.system.n_particles`, `.dims`, `.masses`, `.charges`, `.c` reads with accessor-based access. Zöllner-specific plots (`plot_zollner_energy`, `plot_zollner_force_residual`, `plot_weber_vs_zollner`, `plot_zollner_phase_space`) branch on `has_term(sys, :zollner)`.
- [ext/WeberElectrodynamicsMakieExt.jl](ext/WeberElectrodynamicsMakieExt.jl) — same accessor changes; `animate_weber` otherwise unchanged.

**Acceptance**: `plot_pair_energy`, `plot_zollner_energy`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum_errors`, `animate_weber` all work on a new-API solution and produce visually identical output to captured PNG snapshots.

**Size**: L.

### Phase 5 — Tests, examples, docs (parallel with 4, 6, 7)

**Branch**: `refactor/tests-docs-examples`

**Files modified**:
- [test/test_utils.jl](test/test_utils.jl) — `make_weber_problem`, `make_coulomb_like_problem` rewritten to use `HamiltonianSystem` + `weber_term` + `HamiltonianProblem`.
- [test/test_types.jl](test/test_types.jl) — rewritten for new types; tests for `HamiltonianSystem`, `HamiltonianProblem`, `RegularizedIntegrator` construction, `NamedTerm` queries.
- [test/test_weber_system.jl](test/test_weber_system.jl) — `WeberSystem` tests reframed as `weber_term` + `HamiltonianSystem` tests. Direct `dq_dt_compiled(out, q, p, params)` calls updated to include `t`.
- [test/test_zollner.jl](test/test_zollner.jl) — reframed around `zollner_term` builder. No more `ZollnerOptions`. Tests verify κ computation, term addition, pair decomposition query.
- [test/test_regularization.jl](test/test_regularization.jl) — reframed around `RegularizedIntegrator`. All regularization diagnostics now come from `RegularizationDiagnosticsObserver` callback.
- [test/test_physics.jl](test/test_physics.jl), [test/test_solve.jl](test/test_solve.jl), [test/test_integration.jl](test/test_integration.jl), [test/test_statistics.jl](test/test_statistics.jl) — constructor updates; any raw `params` manipulations rewritten.
- [examples/two_body_reference.ipynb](examples/two_body_reference.ipynb) — rewire setup cells to new API.
- [docs/src/api/system.md](docs/src/api/system.md), [docs/src/api/problem.md](docs/src/api/problem.md), [docs/src/quickstart.md](docs/src/quickstart.md), [docs/src/regularization.md](docs/src/regularization.md), [docs/src/zollner.md](docs/src/zollner.md), [docs/src/internals.md](docs/src/internals.md) — rewritten for new API. `internals.md`'s "Params vector layout" and "Regularization backends" sections replaced.

**Acceptance**:
- `julia --project=. -e 'using Pkg; Pkg.test()'` all green.
- `make.jl` builds docs cleanly.
- Example notebook runs end-to-end.

**Size**: L.

### Phase 6 — `_research/` sweep (parallel with 4, 5, 7)

**Branch**: `refactor/research-sweep` (may split into per-subdirectory sub-branches if review load warrants).

**Scope**: 46 `.jl` files + 15 `.ipynb` notebooks under [_research/](_research/). Most are mechanical `WeberSystem(...)` → `HamiltonianSystem(..., weber_term(...))` and `WeberProblem(...)` → `HamiltonianProblem(...)` replacements. Some scripts directly use `dq_dt_compiled`/`dp_dt_compiled` or raw `params` construction — these need the `t` signature update and the kappa-layout change.

**Strategy**:
- Scripted pass 1: regex-based search/replace for mechanical patterns (constructor renames, kwarg drops).
- Manual pass 2: per-notebook runthroughs for scripts that manipulate `params` or call compiled EOMs directly.
- Representative high-impact files: [_research/notebooks/two_body_regularized/common.jl](_research/notebooks/two_body_regularized/common.jl), [_research/homology/*/*.jl](_research/homology/), [_research/FourBodyTwoPlusTwoMinus/*/*.jl](_research/FourBodyTwoPlusTwoMinus/).

**Acceptance**: every notebook executes end-to-end; every `.jl` script at least parses and imports cleanly (some exploratory scripts may not have reproducible runtime state and are marked as sandbox-only in their header).

**Size**: XL (breadth, not depth).

### Phase 7 — Release prep (parallel with 4, 5, 6)

**Branch**: `refactor/release-prep`

**Files modified**:
- [RELEASENOTES.md](RELEASENOTES.md) — breaking-change entries under `### Breaking changes` (word "breaking" needed for AutoMerge), migration guide section.
- [CHANGELOG.md](CHANGELOG.md) — populated automatically by `release.sh` from RELEASENOTES.
- [CLAUDE.md](CLAUDE.md) — "Params vector layout" section rewritten; "Regularization backends" section reframed as algorithm wrapper; Zöllner mentions updated; code examples in the repo-structure section updated.
- [Project.toml](Project.toml) — minor bump (pre-1.0 breaking convention: 0.X → 0.(X+1).0).
- `/Users/mac/.claude/projects/-Users-mac-dev-Weber-WeberElectrodynamics/memory/MEMORY.md` entries — update "Params Vector Structure", "Zöllner Feature", "Test Conventions" (params now lack kappas), "Regularization Backends" (now an algorithm wrapper), "Collision Bounce" (now a callback). `user_name_privacy.md` and `feedback_branch_before_edits.md` remain unchanged.

**Acceptance**: `./release.sh minor` dry-run shows clean version bump + release notes; CLAUDE.md documents the new API.

**Size**: M.

## Parallelization map

```
Phase 0 (gate) ──► Phase 1 ──► Phase 2 ──► Phase 3 (equivalence proven)
                                                 ├──► Phase 4 (stats)
                                                 ├──► Phase 5 (tests/docs/examples)
                                                 ├──► Phase 6 (_research/)
                                                 └──► Phase 7 (release prep)
```

Phase 4 and 5 have a mild coupling: `test_statistics.jl` in Phase 5 depends on Phase 4 shipping first. Land Phase 4 one commit ahead of `test_statistics.jl` portion of Phase 5, or merge them together.

Phases 6 and 7 are fully independent of 4 and 5 and of each other.

## Numerical equivalence strategy

Before Phase 1 merges, check out `main` in a worktree and run `test/regression/capture.jl` to produce four JLD2 fixtures covering unregularized 2-body, unregularized 3-body, regularized close-approach, and Zöllner off-match. Commit fixtures to `test/regression/fixtures/`.

After every phase from Phase 1 onward, `test/regression/validate.jl` rebuilds each fixture problem using the current (new-API) code, runs it, and asserts `maximum(abs, q_new .- q_old) < 1e-12` and similarly for `p`. Any drift beyond that threshold blocks the phase merge.

This is the single most important invariant in the plan: the refactor must prove numerical equivalence phase-by-phase, so we catch drift immediately rather than at the end.

## Out of scope (explicitly)

- **Actual Lagrange-multiplier / holonomic constraint implementation.** This refactor *prepares* for it via stable interfaces but does not implement it.
- **Non-Weber physics beyond the Zöllner extension.** Future `H = custom_term(...)` workflows will work on day one of Phase 1; exercising them with real physics is out of scope.
- **SciMLBase protocol subtyping.** Decided against it — we mimic the `DynamicalODEProblem` shape but not the type hierarchy.
- **Backwards-compatibility shims / deprecations.** No `@deprecate` calls. Break once, bump minor.
- **[theory/](theory/) and [papers/](papers/) updates.** Confirmed to have zero Julia API surface.

## Critical files (refactor pivots)

- [src/weber_system.jl:75](src/weber_system.jl#L75) — `_pair_index(i, j, n)` survives refactor, migrates to `src/hamiltonian/builders/weber.jl`.
- [src/weber_system.jl:79-131](src/weber_system.jl#L79-L131) — `_build_weber_hamiltonian` body migrates into `weber_term`.
- [src/weber_system.jl:150](src/weber_system.jl#L150) — current `WeberSystem` constructor; replaced by generic `HamiltonianSystem`.
- [src/types.jl:78-97](src/types.jl#L78-L97) — `_compute_zollner_kappas` migrates into `zollner_term`.
- [src/types.jl:271-424](src/types.jl#L271-L424) — `RegularizationBuffers` migrates into `RegularizedCache` (Phase 3).
- [src/types.jl:466-533](src/types.jl#L466-L533) — `WeberProblem` replaced by `HamiltonianProblem` (Phase 1).
- [src/types.jl:621-693](src/types.jl#L621-L693) — `SymmetricProjectionBuffers` migrates into `SymmetricProjection`'s cache (Phase 2).
- [src/solve.jl:4-50](src/solve.jl#L4-L50) — `strang_splitting_flow!` survives with `t` threading.
- [src/solve.jl:52-90](src/solve.jl#L52-L90) — `compute_constraint_residual!` survives with `t` threading.
- [src/solve.jl:92-233](src/solve.jl#L92-L233) — `_projected_cartesian_step!` becomes the body of `SymmetricProjection`'s `step!`.
- [src/solve.jl:964-997](src/solve.jl#L964-L997) — `_step_regularized_dispatch!` becomes the body of `RegularizedIntegrator`'s `step!`.
- [src/solve.jl:1088-1105](src/solve.jl#L1088-L1105) — `_apply_collision_bounces!` becomes the body of `CollisionBounce` callback's `affect!`.
- [src/solve.jl:1118](src/solve.jl#L1118) — `CommonSolve.step!` is the point where callbacks + algorithm dispatch re-enter.
- [src/statistics/energy.jl:110-161](src/statistics/energy.jl#L110-L161) — `compute_pair_weber_components` becomes `weber_term.pair_decomposition(i, j, q, p, params)` (Phase 4).

## Risks & rollback

- **Phase 0 fails the gate.** Pause refactor, revisit whether we need a lazier compile strategy (e.g. per-term JIT) before Phase 1. This is the single biggest "could force a design revision" event.
- **Projection-kernel trait proves insufficient in Phase 3.** Mitigation: widen to expose the full `compute_constraint_residual!` closure; regularization wrapper then owns the residual loop itself.
- **`t` propagation missed in a call site.** Mitigated by the regression fixtures: any missed site will produce numerically drifted output and fail Phase 1 validation.
- **Extension (Plots, Makie) breakage.** Mitigated by accessor functions + `metadata::NamedTuple`. If more fields are needed, add accessors rather than exposing struct internals.
- **Rollback**: each phase is a single mergeable branch off the previous merge commit. A `git revert` of the merge commit restores the prior state. The regression fixtures prove equivalence at each merge point so a revert is always numerically safe.

## Verification

After each phase merge:
1. `julia --project=. -e 'using Pkg; Pkg.test()'` — full suite green (post-Phase-5; subset green earlier).
2. `julia --project=. test/regression/validate.jl` — all four fixtures within 1e-12 (from Phase 1 onward).
3. For phases 4+: `julia --project=docs/ docs/make.jl` — docs build clean.
4. For Phase 5: the example notebook runs end-to-end.
5. For Phase 6: every `.ipynb` under `_research/` runs end-to-end; every `.jl` at least parses + imports.
6. For Phase 7: `./release.sh minor` dry-run prints a clean preview; no git state is modified.

End-to-end final validation:
- Capture a fresh trajectory from a canonical 3-body Zöllner off-match problem using the new API.
- Compare against the Phase 0 `main`-captured fixture — max pointwise error < 1e-12.
- Run the example notebook, inspect all plots and animation output.
- Document migration in RELEASENOTES.md with a code-snippet diff per major change.
