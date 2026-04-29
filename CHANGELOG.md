# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2026-04-29

### Breaking changes

This release completes the layered **System → Problem → Algorithm → Callbacks** architectural refactor. The public API is reshaped end-to-end; no deprecation shims are provided.

- **Renamed top-level types.** `WeberSystem` → `HamiltonianSystem`, `WeberProblem` → `HamiltonianProblem`, `WeberSolution` → `HamiltonianSolution`, `WeberIntegrator` → `HamiltonianIntegrator`, `WeberAlgorithm` → `HamiltonianAlgorithm`.
- **`HamiltonianSystem` is now built from terms.** `HamiltonianSystem(n, dims)` still works (defaults to Weber + Zöllner). A new generic constructor `HamiltonianSystem(H, q, p; params, t)` and term builders `weber_term(...)` and `zollner_term(...)` let users assemble custom Hamiltonians. Each system exposes `system.terms::Vector{NamedTerm}` with per-term `pair_decomposition` closures driving statistics.
- **`HamiltonianProblem` fields dropped.** `masses`, `charges`, `c`, `kappas`, `regularization`, `zollner` are no longer fields on the problem. Construct with kwargs as before; access via exported accessors: `masses(prob)`, `charges(prob)`, `speed_of_light(prob)`, `kappas(prob)`, `regularization(prob)`, `zollner(prob)`, `params(prob)`, `n_particles(sys)`, `dims(sys)`.
- **Regularization is now an algorithm wrapper, not a problem option.** Replace `regularization=RegularizationOptions(...)` on `HamiltonianProblem` with `RegularizedIntegrator(SymmetricProjectionIntegrator(); r_on_factor=..., r_off_factor=..., backend=..., ...)` passed to `solve`.
- **Collision bounce is now a callback.** Replace `regularization=RegularizationOptions(collision_bounce_radius=r)` with `solve(prob, alg; callbacks=CollisionBounce(r))`.
- **Compiled EOM signature includes `t`.** `dq_dt_compiled(out, q, p, t, params)` and `dp_dt_compiled(out, q, p, t, params)` — direct callers must thread the time argument.
- **κ (Zöllner per-pair coupling) removed from `params`.** The parameter vector is now `[m₁…mₙ, q₁…qₙ, c]` (length `2N+1`); κ values live on `HamiltonianProblem.kappas` and are accessed via `kappas(prob)` (same accessor as before, now a direct field read). Compiled EOMs take κ as a separate positional argument: `dq_dt_compiled(out, q, p, t, params, kappas)`, `dp_dt_compiled(out, q, p, t, params, kappas)`, and `hamiltonian_compiled(q, p, t, params, kappas)`. New per-pair accessor `kappa(prob, i, j)` replaces manual `_pair_index` arithmetic in downstream code. Code that spliced κ onto `params` literals, sliced `params[(2N+2):end]`, or called the compiled EOMs directly must be updated.

### Added

- `RegularizedIntegrator{BaseAlg}` algorithm wrapper that delegates the inner step to the wrapped base algorithm via `_step_core!` / `_allocate_cache` hooks, preserving its symplectic projection while interleaving Levi-Civita / KS substeps near close encounters.
- `CollisionBounce(radius)` pre-step callback. Compose callbacks by passing a tuple (or any iterable of `HamiltonianCallback`) to the `callbacks=` kwarg of `solve` / `init`.
- Accessor API on `HamiltonianSystem` and `HamiltonianProblem` to insulate downstream code from struct layout.
- `NamedTerm{S}` with per-term compiled EOMs and optional `pair_decomposition(i, j, q, p, params, kappas)` for pair-wise statistics.
- `kappa(prob, i, j)` per-pair accessor, co-located with `_pair_index` in `hamiltonian_system.jl`.
- Symbolic builders: `weber_term`, `zollner_term`.
- Regression fixtures (`test/regression/fixtures/*.jld2`) with 1e-12 numerical-equivalence validation across every refactor phase.

### Migration guide

```julia
# Before
prob = WeberProblem(sys, (0.0, T), q0, p0;
    masses=ms, charges=qs, c=1.0,
    regularization=RegularizationOptions(enabled=true, backend=:lifted_pair,
                                          collision_bounce_radius=0.02))
sol = solve(prob, SymmetricProjectionIntegrator())

# After
prob = HamiltonianProblem(sys, (0.0, T), q0, p0; masses=ms, charges=qs, c=1.0)
alg  = RegularizedIntegrator(SymmetricProjectionIntegrator(); backend=:lifted_pair)
sol  = solve(prob, alg; callbacks=CollisionBounce(0.02))

# Field access
ms  = masses(prob)          # was: prob.masses
κs  = kappas(prob)          # was: prob.kappas (field access; slicing params[tail] no longer works)
κij = kappa(prob, 1, 2)     # per-pair — preferred over kappas(prob)[_pair_index(...)]
n   = n_particles(sys)      # was: sys.n_particles

# Direct compiled-EOM calls now take kappas as a separate positional arg
sys.dp_dt_compiled(out, q, p, t, params(prob), kappas(prob))
```

## [0.4.3] - 2026-04-18

### Changed

- `animate_weber` visual polish pass. Trajectories now use Wong-palette
  colors with an alpha-graded tail that fades older segments, rounded
  line joins, and stroked particle markers. 3D panels gain perspective
  (`perspectiveness = 0.3`) and a tight `viewmode = :fit`. Phase-space
  sidebar and info labels are restyled for readability. A shared
  `_weber_theme()` (LaTeX fonts, muted grid, clean spines) is applied
  via `with_theme` so all axes render consistently. No API change.

## [0.4.2] - 2026-04-17

### Added
- Research notebook `exhaust_nonzero_radial_velocity_ics.ipynb` exercising
  every IC recipe from `theory/NonZeroRadialVelocityBoundICs.md` across six
  scenarios: two-body mid-flight (2D outbound, 2D inbound, 3D tilted), hot
  binary with cold orbiter (3D), breathing square (N=4), and breathing
  hexagon (N=6). Each scenario builds physical ICs, applies the §6 forward
  map F to canonical momenta, runs the §7 verification checklist, integrates
  for ≥1 natural period, and produces trajectory + energy-error plots. The
  `verify_ic` helper amends the doc's §6.4 boxed identity with the exact
  analytic F-consistency check and also reports the Wesley-vs-Legendre H
  mismatch.
- Research notebook `scenario_3_makie_live.ipynb` demonstrating the
  redesigned `animate_weber` dashboard on the scenario-3 bound non-zero-ṙ
  IC (3D trajectory panel, phase-space sidebar, live energy-error readout).

### Changed
- `animate_weber` dashboard redesigned: the trajectory panel is now the
  dominant view, with the phase-space panel as a sidebar and a compact live
  energy-error readout below. The separate kinetic/potential,
  linear-momentum, and angular-momentum panels have been removed to reduce
  visual clutter and free up room for the trajectory view.
- The trajectory panel now uses `Axis3` when `prob.system.dims == 3`, giving
  a proper 3D view of the simulation alongside the existing 2D `Axis`
  rendering.
- Default `figure_size` reduced from `(1400, 900)` to `(1200, 800)` to
  better suit typical laptop displays.

The `animate_weber(prob)` / `animate_weber(sol)` public API (function names,
keyword arguments, return value) is unchanged.

## [0.4.1] - 2026-04-14

### Changed
- Plot labels in the `Plots` extension now use `LaTeXStrings` for axis labels
  and legend entries across all 11 plot functions (energy, pair energy,
  energy errors, pair forces, phase space, trajectories, momentum errors,
  and the four Zöllner plots). Titles remain prose. README figures
  regenerated from `examples/two_body_reference.ipynb`.
- `plot_energy` is now a single-panel plot of the energy components
  (`T`, `U`, `H`); the redundant relative-error panel has been removed.
- `plot_energy_errors` is now a two-panel plot (local error and relative
  energy error `|ΔE/E₀|`). The percentage scaling has been dropped in favour
  of the raw fraction, and the Hamiltonian-validation panel has been removed.
- `plot_momentum_errors` titles are now plain prose (`Linear Momentum Drift`
  / `Angular Momentum Drift`); the `‖ΔP‖` / `‖ΔL‖` symbols have moved from
  the panel titles to the y-axis labels where they belong.
- Absolute-value and norm labels across all plot functions
  (`plot_energy_errors`, `plot_momentum_errors`, `plot_pair_forces`,
  `plot_zollner_force_decomposition`) now use plain `|…|` / `||…||` bars
  instead of `\lvert…\rvert` / `\lVert…\rVert`, which the Plots.jl GR
  backend silently drops, leaving labels blank.
- README quickstart gained a third figure — `plot_energy_errors` — and the
  `examples/two_body_reference.ipynb` notebook now writes the three
  README figures to `examples/figures/` on every run.

### Added
- `LaTeXStrings` as a weak dependency of the `WeberElectrodynamicsPlotsExt`
  extension (compat `1`).

## [0.4.0] - 2026-04-13

### Breaking changes
- `plot_momentum` removed. Replaced by `plot_momentum_errors(data::MomentumData)`,
  which plots conservation errors as two stacked log-scale panels: linear
  drift `‖P(t) − P(0)‖` on top and angular drift `|Lz(t) − Lz(0)|` (2D) or
  `‖L(t) − L(0)‖` (3D) on the bottom. Each legend reports max absolute drift
  and, when the initial magnitude is nonzero, max relative drift
  `max_t ‖Δ·‖ / ‖·₀‖`. The old raw timeseries plot was visually flat for
  symplectic integrators in the COM frame and conveyed no integrator-quality
  signal. Migration: replace `plot_momentum(mom)` with
  `plot_momentum_errors(mom)`.

### Added
- Aqua.jl quality gate in the test suite — checks stale deps, unbound type
  parameters, undefined exports, compat bounds, and Project.toml formatting.
- CI matrix now includes `macos-latest` and `windows-latest` on Julia 1 in
  addition to the Ubuntu × {1.9, 1} coverage.
- `codecov.yml` with `project: auto` and `patch: 80%` targets, so unrelated
  PRs don't trip spurious coverage-drop failures.

### Changed
- `Project.toml` gains compat bounds for stdlibs (`LinearAlgebra`, `Printf`,
  `Random`, `Test`) and for the `Aqua` test extra, required by Aqua's
  `deps_compat` check.

### Chore
- Ignore `.claude/` (editor scratch directory) in `.gitignore`.

## [0.3.1] - 2026-04-12

### Fixed

- `RegularizationDiagnostics.used_backend` now correctly reports `:disabled` when
  regularization is enabled but no encounter is detected during the run (previously it
  incorrectly reported the effective backend even when no regularized step was ever taken).
- Removed dead 2D chain-mode LC lifts in `_step_regularized_chain!`: per-substep
  Levi-Civita lifts for each chain edge in 2D were overwriting the same buffers on every
  loop iteration and never being read, wasting two compiled-RHS calls per substep per edge.
- Removed always-false `active_count < 2` guard in `_detect_regularization_component!`:
  both anchor particles are unconditionally seeded into the BFS, so the count is always ≥ 2.

### Changed

- `theory/RegularizedIntegrationDesign.md`: corrected and completed the spec with six
  additions — adaptive substep count formula, multi-substep ABA composition details,
  fixed-anchor hysteresis note, `active_steps` diagnostics field, `collision_bounce_radius`
  option section, and clarified `used_backend = :disabled` for zero-encounter enabled runs.

## [0.3.0] - 2026-04-06

### Breaking changes

- `WeberProblem` no longer accepts flat `regularization_*` / `zollner_*` kwargs.
  Pass `RegularizationOptions(...)` and `ZollnerOptions(...)` structs directly:

  ```julia
  # Before (0.2.x)
  WeberProblem(...;
      regularization_enabled = true,
      regularization_backend = :adaptive_cartesian,
      zollner_enabled = true,
      zollner_a = 0.05,
  )

  # After (0.3.0)
  WeberProblem(...;
      regularization = RegularizationOptions(enabled = true, backend = :adaptive_cartesian),
      zollner        = ZollnerOptions(enabled = true, a = 0.05),
  )
  ```

  Both structs default to their disabled state (`RegularizationOptions()`,
  `ZollnerOptions()`), so problems with no regularization or Zöllner options
  require no changes.

### Changed

- `WeberProblem` constructor signature reduced from 20+ kwargs to 6 core kwargs
  plus the two options structs.
- Makie animation extension no longer performs an internal reverse-mapping from
  struct fields back to flat kwargs — it passes the structs through directly.

## [0.2.7] - 2026-04-05

### Fixed
- `_lc_lift!`: catastrophic cancellation in `r + x` when `x ≈ −r` (i.e. when the first LC pre-image coordinate is much smaller than the second). Direct subtraction `r + x` lost ~7 digits of precision; now uses `y²/(r − x)` (and symmetrically `y²/(r + x)` for the `r − x` branch), which is well-conditioned for all inputs. The bug surfaced as a flaky CI failure in the Levi-Civita round-trip identity test (`test_regularization.jl:215`, error ~3×10⁻⁷ vs tolerance 10⁻¹⁰).
- Levi-Civita round-trip test is now seeded (`Random.seed!(42)`) for deterministic reproduction.

## [0.2.6] - 2026-04-05

### Fixed
- Corrected `WeberProblem` docstring: `regularization_enabled` default was documented as `true` but has been `false` since v0.2.0.

### Changed
- Added inline code comments clarifying two known design choices in the regularization integrator: the frozen monitor (`r_eff`) in the Levi-Civita lifted-pair substep, and the one-pass (diagnostic-only) KS constraint projection in the 3D adaptive-Cartesian backend.
- Added note to `_compute_zollner_kappas` that neutral particles (`q = 0`) are treated as unlike any charged particle due to `sign(0.0)` semantics.
- Added a note to `docs/src/regularization.md` documenting that the 3D KS constraint is enforced by a one-pass projection (not iterative), and that `max_constraint_violation` tracks the residual.

## [0.2.5] - 2026-03-30

### Changed
- Language and citation improvements to the Computational Weber Electrodynamics paper (v1.2).

## [0.2.4] - 2026-03-29

### Changed

- Added DOIs to all journal articles (7) and Springer books (2) in the paper bibliography (`references.bib`); all verified via CrossRef.
- Fixed BibTeX title-casing: protected `{H}amiltonian` in Jayawardana2023 and Ohsawa2023 entries, `{WeberElectrodynamics.jl}` and `{Julia}` in the software entry.
- Corrected publisher of `assis-electric-force` from `C. Roy Keys Inc.` to `Apeiron` for consistency with all other Assis entries.
- Corrected year of `assis-weber-vol5` from 2024 to the confirmed publication year.
- Paper version bumped to 1.1.

## [0.2.3] - 2026-03-29

### Changed

- Factored source-code architecture and internal conventions out of `CLAUDE.md` into a new [Developer Guide](docs/src/developer-guide.md) docs page. `CLAUDE.md` now contains concise reminders with links; the full reference lives in the rendered docs.

## [0.2.2] - 2026-03-28

### Added
- Zenodo DOI reference (`10.5281/zenodo.19239678`) added to the paper bibliography (`references.bib`) and cited at the point where the software implementation is introduced.

### Fixed
- Corrected Zenodo DOI badge URL in `README.md` (`badge/doi/` → `badge/DOI/`) so the badge renders correctly on GitHub.

## [0.2.1] - 2026-03-28

### Added
- Four new tests covering previously untested scenarios: 3D `:lifted_pair` →
  `:adaptive_cartesian` fallback; collision bounce with `regularization_enabled = true`;
  Zöllner κ values respected during regularized substeps; `prob.params` tail equals
  `prob.kappas` (explicit params-vector layout verification).
- Diagnostics field reference table in `docs/src/regularization.md` explaining
  when `activation_count`, `min_encounter_distance`, `max_constraint_violation`,
  `backend_fallback_steps`, and `total_substeps` warrant attention.
- Clarification in `docs/src/zollner.md` that κ values are automatically included in
  the parameter vector for regularization sub-steps — no extra configuration needed.

### Fixed
- Added comment in `src/solve.jl` documenting the intentional sign-flip behaviour in
  the 1D regularization lift branch (explains why sign continuity is not enforced
  there, unlike the 2D/3D branches).

## [0.2.0] - 2026-03-27

### Changed
- `regularization_enabled` now defaults to `false` (was `true`). The core
  symplectic integrator runs unregularized by default. Pass
  `regularization_enabled = true` to `WeberProblem` to opt in to Levi-Civita /
  KS regularization for close encounters.
- `RegularizationOptions(enabled = ...)` default flipped to `false` accordingly.

### Docs
- New page **Regularization** (`docs/src/regularization.md`): usage guide,
  backend selection, hysteresis parameters, chain mode, and collision bounce.
- New page **Zöllner Extension** (`docs/src/zollner.md`): theory background,
  usage, and Zöllner-specific plot functions. Marked *Research / Experimental*.
- Doc navigation reorganised into Core → Advanced (Regularization) →
  Research (Zöllner Extension) → Theory tiers.
- `RegularizationOptions`, `RegularizationDiagnostics`, and `ZollnerOptions`
  API docs moved from `api/problem.md` to their dedicated feature pages.
- Added "Optional features" section to Quick Start linking both feature pages.

### Refactored
- Added `# Zöllner Extension` and `# Regularization` section comment headers
  in `src/types.jl`, `src/statistics/energy.jl`, and `src/statistics/forces.jl`
  to make the tier boundaries visible at a glance in the source.

## [0.1.1] - 2026-03-26

### Docs

- Added Documenter.jl scaffold (`docs/`), API reference, quick-start guide, and theory page.
- Added GitHub Actions workflow for automatic doc deployment to GitHub Pages (`Docs.yml`).
- Added CI badges, Docs badges, Coverage (Codecov), PkgEval, Julia compat, and license badges to README.
- Added Codecov coverage upload to `CI.yml`.
- Moved research documents (`theory/`, `exploratory/`, `sub_critical_weber_research/`) from `docs/` to `research/`.

## [0.1.0] - 2026-03-22

### Added
- Initial registered release.
- `WeberSystem`: symbolic Weber Hamiltonian construction via Symbolics.jl, compiling
  `dq_dt`, `dp_dt`, and `hamiltonian` functions.
- `WeberProblem`, `WeberSolution`, `WeberIntegrator`, `SymmetricProjectionIntegrator` types.
- CommonSolve.jl interface (`solve`, `init`, `step!`, `solve!`).
- Symplectic Strang-splitting integrator with symmetric projection via fixed-point iteration
  on Lagrange multipliers for Weber's non-separable Hamiltonian.
- `RegularizationOptions`: close-encounter regularization with Levi-Civita (2D) and
  KS quaternion (3D) backends, adaptive hysteresis switching, chain encounter fallback.
- Collision bounce for sub-critical head-on (ℓ=0) like-charge collisions
  (`regularization_collision_bounce_radius` kwarg).
- `ZollnerOptions`: Zöllner electrogravitational mismatch extension (κ per pair,
  `zollner_enabled` / `zollner_a` kwargs on `WeberProblem`).
- Statistics module: `TrajectoryData`, `EnergyData`, `PairEnergyData`, `PairForceData`,
  `MomentumData` with `compute_*` functions and optional `stride` downsampling.
- Plots.jl weak-dependency extension: `plot_trajectories`, `plot_energy`, `plot_pair_energy`,
  `plot_energy_errors`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum`,
  `plot_zollner_energy`, `plot_zollner_force_residual`, `plot_weber_vs_zollner`,
  `plot_zollner_phase_space`.
- Makie weak-dependency extension (any backend): `animate_weber` with streaming
  (`animate_weber(prob)`) and replay (`animate_weber(sol)`) modes, interactive
  dashboard with trajectory, energy, momentum, angular momentum, and phase-space panels.
