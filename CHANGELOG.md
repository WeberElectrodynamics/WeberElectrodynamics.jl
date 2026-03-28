# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Zenodo DOI reference (`10.5281/zenodo.19239678`) added to the paper bibliography (`references.bib`) and cited at the point where the software implementation is introduced.

### Fixed
- Corrected Zenodo DOI badge URL in `README.md` (`badge/doi/` → `badge/DOI/`) so the badge renders correctly on GitHub.

## [0.2.2] - 2026-03-28

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
