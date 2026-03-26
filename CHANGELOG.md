# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
