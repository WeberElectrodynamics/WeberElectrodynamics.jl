# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [5.1.0] - 2026-03-02

### Added
- Animation viewer interactive sliders: Trail length, Display window, and Speed
  (log-linear 1--1000 steps/frame) — adjustable during playback.
- Display windowing: energy, momentum, and angular momentum plots show only the
  most recent N timesteps (controlled by Window slider), with `track_x` autoscaling
  that follows the advancing time window instead of one-way expansion.
- `_trim_to_window` helper for rolling display of time-series data.
- `_log_linear_range` helper generating 1,2,...,9,10,20,...,90,100,...,1000 speed values.
- Y-axis one-way limits reset automatically when the display window size changes.

### Changed
- `AnimationState.tail_length`, `display_window`, `compute_batch` are now
  `Observable{Int}` for reactive slider binding (was plain `Int`).
- `_start_animation!` reads `compute_batch` from state each frame instead of
  capturing a fixed value at startup.
- Default `compute_batch` changed from 10 to 1 for both streaming and replay modes
  (user adjusts via Speed slider).
- Bottom-row layout reorganized: controls and sliders on the left, phase-space
  selector menu and error display (2x3 grid, larger font) on the right.
- Error labels increased from fontsize 9 to 12 for readability.

## [5.0.1] - 2026-03-02

### Fixed
- Animation viewer now opens a native GLMakie window from Jupyter notebooks
  instead of rendering a static inline PNG (`Makie.inline!(false)` + explicit
  `display`).
- `xlims`/`ylims` getter calls replaced with `_get_xlims`/`_get_ylims` helpers
  reading `ax.finallimits[]`; the getters are not exported in Makie 0.24.
- Unicode emoji in play/pause/reset button labels (`⏸`, `▶`, `↺`) replaced
  with ASCII text; TeX Gyre Heros Makie font cannot render these glyphs.
- One-way auto-scaling: axis limits only expand, never shrink, preventing
  jitter during animation playback.
- Reset button now calls `autolimits!` on all axes to clear one-way limits.

## [5.0.0] - 2026-03-02

### Breaking Changes
- **Animation extension renamed**: `WeberElectrodynamicsGLMakieExt` → `WeberElectrodynamicsMakieExt`.
  The weak dependency is now `Makie` instead of `GLMakie`, so any Makie backend (GLMakie,
  CairoMakie, WGLMakie) triggers the extension. Users who `import`ed the old extension module
  by name must update.

### Changed
- Animation extension depends on abstract `Makie` API instead of `GLMakie` backend directly.
  Users choose their backend in the notebook (`using GLMakie`, `using CairoMakie`, etc.)
  and the extension activates automatically.
- Compat bounds: `Makie = "0.21, 0.22, 0.23, 0.24"` (replaces `GLMakie = "0.9, 0.10"`).
- Tested on Julia 1.12.5 with Symbolics v7.15.3, Makie v0.24.8, GLMakie v0.13.8.

## [4.2.0] - 2026-03-01

### Added
- GLMakie animation extension with streaming and replay modes for interactive
  real-time visualization of Weber electrodynamics simulations.
- `animate_weber(prob)` for live streaming integration with rolling buffer display.
- `animate_weber(sol)` for replaying pre-computed solutions.
- Multi-panel dashboard: particle trajectories, energy & momentum, angular momentum,
  and selectable phase space (pair or per-particle).
- Play/pause/reset controls with live error tracking (energy, momentum, angular momentum).
- Three-body bound state exploration notebooks (`examples/three_body_bound_states/`).

## [4.1.0] - 2026-03-01

### Added
- Collision bounce for sub-critical like-charge oscillation: pre-step reflection of
  relative coordinates through the origin when pair separation drops below a configurable
  radius, analytically continuing C⁰ head-on collisions (ℓ=0) without integrating
  through the r=0 singularity.
- `collision_bounce_radius` field on `RegularizationOptions`.
- `regularization_collision_bounce_radius` keyword argument on `WeberProblem`
  (default `0.0` = disabled).
- General-dimension helpers `_current_pair_r`, `_reflect_pair!`, and
  `_apply_collision_bounces!` in the solver.
- Test for sub-critical like-charge oscillation with collision bounce.
- Critical radius reference notebook (`examples/critical_radius_reference.ipynb`)
  demonstrating molecular oscillation and super-critical scattering regimes.
- Exploratory doc `docs/exploratory/CollisionBounceRegularization.md` with lessons
  learned on regularization approaches for Weber's velocity-dependent force.

### Fixed
- 2D trajectory plot now computes tight axis limits from data instead of letting
  `aspect_ratio=:equal` produce excessive whitespace for small-scale trajectories.

### Removed
- Stale exploratory docs (`ParallelizationResearchReport.md`,
  `weber_dynamics_rust_design.md`).

## [4.0.0] - 2026-02-19

### Breaking Changes
- **`params` vector structure extended**: now `[m₁,…,mₙ, q₁,…,qₙ, c, κ₁₂, κ₁₃, …, κ_{N-1,N}]`
  (length `2N + 1 + N(N-1)/2`). Direct callers of `system.hamiltonian_compiled` or
  `system.dp_dt_compiled` must append per-pair κ values (all `1.0` to reproduce standard Weber).
- `compute_pair_weber_components` return type changed from 3-tuple to 4-tuple: added
  `zollner_extra_potential` as the fourth element.
- `PairEnergyData` struct gained new fields `kappa`, `zollner_extra_potential`, and
  `zollner_extra_force`; code that pattern-matches on the struct must be updated.

### Added
- `ZollnerOptions` struct: `ZollnerOptions(enabled=false, a=0.0)` for configuring the
  Zöllner electrogravitational mismatch parameter.
- `zollner_enabled` and `zollner_a` keyword arguments on `WeberProblem`.
- `kappas::Vector{Float64}` field on `WeberProblem` storing per-pair coupling factors κ_ij.
- `_pair_index(i, j, n)` helper for pair indexing in the κ vector.
- Four new plot functions: `plot_zollner_energy`, `plot_zollner_force_residual`,
  `plot_weber_vs_zollner`, `plot_zollner_phase_space`.
- Comprehensive Zöllner test suite in `test/test_zollner.jl`.

### Fixed
- `plot_weber_vs_zollner`: trajectory arrays were indexed as `[dim, :]` (2 points) instead of
  `[:, dim]` (all timesteps), producing straight line segments instead of full orbits.

### Changed
- Rewrote `examples/four_body_regularized_reference.ipynb`: analytically-derived initial
  conditions (closed-form speed from target energy ratio η = KE/|PE₀|), square geometry with
  alternating (+,−,+,−) charges, four Zöllner runs (a = 0, 0.02, 0.05, 0.10), extensive
  relative-coordinate analysis (pairwise separations, phase portraits, Zöllner scaling).

## [3.0.1] - 2026-02-18

### Fixed
- Prevent regularization dispatch from crashing on single-particle (`n_pairs=0`) systems
- Recompute encounter mode each step so active components can transition between pair and chain modes
- Enforce documented Cartesian fallback when overlapping components are active and `regularization_chain_enabled=false`
- Apply `regularization_constraint_tolerance` to KS residual diagnostics thresholding

## [3.0.0] - 2026-02-12

### Breaking Changes
- Add regularization configuration to `WeberProblem` and enable regularization by default
- Extend `WeberSolution` with regularization diagnostics payload
- Integrator stepping now dispatches between Cartesian and regularized encounter paths

### Added
- `RegularizationOptions` public type with encounter/switching/substep controls
- `RegularizationDiagnostics` public type with activation, mode, substep, and KS-constraint metrics
- Adaptive close-encounter dispatch with hysteresis (`r_on`/`r_off`)
- Pair regularization path (1D square-map scaffold, 2D Levi-Civita lift, 3D KS lift with constraint projection)
- Chain encounter fallback for overlapping close-pair components
- Preallocated regularization buffers integrated into solver workspace
- New theory/design documentation for regularized integration workflow
- Regularization backend selection: `:lifted_pair` and `:adaptive_cartesian`
- Backend fallback policy for unsupported lifted dimensions (1D/3D) with optional warning
- 2D true lifted pair split stepping path (external midpoint + LC lifted midpoint)
- Extended diagnostics for backend usage and fallback accounting
- New 4-body singular encounter reference notebook: `examples/four_body_regularized_reference.ipynb`

### Changed
- Improved 2D Levi-Civita lift branch handling near singular/axis-degenerate states
- Stabilized lifted substep time scaling by freezing the monitor within each substep

## [2.1.0] - 2026-02-03

### Added
- `MomentumData` struct for total linear and angular momentum analysis
- `compute_momentum_timeseries` function for computing momentum conservation
- `plot_momentum` with 2-panel layout (linear momentum components + angular momentum)
- Support for 1D (linear only), 2D (Lz scalar), and 3D (L vector) momentum analysis
- Momentum tests including conservation checks for isolated systems
- Momentum analysis section in example notebook

### Changed
- Remove docstrings from `energy.jl` and `momentum.jl` (documentation in code comments)
- Remove problematic unicode comment from `energy.jl`

## [2.0.0] - 2026-02-01

### Breaking Changes
- Remove `ForceData`, `compute_force_timeseries`, `NewtonsThirdLawData`, `check_newtons_third_law`
- Remove `compute_phase_space_data` standalone function
- Remove `plot_forces` (replaced by `plot_pair_forces`)
- Change `plot_phase_space` to accept `PairForceData` instead of `PhaseSpaceData`

### Added
- `PairForceData` struct for single-pair force analysis with Weber force decomposition
- `ForceStatistics` struct with min, max, mean, range for force magnitude
- `compute_pair_force_timeseries` function for computing forces on a single particle pair
- Weber force vector form decomposition (Coulomb + v·v + r·a + rv² terms)
- Weber force radial form decomposition (Coulomb + ṙ² + r·r̈ terms)
- `PhaseSpaceData` embedded in `PairForceData` to avoid redundant computation
- `plot_pair_forces` with 4-panel vertical layout (magnitude, components, vector form, radial form)

### Changed
- Phase space data (r, ṙ, θ, L) now computed once within force computation loop
- Signed scalar plotting for force decomposition terms (positive = repulsion, negative = attraction)

## [1.0.0] - 2026-02-01

### Breaking Changes
- Make `WeberSystem` symbolic and move physics parameters to `WeberProblem`
- Simplify API: remove extensible solver system and `@hamiltonian` macro
- Remove nonlinear solve extensibility infrastructure
- Hardcode Weber Hamiltonian for performance and simplicity
- Remove internal exports; all API access through notebooks
- Refactor naming conventions for clarity and physics alignment
- Align internal naming with paper notation (Jayawardana-Ohsawa 2023)
- Remove unused `t` parameter from energy functions

### Added
- Compiled Hamiltonian with rewritten energy analysis and optimized statistics
- Diff buffer for improved memory management
- Top-level Makefile with `cwe-` prefixed commands
- JuliaFormatter integration with format command
- VSCode settings for LaTeX builds
- Figure preview workflow for paper development
- Exploratory Rust crate design for Weber dynamics
- Parallelization research report

### Changed
- Improve plot presentation for scientific publication quality
- Optimize integrator by hoisting closures and adding `@inbounds`
- Reorganize paper directory structure for reproducibility

### Fixed
- Fix plot styling issues
- Fix remaining memory allocation issues
- Fix LaTeX formatting in Weber electrodynamics documentation

### Removed
- Remove LaTeXStrings dependency
- Remove docstrings and comments from src/ Julia files (moved to documentation)

### Documentation
- Extend semi-explicit integrator reference documentation
- Clarify sign convention in Weber force derivation
- Extend theory documentation to 3D and add missing equations
- Improve Weber Force section introduction

## [0.1.0] - 2026-01-20

### Added
- Initial release
- `SymmetricProjection` semi-explicit symplectic integrator for non-separable Hamiltonians
- `@hamiltonian` macro and `build_hamiltonian` function for Hamiltonian construction
- `WeberProblem`, `WeberSolution`, `WeberIntegrator` types
- CommonSolve.jl interface (`solve`, `init`, `step!`, `solve!`)
- Statistics module:
  - `TrajectoryData` and `create_trajectory_data`
  - `EnergyData` and `compute_energy_timeseries`
  - `ForceData` and `compute_force_timeseries`
  - `NewtonsThirdLawData` and `check_newtons_third_law`
  - `PhaseSpaceData` and `compute_phase_space_data`
- Plots.jl extension with `plot_trajectories`, `plot_energy`, `plot_forces`, `plot_phase_space`
- `RelaxedFixedPoint` nonlinear solver with configurable relaxation parameter
- Support for SimpleNonlinearSolve.jl algorithms
