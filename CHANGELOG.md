# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
