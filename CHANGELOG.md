# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
