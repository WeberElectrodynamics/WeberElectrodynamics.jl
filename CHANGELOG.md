# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
