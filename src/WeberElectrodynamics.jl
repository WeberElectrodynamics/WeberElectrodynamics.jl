module WeberElectrodynamics

# =============================================================================
# Weber System (n-body Hamiltonian)
# =============================================================================
include("weber_system.jl")
export WeberSystem

# =============================================================================
# Core Types (Problem, Solution, Integrator)
# =============================================================================
include("types.jl")
export WeberProblem
export WeberSolution
export WeberIntegrator
export RegularizationOptions
export RegularizationDiagnostics
export ZollnerOptions

# =============================================================================
# Regularization Helpers
# =============================================================================
include("regularization.jl")

# =============================================================================
# CommonSolve Interface
# =============================================================================
include("solve.jl")
# Re-export CommonSolve functions for convenience
using CommonSolve: solve, init, step!, solve!
export solve, init, step!, solve!

# =============================================================================
# Statistics & Analysis
# =============================================================================
include("statistics/trajectories.jl")
export TrajectoryData
export compute_trajectory_data

include("statistics/energy.jl")
export EnergyData
export PairEnergyData
export EnergyStatistics
export compute_energy_timeseries

include("statistics/forces.jl")
export PairForceData
export ForceStatistics
export PhaseSpaceData
export compute_pair_force_timeseries

include("statistics/momentum.jl")
export MomentumData
export compute_momentum_timeseries

# =============================================================================
# Plotting (Extension)
# =============================================================================
function plot_trajectories end
function plot_energy end
function plot_pair_energy end
function plot_energy_errors end
function plot_pair_forces end
function plot_phase_space end
function plot_momentum end
function plot_zollner_energy end
function plot_zollner_force_residual end
function plot_weber_vs_zollner end
function plot_zollner_phase_space end
export plot_trajectories, plot_energy, plot_pair_energy, plot_energy_errors, plot_pair_forces, plot_phase_space, plot_momentum
export plot_zollner_energy, plot_zollner_force_residual, plot_weber_vs_zollner, plot_zollner_phase_space

end
