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
export compute_energy_timeseries

include("statistics/forces.jl")
export ForceData
export compute_force_timeseries
export NewtonsThirdLawData
export check_newtons_third_law

include("statistics/phase_space.jl")
export PhaseSpaceData
export compute_phase_space_data

# =============================================================================
# Plotting (Extension)
# =============================================================================
function plot_trajectories end
function plot_energy end
function plot_forces end
function plot_phase_space end
export plot_trajectories, plot_energy, plot_forces, plot_phase_space

end
