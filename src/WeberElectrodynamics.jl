module WeberElectrodynamics

# =============================================================================
# Nonlinear Solver (internal)
# =============================================================================
include("nonlinear_solvers.jl")

# =============================================================================
# Core Types
# =============================================================================
include("types.jl")
export WeberAlgorithm, SymmetricProjectionIntegrator
export WeberHamiltonian
export WeberProblem
export WeberSolution
export WeberIntegrator
export NonlinearSolveError

# =============================================================================
# Hamiltonian Construction
# =============================================================================
include("hamiltonian_systems.jl")
export compile_hamiltonian
export generate_phase_space_symbols

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
