module WeberElectrodynamics

# =============================================================================
# Nonlinear Solvers
# =============================================================================
include("nonlinear_solvers.jl")
export RelaxedFixedPoint
export solve_nonlinear!

# =============================================================================
# Core Types
# =============================================================================
include("types.jl")
export WeberAlgorithm, SymmetricProjection
export WeberHamiltonian
export WeberProblem
export WeberSolution
export WeberIntegrator
export IntegratorBuffers  # For advanced users
export NonlinearSolveError
export NewtonConvergenceError  # Deprecated alias

# =============================================================================
# Hamiltonian Construction
# =============================================================================
include("hamiltonian_systems.jl")
export @hamiltonian, build_hamiltonian
export create_phase_space_variables

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
export create_trajectory_data

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
