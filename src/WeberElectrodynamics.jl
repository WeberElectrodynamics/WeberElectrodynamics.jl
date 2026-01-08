module WeberElectrodynamics

include("hamiltonian_systems.jl")
export SymbolicHamiltonian, HamiltonianVectorField
export symbolize, compile, create_phase_space_variables

include("symplectic_integrator.jl")
export PhaseSpacePoint, IntegratorSettings, IntegratorState, TimeSpan, IntegratorSolution
export NewtonConvergenceError
export step!, integrate

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

# Plots.jl plotting for Jupyter notebooks
include("plotting.jl")
export plot_trajectories, plot_energy, plot_forces, plot_phase_space

end
