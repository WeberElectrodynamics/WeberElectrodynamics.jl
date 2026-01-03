module WeberElectrodynamics

include("hamiltonian_systems.jl")
export SymbolicHamiltonian, HamiltonianVectorField
export symbolize, compile, create_phase_space_variables

include("symplectic_integrator.jl")
export PhaseSpacePoint, IntegratorSettings, IntegratorState, TimeSpan, IntegratorSolution
export step!, integrate

include("statistics/trajectory_data.jl")
export TrajectoryData
export create_trajectory_data

include("statistics/energy_data.jl")
export EnergyData
export compute_energy_timeseries

include("statistics/force_data.jl")
export ForceData
export compute_force_timeseries
export NewtonsThirdLawData
export check_newtons_third_law

include("visualization.jl")
export plot_trajectories_2d, plot_energy, plot_newtons_third_law

end
