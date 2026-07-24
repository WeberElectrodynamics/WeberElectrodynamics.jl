# Statistics

For generic custom Hamiltonians, `compute_energy_timeseries` evaluates the
compiled total Hamiltonian and kinetic split. Full pair decompositions are
available when the system includes built-in Weber/Zöllner `NamedTerm`
decomposition closures.

## Trajectories

```@docs
TrajectoryData
compute_trajectory_data
```

## Energy

```@docs
EnergyData
PairEnergyData
EnergyStatistics
compute_energy_timeseries
conservation_summary
```

## Forces

```@docs
PairForceData
ForceStatistics
PhaseSpaceData
compute_pair_force_timeseries
```

## Momentum

```@docs
MomentumData
compute_momentum_timeseries
```
