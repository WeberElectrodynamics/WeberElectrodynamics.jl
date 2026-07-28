# Statistics

All reported velocities, radial velocities, and kinetic energies are
**physical** quantities recovered from the canonical momenta — see
[The Weber Hamiltonian](@ref). In particular:

- `EnergyData.kinetic_energy` is `Σ ½ mᵢ|vᵢ|²`, not `Σ |pᵢ|²/(2mᵢ)`.
- `PairEnergyData.radial_velocity` is the physical `ṙ`, not
  `r̂·(pᵢ/mᵢ − pⱼ/mⱼ)`.
- `PairForceData` velocities and accelerations come from `physical_velocities`.

Because the manual decomposition `Σ ½ mᵢ|vᵢ|² + Σ Uᵢⱼ` equals the compiled
canonical `H(q, p)` by the Legendre transform,
`EnergyData.hamiltonian_validation_error` is a genuine independent cross-check
of the compiled Hamiltonian rather than a restatement of it.

For generic custom Hamiltonians, `compute_energy_timeseries` evaluates the
compiled total Hamiltonian and falls back to the canonical kinetic split. Full
pair decompositions are available when the system includes the built-in Weber
`NamedTerm` closures.

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
