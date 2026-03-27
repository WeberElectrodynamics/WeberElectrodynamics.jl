# Quick Start

## Two-body simulation

```julia
using WeberElectrodynamics

# Build the system (symbolic differentiation happens here — takes a few seconds)
sys = WeberSystem(2, 2)   # 2 particles, 2D

# Initial conditions: particles at (±1, 0) with transverse momenta
q0 = [1.0, 0.0, -1.0, 0.0]
p0 = [0.0, 0.5,  0.0, -0.5]

prob = WeberProblem(
    sys, (0.0, 10.0), q0, p0;
    masses  = [1.0, 1.0],
    charges = [1.0, -1.0],
    c       = 10.0,
    dt      = 0.01,
)

sol = solve(prob, SymmetricProjectionIntegrator())
println(sol)
# WeberSolution with 1001 timesteps (retcode: Success)
```

## Analysing the result

```julia
# Energy conservation
energy = compute_energy_timeseries(sol)
println("Max energy drift: ", energy.statistics.global_error_percent_max, " %")

# Trajectory data
traj = compute_trajectory_data(sol, 2, 2)

# Plotting (requires Plots.jl)
using Plots
plot_trajectories(traj)
plot_energy(energy)
```

## Step-by-step integration

```julia
integrator = init(prob, SymmetricProjectionIntegrator())
while step!(integrator)
    # inspect integrator.q, integrator.p, integrator.t at each step
end
sol = solve!(integrator)
```

## Static plots (requires Plots.jl)

```julia
using Plots

traj     = compute_trajectory_data(sol, 2, 2)
energy   = compute_energy_timeseries(sol)
momentum = compute_momentum_timeseries(sol)

plot_trajectories(traj)     # particle paths in configuration space
plot_energy(energy)         # Hamiltonian drift over time
plot_momentum(momentum)     # total linear momentum components
```

## Interactive animation (requires GLMakie)

```julia
using GLMakie
animate_weber(prob)   # live streaming
animate_weber(sol)    # replay
```

## Optional features

The integrator above runs the core unregularized symplectic method. Two optional
extensions can be enabled explicitly:

- **[Regularization](regularization.md)** — Levi-Civita / KS handling for
  close encounters. Pass `regularization_enabled = true` to `WeberProblem`.
- **[Zöllner Extension](zollner.md)** — research feature implementing Zöllner's
  electrogravitational mismatch hypothesis. Pass `zollner_enabled = true` and
  `zollner_a = <value>` to `WeberProblem`.
