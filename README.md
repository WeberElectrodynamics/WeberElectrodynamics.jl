# WeberElectrodynamics.jl

A Julia package for simulating charged particle dynamics using Weber's velocity-dependent electrodynamics.

## Installation

```julia
using Pkg
Pkg.add("WeberElectrodynamics")
```

## Quick Start

```julia
using WeberElectrodynamics

# Define a two-body Hamiltonian
H = build_hamiltonian(2, 2; param_names=[:m1, :m2, :k, :c]) do q, p, params
    m1, m2, k, c = params
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    # Kinetic energy
    KE = (px1^2 + py1^2)/(2m1) + (px2^2 + py2^2)/(2m2)

    # Weber potential
    dx, dy = x1 - x2, y1 - y2
    r = sqrt(dx^2 + dy^2)
    vx1, vy1 = px1/m1, py1/m1
    vx2, vy2 = px2/m2, py2/m2
    rdot = (dx*(vx1-vx2) + dy*(vy1-vy2)) / r
    PE = k / r * (1 - rdot^2 / (2*c^2))

    KE + PE
end

# Set up initial conditions (circular orbit)
m1, m2, k, c = 1.0, 0.1, -0.1, 4.0
r0 = 2.0
M = m1 + m2
v = sqrt(abs(k) * M / (m1 * m2 * r0))

q₀ = [-m2/M * r0, 0.0, m1/M * r0, 0.0]
p₀ = [0.0, -m1 * m2/M * v, 0.0, m2 * m1/M * v]

# Create and solve problem
prob = WeberProblem(H, (0.0, 10.0), q₀, p₀; params=[m1, m2, k, c], dt=0.001)
sol = solve(prob)

# Access solution
println("Solved $(length(sol)) timesteps")
for (t, q, p) in sol
    # Process each timestep...
end
```

## Features

- **SymmetricProjection integrator**: Semi-explicit symplectic integrator for non-separable Hamiltonians
- **CommonSolve.jl interface**: Standard `solve`, `init`, `step!`, `solve!` functions
- **Analysis tools**: Energy timeseries, force computation, Newton's third law verification, phase space data
- **Plots.jl extension**: `plot_trajectories`, `plot_energy`, `plot_forces`, `plot_phase_space`

## API

### Hamiltonian Construction

```julia
# Function API
H = build_hamiltonian(H_func, n_particles, dims; param_names=[:m, :k])

# Macro API
H = @hamiltonian 2 2 [:m1, :m2, :k] (q, p, params) -> begin
    # Hamiltonian expression
end
```

### Problem & Solution

```julia
prob = WeberProblem(H, tspan, q₀, p₀; params, dt, tolerance=1e-13, max_iterations=100)
sol = solve(prob)  # Returns WeberSolution

# Access solution data
sol.t        # Time points
sol.q        # Position vectors
sol.p        # Momentum vectors
sol.retcode  # :Success or :Failure
```

### Stepped Integration

```julia
integrator = init(prob)           # Initialize
while step!(integrator)           # Advance one step
    # Access integrator.t, integrator.q, integrator.p
end
sol = solve!(integrator)          # Finalize solution
```

### Statistics & Analysis

```julia
# Trajectory data
traj = create_trajectory_data(sol, n_particles, dims; stride=1)

# Energy analysis
energy = compute_energy_timeseries(sol, E_func, KE_func, PE_func, params)

# Force analysis
forces = compute_force_timeseries(sol, n_particles, dims, masses, charges, c)
n3_law = check_newtons_third_law(forces)

# Phase space
ps = compute_phase_space_data(sol, n_particles, dims, masses)
```

### Plotting (requires Plots.jl)

```julia
using Plots
plot_trajectories(traj)
plot_energy(energy)
plot_forces(forces)
plot_phase_space(ps)
```

## License

MIT
