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

## Setting up initial conditions

The Weber Hamiltonian depends on the radial velocity ṙᵢⱼ between every pair, which complicates
choosing q0/p0 for a target energy. The standard approach avoids this:

**Zero-radial-velocity principle.** If ṙᵢⱼ(0) = 0 for every pair, the Weber potential reduces
to the Coulomb potential at t = 0. This is guaranteed by setting momenta tangential to the
separation vectors (2D), or by using rigid-rotation momenta (3D). All constructions below
satisfy this principle.

**Centre-of-mass frame.** Work with the COM at the origin:

```
∑ mᵢ rᵢ(0) = 0      (position)
∑ pᵢ(0)    = 0      (momentum)
```

Tangential and rigid-rotation constructions satisfy both conditions automatically when
particles are placed symmetrically around the origin.

### Two-body orbits

For masses `m1`, `m2`, charges `q1`, `q2`, and initial separation `r0`:

| Quantity | Formula |
|----------|---------|
| Reduced mass | `μ = m1*m2/(m1+m2)` |
| Pair coupling | `k = κ₁₂ * q1 * q2`  (κ=1 standard Weber; κ=1+a for Zöllner unlike pairs) |
| Circular momentum | `p_circ = √(μ * \|k\| / r0)` |

Scale `p_circ` by a dimensionless factor η_v to select the orbit type:

| η_v | Orbit type |
|-----|-----------|
| 0 | radial free fall |
| (0, 1) | elliptical, sub-circular (bound) |
| 1 | circular |
| (1, √2) | elliptical, super-circular (bound) |
| √2 | parabolic escape |
| > √2 | hyperbolic (unbound) |

From a target energy `E < 0`: `η_v = √(2 * (1 + E * r0 / |k|))`.

```julia
# Two-body unlike-charge orbit in 2D
m1, m2 = 1.0, 1.0
q1, q2 = 1.0, -1.0
r0     = 2.0

μ       = m1 * m2 / (m1 + m2)
k       = q1 * q2                        # < 0 for unlike charges
p_circ  = sqrt(μ * abs(k) / r0)          # circular orbit momentum magnitude

# Place along x-axis, COM at origin
q0 = [-m2/(m1+m2)*r0, 0.0,  m1/(m1+m2)*r0, 0.0]
p0 = [0.0, p_circ, 0.0, -p_circ]

# Scale by η_v to select orbit type
η_v = 1.0           # circular; use 0.8 for sub-circular ellipse, √2 for escape
p0  = p0 .* η_v
```

### N-body symmetric polygon

Place N particles at regular N-gon vertices of circumradius `R` with alternating charges ±Q.
Assign tangential momenta with common speed `v`:

```julia
using LinearAlgebra

N, R, Q, m = 4, 2.0, 1.0, 1.0
angles = [2π*(i-1)/N for i in 1:N]
q0 = vcat([[R*cos(θ), R*sin(θ)] for θ in angles]...)
charges = [(-1)^(i+1) * Q for i in 1:N]

# Compute initial potential energy U0 (pure Coulomb at t=0)
U0 = sum(
    charges[i]*charges[j] / norm(q0[2i-1:2i] - q0[2j-1:2j])
    for i in 1:N for j in i+1:N
)

η  = 0.5                               # T0/|U0|; η < 1 → bound orbit
v  = sqrt(2η * abs(U0) / (N * m))
p0 = vcat([m*v .* [-sin(θ), cos(θ)] for θ in angles]...)
```

See [theory/InitialConditions.md](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/InitialConditions.md)
for full derivations including polygon U₀ formulas, 3D rigid-rotation construction, and the
energy–angular-momentum parameterisation.

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
plot_energy(energy)         # T, U, H components over time
plot_energy_errors(energy)  # local + relative energy error (log scale)
plot_momentum_errors(momentum)  # drift of linear & angular momentum from initial
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
  close encounters. Pass `regularization = RegularizationOptions(enabled = true)`
  to `WeberProblem`.
- **[Zöllner Extension](zollner.md)** — research feature implementing Zöllner's
  electrogravitational mismatch hypothesis. Pass
  `zollner = ZollnerOptions(enabled = true, a = <value>)` to `WeberProblem`.
