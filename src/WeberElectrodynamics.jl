module WeberElectrodynamics

# Stdlib re-exports (needed by package extension)
using Printf: @sprintf
using LinearAlgebra: norm, dot
export @sprintf, norm, dot

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
export RegularizationOptions
export RegularizationDiagnostics
export ZollnerOptions
export SymmetricProjectionIntegrator

# =============================================================================
# Regularization Helpers
# =============================================================================
include("regularization.jl")

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
export PairEnergyData
export EnergyStatistics
export compute_energy_timeseries

include("statistics/forces.jl")
export PairForceData
export ForceStatistics
export PhaseSpaceData
export compute_pair_force_timeseries

include("statistics/momentum.jl")
export MomentumData
export compute_momentum_timeseries

# =============================================================================
# Plotting (Extension)
# =============================================================================

"""
    plot_trajectories(data::TrajectoryData) -> Plot

Plot particle trajectories in 1D, 2D, or 3D based on `data.dims`.

- **1D**: Position vs. time for each particle.
- **2D**: x–y scatter with equal aspect ratio and start/end markers.
- **3D**: x–y–z line plot with start/end markers.

Requires `using Plots` to activate the Plots.jl extension.

# Arguments
- `data::TrajectoryData`: Trajectory data from `compute_trajectory_data`.

# Returns
- A `Plots.Plot` object.
"""
function plot_trajectories end

"""
    plot_energy(data::EnergyData) -> Plot

Plot total n-body energy timeseries with comprehensive error visualization.

Two-panel layout: energy components (Total, Kinetic, Potential) and relative
error on a log scale with max/avg reference lines.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_energy end

"""
    plot_pair_energy(data::EnergyData, pair::Tuple{Int,Int}) -> Plot

Plot the Weber potential decomposition for a single particle pair.

Two-panel layout: Coulomb term, velocity correction, and total pair potential vs.
time; and radial velocity ṙ with approach/separation shading.

Requires `using Plots` to activate the Plots.jl extension.

# Arguments
- `data::EnergyData`: Energy data from `compute_energy_timeseries`.
- `pair::Tuple{Int,Int}`: Particle pair indices, e.g., `(1, 2)`.
"""
function plot_pair_energy end

"""
    plot_energy_errors(data::EnergyData) -> Plot

Plot comprehensive energy error analysis (local error, global drift, and
Hamiltonian validation) in a three-panel log-scale layout.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_energy_errors end

"""
    plot_pair_forces(data::PairForceData) -> Plot

Plot the Weber force decomposition for a particle pair in a four-panel layout:
force magnitude, force components, vector form decomposition, and radial form
decomposition.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_pair_forces end

"""
    plot_phase_space(data::PairForceData) -> Plot

Plot the (r, ṙ) phase-space portrait for a particle pair with marked initial
and final states.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_phase_space end

"""
    plot_momentum_errors(data::MomentumData) -> Plot

Plot conservation errors for linear and angular momentum as two stacked
log-scale panels (top: linear drift `‖ΔP‖`, bottom: angular drift `|ΔLz|`
in 2D or `‖ΔL‖` in 3D). Each legend reports max absolute drift and — when
the initial magnitude is nonzero — max relative drift `max_t ‖Δ·‖ / ‖·₀‖`.
1D systems show only the linear panel.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_momentum_errors end

"""
    plot_zollner_energy(data::EnergyData) -> Plot

Visualise Zöllner potential contributions against standard Weber energies.

Two-panel layout: total potential with per-pair breakdown and Zöllner extra terms
(κ values shown); Zöllner gravitational residual ΣΔV_Z compared to total energy.
Only meaningful when `prob.zollner.enabled == true`.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_zollner_energy end

"""
    plot_zollner_force_residual(data::PairForceData) -> Plot

Two-panel visualization of the Zöllner extra force for one particle pair:
total force magnitude alongside the Zöllner extra magnitude |(κ−1)·F_Coulomb|,
and the ratio of extra Zöllner force to total force.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_zollner_force_residual end

"""
    plot_weber_vs_zollner(sol1, sol2; labels=["Weber","Zöllner"]) -> Plot

Overlay trajectories from two solutions (same initial conditions, different κ)
to visualise the orbital divergence introduced by the Zöllner mismatch.
Solid lines = `sol1`, dashed lines = `sol2`. Requires 2D or 3D.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_weber_vs_zollner end

"""
    plot_zollner_phase_space(data1, data2; labels=["Weber","Zöllner"]) -> Plot

Overlay (r, ṙ) phase portraits for the same pair from two simulations, showing
how the Zöllner mismatch shifts the phase-space orbit.

Requires `using Plots` to activate the Plots.jl extension.
"""
function plot_zollner_phase_space end

export plot_trajectories, plot_energy, plot_pair_energy, plot_energy_errors, plot_pair_forces, plot_phase_space, plot_momentum_errors
export plot_zollner_energy, plot_zollner_force_residual, plot_weber_vs_zollner, plot_zollner_phase_space

# =============================================================================
# Animation (Makie Extension)
# =============================================================================

"""
    animate_weber(prob_or_sol; kwargs...) -> screen

Launch an interactive animated dashboard for a Weber simulation.

Two modes:
- `animate_weber(prob::WeberProblem; ...)` — **streaming**: runs the integrator
  live and displays data in real time. Requires a windowed backend (GLMakie or WGLMakie).
- `animate_weber(sol::WeberSolution; ...)` — **replay**: replays a completed solution.

The dashboard shows particle trajectories, energy, momentum, and phase-space panels.
Requires at least 2D (`dims ≥ 2`).

**Keywords (both modes)**:
- `buffer_size=2000`: Rolling display window size.
- `tail_length=200`: Visible trajectory tail length (≤ `buffer_size`).
- `compute_batch=1`: Steps computed or advanced per animation frame.
- `initial_pair=(1,2)`: Default pair for the phase-space panel.
- `phase_mode=:pair`: `:pair` or `:particle` phase-space mode.
- `initial_particle=1`, `initial_component=1`: Phase-space selection.
- `figure_size=(1400,900)`: Window size in pixels.

**Streaming-only keywords**:
- `alg=SymmetricProjectionIntegrator()`: Integrator algorithm.

**Replay-only keywords**:
- `stride=1`: Skip every `stride`-th stored timestep during replay.

Requires `using GLMakie` (or another Makie backend) to activate the extension.

# Returns
- A Makie screen handle; close the window to stop.
"""
function animate_weber end
export animate_weber

end
