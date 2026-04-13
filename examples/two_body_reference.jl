#!/usr/bin/env julia
# Canonical two-body Weber reference problem.
#
# See two_body_reference.md for the specification and
# two_body_reference.ipynb for an annotated tutorial.
#
# Usage:
#   julia --project=. examples/two_body_reference.jl            # save PNGs only
#   WEBER_DISPLAY=1 julia --project=. examples/two_body_reference.jl  # also open plot windows

using WeberElectrodynamics
using Plots
using Printf

const m1, m2 = 1.0, 1.0
const q1, q2 = 1.0, -1.0
const c = 4.0

const q0 = [-1.0, 0.0, 1.0, 0.0]
const p0 = [ 0.0, -0.25, 0.0, 0.25]

const mu = m1 * m2 / (m1 + m2)
const a_semi = 8 / 7
const T_orbit = 2π * sqrt(a_semi^3 * mu / abs(q1 * q2))

system = WeberSystem(2, 2)
prob = WeberProblem(
    system, (0.0, 5 * T_orbit), q0, p0;
    masses = [m1, m2], charges = [q1, q2], c = c, dt = 0.001,
)

sol = solve(prob)
@assert sol.retcode == :Success "solve failed: $(sol.retcode)"

traj   = compute_trajectory_data(sol, 2, 2; stride = 10)
energy = compute_energy_timeseries(sol; stride = 10)
forces = compute_pair_force_timeseries(sol, (1, 2), 2, 2, [m1, m2], [q1, q2], c; stride = 10)
mom    = compute_momentum_timeseries(sol; stride = 10)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

figures = [
    ("trajectories.png",  plot_trajectories(traj)),
    ("energy.png",        plot_energy(energy)),
    ("pair_energy.png",   plot_pair_energy(energy, (1, 2))),
    ("energy_errors.png", plot_energy_errors(energy)),
    ("pair_forces.png",   plot_pair_forces(forces)),
    ("phase_space.png",   plot_phase_space(forces)),
    ("momentum_errors.png", plot_momentum_errors(mom)),
]

const display_plots = get(ENV, "WEBER_DISPLAY", "") == "1"
for (name, fig) in figures
    path = joinpath(outdir, name)
    savefig(fig, path)
    println("wrote ", path)
    display_plots && display(fig)
end

@printf("\nenergy drift (max):  %.2e %%\n", energy.statistics.global_error_percent_max)
# Momentum drift is reported in the plot_momentum_errors legend; see momentum_errors.png.
