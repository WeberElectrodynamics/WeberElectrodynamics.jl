#!/usr/bin/env julia
#
# Reproduces examples/two_body_reference.ipynb and exports the numerical
# data backing the paper's two-body reference figures as tab-separated
# .dat files consumed by pgfplots. Run from the repository root:
#
#   julia --project=. papers/Computational-Weber-Electrodynamics/scripts/export_two_body_data.jl
#
# Outputs (relative to this script):
#   ../data/two_body_trajectories.dat
#   ../data/two_body_energy.dat
#   ../data/two_body_energy_errors.dat

using WeberElectrodynamics
using Printf

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
isdir(DATA_DIR) || mkpath(DATA_DIR)

# --- Canonical two-body parameters (see two_body_reference.ipynb) ---
const m1, m2 = 1.0, 1.0
const q1, q2 = 1.0, -1.0
const c = 4.0

const q0 = [-1.0, 0.0, 1.0, 0.0]
const p0 = [0.0, -0.25, 0.0, 0.25]

const mu = m1 * m2 / (m1 + m2)
const a_semi = 8 / 7
const T_orbit = 2π * sqrt(a_semi^3 * mu / abs(q1 * q2))
const tspan = (0.0, 5 * T_orbit)

const STRIDE = 10

const PARAM_HEADER = """
# Two-body Weber reference run
# m1=$(m1)  m2=$(m2)  q1=$(q1)  q2=$(q2)  c=$(c)
# q0=$(q0)
# p0=$(p0)
# tspan=$(tspan)  (5 Kepler periods, T=$(round(T_orbit, digits=6)))
# dt=1e-3  stride=$(STRIDE)
"""

function write_dat(path::AbstractString, header::AbstractString,
                   colnames::Vector{String}, cols::Vector{<:AbstractVector})
    n = length(first(cols))
    @assert all(length(c) == n for c in cols) "column length mismatch"
    open(path, "w") do io
        print(io, PARAM_HEADER)
        println(io, "# ", header)
        println(io, join(colnames, "\t"))
        for i in 1:n
            row = join((@sprintf("%.10e", cols[k][i]) for k in eachindex(cols)), "\t")
            println(io, row)
        end
    end
    @info "wrote" path n_rows=n
end

function main()
    system = WeberSystem(2, 2)
    prob = WeberProblem(system, tspan, q0, p0;
        masses=[m1, m2], charges=[q1, q2], c=c, dt=0.001)

    sol = solve(prob)
    @assert sol.retcode == :Success "solver did not succeed: $(sol.retcode)"
    @info "solved" retcode=sol.retcode n_steps=length(sol.t)

    # --- Trajectories ---
    traj = compute_trajectory_data(sol, 2, 2; stride=STRIDE)
    t_sampled = sol.t[1:STRIDE:end]
    x1 = traj.trajectories[1][:, 1]
    y1 = traj.trajectories[1][:, 2]
    x2 = traj.trajectories[2][:, 1]
    y2 = traj.trajectories[2][:, 2]
    write_dat(joinpath(DATA_DIR, "two_body_trajectories.dat"),
        "particle positions sampled every $(STRIDE) steps",
        ["t", "x1", "y1", "x2", "y2"],
        [t_sampled, x1, y1, x2, y2])

    # --- Energy ---
    energy = compute_energy_timeseries(sol; stride=STRIDE)
    H = energy.total_energy
    T = energy.kinetic_energy
    U = energy.total_potential_energy
    write_dat(joinpath(DATA_DIR, "two_body_energy.dat"),
        "total Hamiltonian H = T + U with kinetic T and potential U components",
        ["t", "H", "T", "U"],
        [energy.t, H, T, U])

    # --- Energy errors ---
    # Local error: |E_t - E_{t-1}|, defined 0 at first sample
    local_err = similar(H)
    local_err[1] = 0.0
    @inbounds for i in 2:length(H)
        local_err[i] = abs(H[i] - H[i-1])
    end
    # Global error percent: |(E_t - E_0) / E_0| * 100
    global_err_pct = similar(H)
    E0 = H[1]
    @inbounds for i in eachindex(H)
        global_err_pct[i] = abs((H[i] - E0) / E0) * 100
    end
    write_dat(joinpath(DATA_DIR, "two_body_energy_errors.dat"),
        "local error |E_t - E_{t-1}| and global error percent |(E_t - E_0)/E_0| * 100",
        ["t", "local_err", "global_err_percent"],
        [energy.t, local_err, global_err_pct])

    @info "statistics" local_max=energy.statistics.local_error_max global_pct_max=energy.statistics.global_error_percent_max
end

main()
