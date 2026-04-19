"""
Agent 8 — Sub-Weber-radius dynamics for the 4-body 2+/2- Weber Hamiltonian.

Convention: m=q=c=1, so for like-charge pairs the critical radius is
    rho = q_i q_j / (mu_ij c^2) = 1 / (1/2) = 2.
Particles 1,2 are positive; 3,4 are negative. The 6 pair indices in the order
returned by SharedMetrics.pair_distances are:
    (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)

Experiments
-----------
A: like-pair (+,+) sub-critical "nucleus" with (-,-) far-away orbiters.
B: two simultaneous sub-critical nuclei  (+,+) and (-,-).
C: tight (+,-) Coulomb dimer plus far (+,-) orbiter.
D: 5x5 collision survival map over (r0_like, p0_like) for the (+,+) pair,
   spectator (-,-) at r=1 placed far away. Output: phase_diagram.png.

Honest reporting: per the Frauenfelder-Weber 2024 result, sub-rho spirals
with l != 0 are non-regularizable. Many of these runs are expected to fail
with retcode != :Success. We log the failures as data.

Budget: 4 named experiments + 25 grid cells + a few warmups; tmax <= 50.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using WeberElectrodynamics
using LinearAlgebra
using Printf
using Plots

include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "..", "shared", "metrics.jl"))

const N = 4
const D = 2
const C_LIGHT = 1.0
const Q = 1.0
const M = 1.0
const RHO_LIKE = Q * Q / ((M / 2) * C_LIGHT^2)   # = 2.0

const OUTDIR = @__DIR__
isdir(joinpath(OUTDIR, "figures")) || mkdir(joinpath(OUTDIR, "figures"))

# ----------------------------------------------------------------------
# IC builders specific to Agent 8
# ----------------------------------------------------------------------

"""
Build a 2D 4-body IC (flat vectors). All particles in the xy-plane.
positions: 4x2 matrix; momenta: 4x2 matrix.
"""
function _flat(X, P)
    return (vec(permutedims(X)), vec(permutedims(P)))
end

"""
Experiment A: like (+,+) sub-critical pair as nucleus.
- (+,+) at (-r/2,0) and (+r/2,0), small head-on momenta +/- p_head along x.
- (-,-) far away at (-R/2, R0) and (+R/2, R0) with small tangential kicks.
"""
function ic_A(; r_like = 0.5, p_head = 0.2, R = 5.0, p_orb = 0.1)
    X = [
        -r_like/2  0.0;
        r_like/2  0.0;
        -R/2       R;
        R/2       R;
    ]
    P = [
        p_head    0.0;
        -p_head   0.0;
        0.0       p_orb;
        0.0      -p_orb;
    ]
    masses = [M, M, M, M]
    charges = [Q, Q, -Q, -Q]
    # com / total p removal
    com = (masses' * X) ./ sum(masses)
    X .-= com
    P .-= sum(P, dims = 1) ./ N
    return _flat(X, P)..., masses, charges
end

"""
Experiment B: two sub-critical like nuclei, well separated along x.
"""
function ic_B(; r_like = 0.5, p_head = 0.15, R = 6.0)
    X = [
        -R/2 - r_like/2  0.0;
        -R/2 + r_like/2  0.0;
        R/2 - r_like/2  0.0;
        R/2 + r_like/2  0.0;
    ]
    P = [
        p_head  0.0;
        -p_head 0.0;
        p_head  0.0;
        -p_head 0.0;
    ]
    masses = [M, M, M, M]; charges = [Q, Q, -Q, -Q]
    com = (masses' * X) ./ sum(masses); X .-= com
    P .-= sum(P, dims = 1) ./ N
    return _flat(X, P)..., masses, charges
end

"""
Experiment C: tight (+,-) inner dimer + far (+,-) orbiter.
Pair (1,3) acts as the tight dimer at r0_in. Pair (2,4) is the far dimer.
"""
function ic_C(; r0_in = 0.3, R = 4.0, v_orb = 0.2)
    X = [
        -R/2  r0_in/2;       # +1  (inner dimer top)
        R/2  0.0;            # +2  (far)
        -R/2 -r0_in/2;       # -1  (inner dimer bottom)
        R/2  1.0;            # -2  (far, offset)
    ]
    # Tangential momenta for inner dimer (1,3) about its own COM
    mu = M / 2
    v_in = sqrt(0.5 * Q * Q / (mu * r0_in))   # half of circular Coulomb estimate
    P = [
        M * v_in        0.0;
        0.0       M * v_orb;
        -M * v_in       0.0;
        0.0      -M * v_orb;
    ]
    masses = [M, M, M, M]; charges = [Q, Q, -Q, -Q]
    com = (masses' * X) ./ sum(masses); X .-= com
    P .-= sum(P, dims = 1) ./ N
    return _flat(X, P)..., masses, charges
end

"""
Experiment D grid cell: like (+,+) pair at r0 with head-on momenta +/- p0,
spectator (-,-) at r=1 fixed, far away at (0, R_far).
"""
function ic_D(; r0_like = 0.5, p0 = 0.2, R_far = 8.0)
    X = [
        -r0_like/2  0.0;
        r0_like/2  0.0;
        -0.5      R_far;
        0.5       R_far;
    ]
    P = [
        p0    0.0;
        -p0   0.0;
        0.0   0.0;
        0.0   0.0;
    ]
    masses = [M, M, M, M]; charges = [Q, Q, -Q, -Q]
    com = (masses' * X) ./ sum(masses); X .-= com
    P .-= sum(P, dims = 1) ./ N
    return _flat(X, P)..., masses, charges
end

# ----------------------------------------------------------------------
# Diagnostics
# ----------------------------------------------------------------------

function pair_distance_timeseries(sol)
    nsteps = length(sol.t)
    out = zeros(nsteps, 6)
    for k in 1:nsteps
        out[k, :] = pair_distances(sol.q[k], N, D)
    end
    return out
end

function summarize_run(label, sol; bounce_r=0.0)
    en = compute_energy_timeseries(sol)
    drift = en.statistics.global_error_percent_max
    pdist = pair_distance_timeseries(sol)
    rmin = minimum(pdist, dims = 1)
    rmax = maximum(pdist, dims = 1)
    bound = bound_indicator(sol; escape_radius = 25.0)
    @printf("[%s] retcode=%s  steps=%d  t_end=%.3f  E_drift_max=%.3g%%  bound=%s\n",
            label, sol.retcode, length(sol.t), sol.t[end], drift, bound)
    @printf("    pair r_min: ")
    for v in rmin; @printf("%.3f ", v); end; println()
    @printf("    pair r_max: ")
    for v in rmax; @printf("%.3f ", v); end; println()
    return (label=label, retcode=sol.retcode, drift=drift, bound=bound,
            pdist=pdist, t=sol.t)
end

function classify_grid_cell(sol; like_pair_idx=1, escape_thresh=20.0)
    if sol.retcode != :Success
        return :Failure
    end
    pdist = pair_distance_timeseries(sol)
    # Unbound test: any pair distance grows beyond escape_thresh
    if maximum(pdist) > escape_thresh
        return :Unbound
    end
    return :Success
end

# ----------------------------------------------------------------------
# Run experiments
# ----------------------------------------------------------------------

const TMAX = 50.0
const DT = 5e-4
const BOUNCE_R = 0.05

println("# Sub-Weber-radius dynamics — Agent 8")
println("# rho_like = $(RHO_LIKE)   tmax=$(TMAX)  dt=$(DT)  bounce_r=$(BOUNCE_R)\n")

results = Dict{String,Any}()

# Experiment A
let
    q0,p0,m,ch = ic_A(r_like=0.5, p_head=0.2, R=5.0, p_orb=0.1)
    sol = SharedSurvey.run(q0,p0,m,ch; tmax=TMAX, dt=DT, c=C_LIGHT, dims=D, bounce_r=BOUNCE_R)
    results["A"] = summarize_run("A: ++ nucleus + -- orbiters", sol)
end

# Experiment B
let
    q0,p0,m,ch = ic_B(r_like=0.5, p_head=0.15, R=6.0)
    sol = SharedSurvey.run(q0,p0,m,ch; tmax=TMAX, dt=DT, c=C_LIGHT, dims=D, bounce_r=BOUNCE_R)
    results["B"] = summarize_run("B: two sub-critical nuclei", sol)
end

# Experiment C
let
    q0,p0,m,ch = ic_C(r0_in=0.3, R=4.0, v_orb=0.2)
    sol = SharedSurvey.run(q0,p0,m,ch; tmax=TMAX, dt=DT, c=C_LIGHT, dims=D, bounce_r=0.02)
    results["C"] = summarize_run("C: tight (+,-) dimer + far (+,-)", sol)
end

# Plot pair distances for A,B,C
for label in ("A","B","C")
    r = results[label]
    plt = plot(r.t, r.pdist, label=["(1,2)" "(1,3)" "(1,4)" "(2,3)" "(2,4)" "(3,4)"],
               xlabel="t", ylabel="pair distance", title="Experiment $label   retcode=$(r.retcode)",
               lw=1.2, legend=:outerright)
    hline!(plt, [RHO_LIKE], ls=:dash, c=:red, label="rho_like=$(RHO_LIKE)")
    savefig(plt, joinpath(OUTDIR, "figures", "pairs_$(label).png"))
end

# Experiment D — 5x5 phase grid
println("\n# Experiment D — 5x5 phase grid (r0_like, p0)")
const R0_GRID = collect(range(0.2, stop=1.5, length=5))   # all sub-rho (rho=2)
const P0_GRID = collect(range(0.05, stop=0.6, length=5))
labels = fill(:Failure, length(R0_GRID), length(P0_GRID))
drifts = fill(NaN, size(labels))

for (i, r0) in enumerate(R0_GRID), (j, p0_) in enumerate(P0_GRID)
    q0,p0v,m,ch = ic_D(r0_like=r0, p0=p0_, R_far=8.0)
    local sol
    try
        sol = SharedSurvey.run(q0,p0v,m,ch; tmax=TMAX, dt=DT, c=C_LIGHT, dims=D,
                               bounce_r=BOUNCE_R)
    catch e
        @printf("  cell (r0=%.2f, p0=%.2f) threw: %s\n", r0, p0_, sprint(showerror, e))
        labels[i,j] = :Failure
        continue
    end
    cls = classify_grid_cell(sol)
    labels[i,j] = cls
    if sol.retcode == :Success
        en = compute_energy_timeseries(sol)
        drifts[i,j] = en.statistics.global_error_percent_max
    end
    @printf("  cell (r0=%.2f, p0=%.2f) -> %s  (retcode=%s, drift=%.2g%%)\n",
            r0, p0_, cls, sol.retcode, drifts[i,j])
end

# Phase diagram heatmap
code_map = Dict(:Success => 2, :Unbound => 1, :Failure => 0)
Z = [code_map[labels[i,j]] for i in 1:size(labels,1), j in 1:size(labels,2)]
heat = heatmap(P0_GRID, R0_GRID, Z;
               xlabel = "initial |p| of (+,+) pair",
               ylabel = "initial r of (+,+) pair",
               title  = "Sub-Weber-radius survival map (2+/2-)\nblack=Failure  blue=Unbound  yellow=Success",
               clims  = (0, 2),
               c      = :viridis,
               aspect_ratio = :auto)
hline!(heat, [RHO_LIKE], ls=:dash, c=:red, label="rho_like=$(RHO_LIKE)")
# Annotate cells
for i in 1:size(labels,1), j in 1:size(labels,2)
    annotate!(heat, P0_GRID[j], R0_GRID[i], text(string(labels[i,j])[1:1], 8, :white))
end
savefig(heat, joinpath(OUTDIR, "phase_diagram.png"))

# Save grid summary CSV
open(joinpath(OUTDIR, "grid_summary.csv"), "w") do io
    println(io, "r0_like,p0,label,drift_pct")
    for (i, r0) in enumerate(R0_GRID), (j, p0_) in enumerate(P0_GRID)
        @printf(io, "%.4f,%.4f,%s,%.6g\n", r0, p0_, labels[i,j], drifts[i,j])
    end
end

println("\nDone. Outputs in $(OUTDIR).")
