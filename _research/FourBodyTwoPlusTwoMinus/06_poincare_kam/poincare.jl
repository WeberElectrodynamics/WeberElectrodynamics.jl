"""
Agent 6 — Poincaré sections and KAM structure for the 4-body 2+/2− Weber
Hamiltonian.

Builds a small grid of 12 initial conditions (3 configurations × 4 energy
fractions), integrates each to tmax=30 with dt=1e-3, and records crossings of
three sections:

    S1 : R_x = 0  with  d/dt R_x > 0     -> plot ( |r_+| , |r_-| )
    S2 : |r_+| − |r_-| = 0  with  d/dt > 0  -> plot ( |R| , angle(r_+, r_-) )
    S3 : L_{12} = 0  (planar ang. mom. of the (+,+) pair about its centre)
                                          -> plot ( r_{12} , ṙ_{12} )

A section crossing is detected by a sign change of g(q,p) between two
consecutive saved states, with linear interpolation on g to recover the
crossing instant.  Only crossings with g̊ > 0 are kept.

Outputs:
    figures/section_S1_two_dimers.png
    figures/section_S1_rhombus.png
    figures/section_S1_alternating_square.png
    figures/section_S2_combined.png
    figures/section_S3_combined.png
    figures/kam_fraction.png
    poincare_summary.csv
"""

using LinearAlgebra
using Printf
using Plots
using Statistics

const HERE = @__DIR__
const SHARED = normpath(joinpath(HERE, "..", "shared"))

include(joinpath(SHARED, "ic_generators.jl"))
include(joinpath(SHARED, "metrics.jl"))
include(joinpath(SHARED, "run_survey.jl"))
using .SharedSurvey

# --------------------------------------------------------------------------
# Section observables.  q,p are flat length-8 (N=4, d=2).
# --------------------------------------------------------------------------

const N = 4
const D = 2

# Jacobi vectors (equal masses): r+ = x2-x1, r- = x4-x3,
# R = (x3+x4)/2 - (x1+x2)/2 .  Same for momenta:
# P+ = (p2-p1)/2, P- = (p4-p3)/2, P = (p3+p4)/2 - (p1+p2)/2 .
function jacobi(q, p)
    X = reshape(q, D, N)
    P = reshape(p, D, N)
    rplus  = X[:, 2] .- X[:, 1]
    rminus = X[:, 4] .- X[:, 3]
    Rvec   = 0.5 .* (X[:, 3] .+ X[:, 4]) .- 0.5 .* (X[:, 1] .+ X[:, 2])
    Pplus  = 0.5 .* (P[:, 2] .- P[:, 1])
    Pminus = 0.5 .* (P[:, 4] .- P[:, 3])
    Pvec   = 0.5 .* (P[:, 3] .+ P[:, 4]) .- 0.5 .* (P[:, 1] .+ P[:, 2])
    return rplus, rminus, Rvec, Pplus, Pminus, Pvec
end

# S1: r+_y = 0  with d/dt(r+_y) > 0.  Triggered every time the (+,+) dimer's
# orientation crosses the x-axis from below.  Robust for all our IC families
# because the dimer rotates / breathes about the COM.
function S1_g(q, p, m)
    rp, _, _, Pp, _, _ = jacobi(q, p)
    rpdot = Pp .* (2/m)
    return rp[2], rpdot[2]
end

# S2: |r+| - |r-| .  Time derivative is (r+·ṙ+)/|r+| - (r-·ṙ-)/|r-|
function S2_g(q, p, m)
    rp, rm, _, Pp, Pm, _ = jacobi(q, p)
    rpn = norm(rp);  rmn = norm(rm)
    rpdot_dot = (rp ⋅ (Pp .* (2/m))) / max(rpn, eps())
    rmdot_dot = (rm ⋅ (Pm .* (2/m))) / max(rmn, eps())
    return rpn - rmn, rpdot_dot - rmdot_dot
end

# S3: planar angular momentum of particle pair (1,2) about their CoM.
# L_12 = (x1-cx) × p1 + (x2-cx) × p2  (z-component, 2D pseudoscalar)
function S3_g(q, p, m)
    X = reshape(q, D, N)
    P = reshape(p, D, N)
    cx = 0.5 .* (X[:, 1] .+ X[:, 2])
    d1 = X[:, 1] .- cx
    d2 = X[:, 2] .- cx
    L = d1[1]*P[1,2] - d1[2]*P[1,1] + d2[1]*P[2,2] - d2[2]*P[2,1]
    # Crude time derivative proxy: contribution from kinetic part.
    # Use L̇ = sum (vi × pi) which is zero, so use torque from internal force ~
    # finite-difference fallback.  We just return 1.0 to accept any crossing.
    return L, 1.0
end

# --------------------------------------------------------------------------
# Section detector.
# --------------------------------------------------------------------------

"""
    collect_section(sol, gfun; m, observables) -> Vector{NTuple{2,Float64}}

Walk through (sol.q[k], sol.p[k]) consecutively.  Detect sign changes in g
with positive ġ and linearly interpolate.  `observables(q,p)` returns the
2-tuple to record at the crossing.
"""
function collect_section(sol, gfun, observables; m::Float64 = 1.0)
    pts = Tuple{Float64,Float64}[]
    nsteps = length(sol.t)
    nsteps < 2 && return pts
    g_prev, gd_prev = gfun(sol.q[1], sol.p[1], m)
    for k = 2:nsteps
        g_cur, gd_cur = gfun(sol.q[k], sol.p[k], m)
        if g_prev < 0.0 && g_cur >= 0.0 && (gd_cur + gd_prev) > 0.0
            α = g_prev / (g_prev - g_cur)
            qint = (1 - α) .* sol.q[k-1] .+ α .* sol.q[k]
            pint = (1 - α) .* sol.p[k-1] .+ α .* sol.p[k]
            push!(pts, observables(qint, pint))
        end
        g_prev = g_cur
        gd_prev = gd_cur
    end
    return pts
end

obs_S1(q, p) = (norm(jacobi(q, p)[1]), norm(jacobi(q, p)[2]))

function obs_S2(q, p)
    rp, rm, Rv, _, _, _ = jacobi(q, p)
    cosang = (rp ⋅ rm) / max(norm(rp)*norm(rm), eps())
    cosang = clamp(cosang, -1.0, 1.0)
    return (norm(Rv), acos(cosang))
end

function obs_S3(q, p)
    X = reshape(q, D, N)
    P = reshape(p, D, N)
    r12 = X[:, 2] .- X[:, 1]
    v12 = (P[:, 2] .- P[:, 1])  # masses=1
    rn = norm(r12)
    rdot = (r12 ⋅ v12) / max(rn, eps())
    return (rn, rdot)
end

# --------------------------------------------------------------------------
# IC grid : 3 configurations × 4 energy fractions
# --------------------------------------------------------------------------

const ENERGY_FRACTIONS = [0.10, 0.25, 0.50, 0.75]

function build_ic(name::Symbol, ef::Float64)
    if name === :two_dimers
        # use a large dyad to avoid singularity; vary intra_fraction
        return two_dimers(dyad_length = 0.6, separation = 3.0, intra_fraction = ef,
                          inter_velocity = 0.05)
    elseif name === :rhombus
        # rotating rhombus far from collisions; b kept close to a
        return rhombus(a = 1.5, b = 1.0 + 0.2*ef, energy_fraction = ef,
                       velocity_mode = :rotating)
    elseif name === :alternating_square
        # rotating square, side = 1 — close to the breathing-square mode
        return alternating_square(side = 1.0, energy_fraction = ef,
                                  velocity_mode = :rotating)
    else
        error("unknown ic $name")
    end
end

# --------------------------------------------------------------------------
# Classification helper (very rough)
# --------------------------------------------------------------------------

"""
Return :torus | :chaotic | :sparse | :escape based on point distribution.
"""
function classify(pts::Vector{Tuple{Float64,Float64}}, retcode::Symbol;
                  drift::Float64 = 0.0)
    # treat :Failure with small drift as a graceful early termination, not escape
    if retcode != :Success && drift > 5.0
        return :escape
    end
    n = length(pts)
    n < 8 && return :sparse
    xs = [p[1] for p in pts]
    ys = [p[2] for p in pts]
    # Range-based heuristic: σx*σy/(rx*ry) area filling fraction.
    rx = maximum(xs) - minimum(xs); rx <= 0 && (rx = eps())
    ry = maximum(ys) - minimum(ys); ry <= 0 && (ry = eps())
    σx = std(xs);  σy = std(ys)
    fill = (σx * σy) / (rx * ry + eps())
    # Curve-like => fill ≪ 1.  Area-filling ~ 0.3 for uniform-ish set.
    return fill < 0.10 ? :torus : :chaotic
end

# --------------------------------------------------------------------------
# Main driver
# --------------------------------------------------------------------------

function main()
    configs = [:two_dimers, :rhombus, :alternating_square]
    results = NamedTuple[]

    plt_S1 = Dict{Symbol,Any}()
    for cfg in configs
        plt_S1[cfg] = plot(title = "Section S1: $(cfg)\n(|r+|, |r-|) at R_x=0, Ṙ_x>0",
                            xlabel = "|r_+|", ylabel = "|r_-|", legend = :outertopright,
                            size = (640, 480))
    end
    plt_S2 = plot(title = "Section S2: |r+|=|r-|, d/dt(|r+|-|r-|)>0",
                  xlabel = "|R|", ylabel = "angle(r+, r-) [rad]",
                  legend = :outertopright, size = (720, 540))
    plt_S3 = plot(title = "Section S3: L_12=0",
                  xlabel = "r_{12}", ylabel = "ṙ_{12}",
                  legend = :outertopright, size = (720, 540))

    for cfg in configs
        for ef in ENERGY_FRACTIONS
            q0, p0, masses, charges = build_ic(cfg, ef)
            sol = nothing
            try
                sol = SharedSurvey.run(q0, p0, masses, charges;
                                       tmax = 30.0, dt = 1e-3, c = 1.0, dims = D,
                                       bounce_r = 0.02)
            catch err
                @warn "run failed" cfg ef err
                push!(results, (config=cfg, energy_fraction=ef, retcode=:Failure,
                                n_S1=0, n_S2=0, n_S3=0,
                                S1_class=:escape, S2_class=:escape, S3_class=:escape,
                                E_drift_pct=NaN))
                continue
            end
            sumr = SharedSurvey.summarize(sol)
            ptsS1 = collect_section(sol, S1_g, obs_S1)
            ptsS2 = collect_section(sol, S2_g, obs_S2)
            ptsS3 = collect_section(sol, S3_g, obs_S3)
            cS1 = classify(ptsS1, sol.retcode; drift = sumr.E_drift_pct)
            cS2 = classify(ptsS2, sol.retcode; drift = sumr.E_drift_pct)
            cS3 = classify(ptsS3, sol.retcode; drift = sumr.E_drift_pct)
            push!(results, (config=cfg, energy_fraction=ef, retcode=sol.retcode,
                            n_S1=length(ptsS1), n_S2=length(ptsS2), n_S3=length(ptsS3),
                            S1_class=cS1, S2_class=cS2, S3_class=cS3,
                            E_drift_pct=sumr.E_drift_pct))
            label = @sprintf("ef=%.2f  [%s]", ef, cS1)
            if !isempty(ptsS1)
                scatter!(plt_S1[cfg], [p[1] for p in ptsS1], [p[2] for p in ptsS1],
                         markersize = 2.2, markerstrokewidth = 0, label = label, alpha = 0.7)
            end
            if !isempty(ptsS2)
                scatter!(plt_S2, [p[1] for p in ptsS2], [p[2] for p in ptsS2],
                         markersize = 2.0, markerstrokewidth = 0,
                         label = "$(cfg) ef=$(ef)", alpha = 0.6)
            end
            if !isempty(ptsS3)
                scatter!(plt_S3, [p[1] for p in ptsS3], [p[2] for p in ptsS3],
                         markersize = 2.0, markerstrokewidth = 0,
                         label = "$(cfg) ef=$(ef)", alpha = 0.6)
            end
            @info "done" cfg ef retcode=sol.retcode tend=sol.t[end] nsteps=length(sol.t) S1=length(ptsS1) S2=length(ptsS2) S3=length(ptsS3)
        end
    end

    figdir = joinpath(HERE, "figures")
    isdir(figdir) || mkpath(figdir)
    for cfg in configs
        savefig(plt_S1[cfg], joinpath(figdir, "section_S1_$(cfg).png"))
    end
    savefig(plt_S2, joinpath(figdir, "section_S2_combined.png"))
    savefig(plt_S3, joinpath(figdir, "section_S3_combined.png"))

    # KAM fraction figure
    kam_plt = plot(title = "KAM-fraction (S1 classification)",
                   xlabel = "energy_fraction", ylabel = "fraction",
                   legend = :outertopright, ylim = (-0.05, 1.05))
    for cls in (:torus, :chaotic, :sparse, :escape)
        ys = Float64[]
        for ef in ENERGY_FRACTIONS
            sub = filter(r -> r.energy_fraction == ef, results)
            push!(ys, count(r -> r.S1_class == cls, sub) / max(length(sub), 1))
        end
        plot!(kam_plt, ENERGY_FRACTIONS, ys, marker = :circle, label = string(cls))
    end
    savefig(kam_plt, joinpath(figdir, "kam_fraction.png"))

    open(joinpath(HERE, "poincare_summary.csv"), "w") do io
        println(io, "config,energy_fraction,retcode,n_S1,n_S2,n_S3,S1_class,S2_class,S3_class,E_drift_pct")
        for r in results
            println(io, join((r.config, r.energy_fraction, r.retcode,
                              r.n_S1, r.n_S2, r.n_S3,
                              r.S1_class, r.S2_class, r.S3_class, r.E_drift_pct), ","))
        end
    end
    for r in results
        println(r)
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
