# Agent 5 — Periodic orbit search & continuation for the 4-body 2+/2− Weber Hamiltonian.
#
# Strategy:
#   - T-brake shooting: an orbit is symmetric-periodic if a brake config (all p=0)
#     returns to a brake at t=T/2 (then T-symmetry mirrors second half).
#   - Try four candidate families: (1) breathing alternating square,
#     (2) rotating-square offset, (3) dimer-dimer orbital, (4) collinear ABAB.
#   - For each candidate compute period T (first brake or pair-distance return),
#     energy drift, and (for promising ones) finite-difference monodromy +
#     Floquet multipliers.
#
# Run from any working dir:
#     julia --project=. research/FourBodyTwoPlusTwoMinus/05_periodic_orbits/search.jl
#
# Outputs:
#   - found_orbits.csv (next to this file)
#   - figures/<id>.pdf for promising candidates

using LinearAlgebra
using Printf
using Plots

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator

const REPO    = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SHARED  = joinpath(REPO, "research", "FourBodyTwoPlusTwoMinus", "shared")
const OUTDIR  = @__DIR__
const FIGDIR  = joinpath(OUTDIR, "figures")
isdir(FIGDIR) || mkpath(FIGDIR)

include(joinpath(SHARED, "ic_generators.jl"))
include(joinpath(SHARED, "metrics.jl"))
include(joinpath(SHARED, "run_survey.jl"))
using .SharedSurvey

const N    = 4
const D    = 2
const DOF  = N * D
const SYS  = HamiltonianSystem(N, D)

# ----------------------------------------------------------------------------
# Direct integration (without going through SharedSurvey, so we keep dt control)
# ----------------------------------------------------------------------------
function integrate(q0, p0, masses, charges; tmax = 20.0, dt = 1e-3, c = 1.0,
                   bounce_r = 0.0)
    prob = HamiltonianProblem(
        SYS, (0.0, tmax), q0, p0;
        masses = masses, charges = charges, c = c, dt = dt,
    )
    return bounce_r > 0 ?
        solve(prob, SymmetricProjectionIntegrator(); callbacks = CollisionBounce(bounce_r)) :
        solve(prob, SymmetricProjectionIntegrator())
end

# ----------------------------------------------------------------------------
# Brake-return detector. Returns (i_return, T) or (0, NaN).
# ----------------------------------------------------------------------------
function brake_return(sol; ptol = 5e-2, skip = 50)
    isempty(sol.p) && return (0, NaN)
    pmag(p) = maximum(abs, p)
    n = length(sol.t)
    n < skip + 5 && return (0, NaN)
    # find first local minimum of |p|_∞ after `skip` that drops below tol
    best_i = 0
    best_v = Inf
    for i = (skip+1):(n-1)
        v = pmag(sol.p[i])
        if v < ptol && v < pmag(sol.p[i-1]) && v < pmag(sol.p[i+1])
            if v < best_v
                best_v = v
                best_i = i
                break
            end
        end
    end
    best_i == 0 && return (0, NaN)
    return (best_i, 2 * (sol.t[best_i] - sol.t[1]))
end

# ----------------------------------------------------------------------------
# Configuration return detector — pair-distance vector returns close to IC
# ----------------------------------------------------------------------------
function config_return(sol; rtol = 5e-2, skip = 50)
    n = length(sol.t)
    n < skip + 5 && return (0, NaN)
    d0 = pair_distances(sol.q[1], N, D)
    nd = norm(d0)
    best_i = 0
    best_v = Inf
    for i = (skip+1):(n-1)
        di = pair_distances(sol.q[i], N, D)
        v  = norm(di .- d0) / nd
        if v < rtol
            if v < best_v
                best_v = v; best_i = i
            else
                break
            end
        end
    end
    best_i == 0 && return (0, NaN)
    return (best_i, sol.t[best_i] - sol.t[1])
end

# ----------------------------------------------------------------------------
# Boundedness check
# ----------------------------------------------------------------------------
function bounded(sol)
    sol.retcode == :Success || return false
    diam0 = maximum(pair_distances(sol.q[1], N, D))
    threshold = 10 * diam0
    for q in sol.q
        maximum(pair_distances(q, N, D)) > threshold && return false
    end
    return true
end

function energy_drift_pct(sol)
    en = compute_energy_timeseries(sol)
    return en.statistics.global_error_percent_max
end

# ----------------------------------------------------------------------------
# Finite-difference monodromy
# ----------------------------------------------------------------------------
function monodromy(q0, p0, masses, charges, T; eps = 1e-5, dt = 1e-3, c = 1.0)
    nsteps = max(50, round(Int, T / dt))
    dt_use = T / nsteps
    z0 = vcat(q0, p0)
    flow(z) = begin
        s = integrate(z[1:DOF], z[DOF+1:end], masses, charges;
                      tmax = T, dt = dt_use, c = c)
        s.retcode == :Success || return fill(NaN, 2 * DOF)
        return vcat(s.q[end], s.p[end])
    end
    z_ref = flow(z0)
    any(isnan, z_ref) && return fill(NaN, 2 * DOF, 2 * DOF)
    M = zeros(2 * DOF, 2 * DOF)
    for k = 1:(2*DOF)
        zp = copy(z0); zp[k] += eps
        zm = copy(z0); zm[k] -= eps
        zp_T = flow(zp); zm_T = flow(zm)
        if any(isnan, zp_T) || any(isnan, zm_T)
            return fill(NaN, 2 * DOF, 2 * DOF)
        end
        M[:, k] = (zp_T .- zm_T) ./ (2 * eps)
    end
    return M
end

# ----------------------------------------------------------------------------
# Result row helper
# ----------------------------------------------------------------------------
mutable struct OrbitResult
    id::String
    family::String
    q0::Vector{Float64}
    p0::Vector{Float64}
    period::Float64
    energy::Float64
    max_floquet::Float64
    bound::Bool
    drift_pct::Float64
    notes::String
end

results = OrbitResult[]

function record!(id, family, q0, p0, sol, T; do_floquet = false, c = 1.0, masses = nothing, charges = nothing)
    bnd = bounded(sol)
    drift = bnd ? energy_drift_pct(sol) : NaN
    en = compute_energy_timeseries(sol)
    E0 = en.total_energy[1]
    floq = NaN
    if do_floquet && isfinite(T) && bnd && drift < 5.0
        try
            M = monodromy(q0, p0, masses, charges, T; c = c)
            if !any(isnan, M)
                ev = eigvals(M)
                floq = maximum(abs, ev)
            end
        catch e
            @warn "monodromy failed for $id" exception=e
        end
    end
    push!(results, OrbitResult(id, family, copy(q0), copy(p0),
        isnan(T) ? NaN : T, E0, floq, bnd, isnan(drift) ? NaN : drift, ""))
    @printf("  %-30s  ret=%-8s  T=%-8s  E0=%+.4f  drift=%-7s  bound=%s  |λ|max=%s\n",
        id, string(sol.retcode),
        isnan(T) ? "  -   " : @sprintf("%6.3f", T),
        E0,
        isnan(drift) ? "   -   " : @sprintf("%6.3f%%", drift),
        bnd, isnan(floq) ? "  -" : @sprintf("%.3f", floq))
end

function plot_traj!(id, sol)
    sol.retcode == :Success || return
    X = [unflatten(q, N, D) for q in sol.q]
    plt = plot(legend = :outerright, aspect_ratio = :equal,
               title = "trajectory $id", xlabel = "x", ylabel = "y")
    colors = [:red, :red, :blue, :blue]
    styles = [:solid, :dash, :solid, :dash]
    for i = 1:N
        xs = [X[k][i, 1] for k in eachindex(X)]
        ys = [X[k][i, 2] for k in eachindex(X)]
        plot!(plt, xs, ys, label = "p$i ($(i<=2 ? "+" : "-"))",
              color = colors[i], linestyle = styles[i])
    end
    savefig(plt, joinpath(FIGDIR, "$(id).pdf"))
end

# ----------------------------------------------------------------------------
# Family 1: Breathing alternating square
# ----------------------------------------------------------------------------
println("\n=== Family 1: breathing alternating square ===")
# Pure brake (zero momenta) → L=0 radial collapse → diagonal collision. Singular.
# Add a small *outward* radial velocity to make this a brake-at-min IC instead;
# alternatively combine with rotation to lift L from 0.
for s in (0.5, 1.0, 2.0)
    q0, _, masses, charges = alternating_square(side = s, velocity_mode = :static)
    sol = integrate(q0, zeros(DOF), masses, charges; tmax = 5.0, dt = 5e-4)
    i_ret, T = brake_return(sol; ptol = 5e-2, skip = 50)
    id = @sprintf("F1a_pure_brake_s%.2f", s)
    record!(id, "breathing_square_L0", q0, zeros(DOF), sol, T)

    # add small outward breathing velocity (still L=0, still singular but
    # potentially survives one period if the collapse is slow enough)
    for vrad in (0.1, 0.3, 0.4, 0.5, 0.6, 0.7)
        p_out = zeros(N, D)
        X = unflatten(q0, N, D)
        for i = 1:N
            r = X[i, :]; nrm = norm(r)
            p_out[i, :] = masses[i] * vrad * r / nrm
        end
        p0 = vec(reshape(p_out', :))
        # ^^ row-major: x1,y1,x2,y2,...
        p0_flat = Float64[]
        for i = 1:N, k = 1:D
            push!(p0_flat, p_out[i, k])
        end
        sol = integrate(q0, p0_flat, masses, charges; tmax = 10.0, dt = 5e-4)
        i_ret, T = brake_return(sol; ptol = 5e-2, skip = 100)
        id = @sprintf("F1b_outbreath_s%.2f_v%.2f", s, vrad)
        record!(id, "breathing_square_L0", q0, p0_flat, sol, T;
                do_floquet = !isnan(T) && bounded(sol),
                masses = masses, charges = charges)
        bounded(sol) && plot_traj!(id, sol)
    end
end

# ----------------------------------------------------------------------------
# Family 2: Rotating square + offset
# ----------------------------------------------------------------------------
println("\n=== Family 2: rotating-square near orbits ===")
# The square layout used by the generator is:
#   1=(s/2,s/2)+ , 2=(-s/2,-s/2)+ , 3=(-s/2,s/2)- , 4=(s/2,-s/2)- .
# 1,2 are diagonally opposite (++), 3,4 diagonally opposite (--). Each
# particle's distance from center is R=s/√2; edges (unlike) at distance s,
# diagonals (like) at distance s√2 = 2R. Force balance gives the same
# ω²(R) = (2√2−1)/(4R³) of Agent 4 (with our R = s/√2).
let s = 1.0
    R       = s / sqrt(2)
    ω_kepl  = sqrt((2*sqrt(2) - 1) / (4 * R^3))   # m=1
    masses  = fill(1.0, N)
    charges = [1.0, 1.0, -1.0, -1.0]
    # Use the generator to get q0 (zero velocities), then build p0 by hand.
    q0_static, _, _, _ = alternating_square(side = s, velocity_mode = :static)
    X = unflatten(q0_static, N, D)
    P = zeros(N, D)
    for i = 1:N
        r = X[i, :]
        P[i, :] = masses[i] * ω_kepl * [-r[2], r[1]]
    end
    p0_static = Float64[]
    for i = 1:N, k = 1:D
        push!(p0_static, P[i, k])
    end

    for δ in (0.0, 0.01, 0.05, -0.01, -0.05)
        q = copy(q0_static); p = copy(p0_static)
        # radial perturbation (breathing)
        for i = 1:N
            ix = (i - 1) * D + 1; iy = ix + 1
            nrm = hypot(q[ix], q[iy])
            q[ix] += δ * q[ix] / nrm
            q[iy] += δ * q[iy] / nrm
        end
        sol = integrate(q, p, masses, charges; tmax = 20.0, dt = 5e-5)
        i_ret, T = config_return(sol; rtol = 8e-2, skip = 10000)
        id = @sprintf("F2_rotsq_d%+.3f", δ)
        record!(id, "rotating_square", q, p, sol, T;
                do_floquet = !isnan(T) && bounded(sol),
                masses = masses, charges = charges)
        bounded(sol) && plot_traj!(id, sol)
    end
end

# ----------------------------------------------------------------------------
# Family 3: Dimer-dimer orbital
# ----------------------------------------------------------------------------
println("\n=== Family 3: dimer-dimer orbital ===")
for v_inter in (0.05, 0.1, 0.2, 0.3)
    q0, p0, masses, charges = two_dimers(
        dyad_length = 0.4, separation = 3.0,
        intra_fraction = 0.4, inter_velocity = v_inter,
    )
    sol = integrate(q0, p0, masses, charges; tmax = 50.0, dt = 2e-4)
    i_ret, T = config_return(sol; rtol = 8e-2, skip = 5000)
    id = @sprintf("F3_dimer_v%.2f", v_inter)
    record!(id, "dimer_dimer", q0, p0, sol, T;
            do_floquet = !isnan(T) && bounded(sol),
            masses = masses, charges = charges)
    bounded(sol) && plot_traj!(id, sol)
end

# ----------------------------------------------------------------------------
# Family 4: Collinear ABAB
# ----------------------------------------------------------------------------
println("\n=== Family 4: collinear ABAB ===")
for sp in (0.8, 1.0, 1.2)
    q0, p0, masses, charges = linear_abab(spacing = sp, transverse_kick = 0.0)
    sol = integrate(q0, p0, masses, charges; tmax = 20.0, dt = 2e-4)
    i_ret, T = brake_return(sol; ptol = 5e-2, skip = 500)
    id = @sprintf("F4_abab_s%.2f", sp)
    record!(id, "collinear_abab", q0, p0, sol, T;
            do_floquet = !isnan(T) && bounded(sol),
            masses = masses, charges = charges)
    bounded(sol) && plot_traj!(id, sol)
end

# also small transverse kick
for sp in (1.0,)
    for kick in (0.02, 0.05)
        q0, p0, masses, charges = linear_abab(spacing = sp, transverse_kick = kick)
        sol = integrate(q0, p0, masses, charges; tmax = 20.0, dt = 2e-4)
        i_ret, T = brake_return(sol; ptol = 5e-2, skip = 500)
        id = @sprintf("F4_abab_s%.2f_k%.2f", sp, kick)
        record!(id, "collinear_abab_kick", q0, p0, sol, T;
                do_floquet = false)
        bounded(sol) && plot_traj!(id, sol)
    end
end

# ----------------------------------------------------------------------------
# CSV output
# ----------------------------------------------------------------------------
csv_path = joinpath(OUTDIR, "found_orbits.csv")
open(csv_path, "w") do io
    println(io, "id,family,period,energy,max_floquet_mag,bound,drift_pct,q0,p0")
    for r in results
        q0s = join(r.q0, ";")
        p0s = join(r.p0, ";")
        @printf(io, "%s,%s,%.6f,%.6f,%.6f,%s,%.4f,%s,%s\n",
            r.id, r.family,
            isnan(r.period) ? -1.0 : r.period,
            r.energy,
            isnan(r.max_floquet) ? -1.0 : r.max_floquet,
            r.bound,
            isnan(r.drift_pct) ? -1.0 : r.drift_pct,
            q0s, p0s)
    end
end
println("\nWrote $(length(results)) results to $csv_path")
