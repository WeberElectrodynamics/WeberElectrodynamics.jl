"""
Agent 12 -- High-resolution Poincaré section analysis for the 4-body 2+/2-
Weber Hamiltonian.

Two systems studied:
  1. Rhombus KAM basin: grid around (a=1.5, b=1.45, eta=0.75, rotating)
  2. Symmetric double-orbiter (planar): grid around (r_pp=4, R=4, orb=1.3)

Section crossings are detected by sign changes with linear interpolation.
Classification: torus (closed curve), chaotic (area-filling), sparse, escape.
"""

using WeberElectrodynamics
using LinearAlgebra
using Printf
using Statistics

const HERE = @__DIR__
const SHARED_4B = normpath(joinpath(HERE, "..", "..", "FourBodyTwoPlusTwoMinus", "shared"))

include(joinpath(SHARED_4B, "ic_generators.jl"))
include(joinpath(SHARED_4B, "metrics.jl"))
include(joinpath(SHARED_4B, "run_survey.jl"))
using .SharedSurvey

const N = 4
const D = 2

# ============================================================================
# Jacobi vectors (equal masses): r+ = x2-x1, r- = x4-x3,
# R = (x3+x4)/2 - (x1+x2)/2
# ============================================================================
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

# ============================================================================
# Pair distance helpers
# ============================================================================
function pair_dist(q, i, j)
    X = reshape(q, D, N)
    return norm(X[:, i] .- X[:, j])
end

function pair_rdot(q, p, i, j)
    X = reshape(q, D, N)
    P = reshape(p, D, N)
    dx = X[:, i] .- X[:, j]
    dv = P[:, i] .- P[:, j]  # masses=1
    r = norm(dx)
    return (dx ⋅ dv) / max(r, eps())
end

# ============================================================================
# Angular section helpers: theta = atan(dy, dx) for pair (i,j)
# We use sin(theta) as the section function (crosses zero when theta = 0 or pi).
# For rotating orbits this gives ~1 crossing per orbit.
# ============================================================================
function pair_angle_y(q, i, j)
    X = reshape(q, D, N)
    return X[2, j] - X[2, i]   # y-component of rij
end

# ============================================================================
# Section detector (generic): g(q,p) crosses zero from below (g<0 -> g>=0)
# with positive derivative. Returns interpolated crossings.
# ============================================================================
function collect_crossings(sol, g_func, obs_func)
    crossings = NamedTuple{(:t, :x, :y), Tuple{Float64,Float64,Float64}}[]
    nsteps = length(sol.t)
    nsteps < 2 && return crossings
    g_prev = g_func(sol.q[1], sol.p[1])
    for k in 2:nsteps
        g_cur = g_func(sol.q[k], sol.p[k])
        if g_prev < 0.0 && g_cur >= 0.0
            alpha = -g_prev / (g_cur - g_prev)
            t_cross = sol.t[k-1] + alpha * (sol.t[k] - sol.t[k-1])
            q_cross = (1 - alpha) .* sol.q[k-1] .+ alpha .* sol.q[k]
            p_cross = (1 - alpha) .* sol.p[k-1] .+ alpha .* sol.p[k]
            ox, oy = obs_func(q_cross, p_cross)
            push!(crossings, (t=t_cross, x=ox, y=oy))
        end
        g_prev = g_cur
    end
    return crossings
end

# ============================================================================
# Classification heuristic
# ============================================================================
function classify_crossings(crossings, retcode; drift_pct=0.0)
    if retcode != :Success && drift_pct > 5.0
        return :escape
    end
    n = length(crossings)
    n < 10 && return :sparse
    xs = [c.x for c in crossings]
    ys = [c.y for c in crossings]
    rx = maximum(xs) - minimum(xs); rx <= 0 && (rx = eps())
    ry = maximum(ys) - minimum(ys); ry <= 0 && (ry = eps())
    sx = std(xs); sy = std(ys)
    fill = (sx * sy) / (rx * ry + eps())
    # Additionally: check if points lie near a 1D curve by looking at
    # the ratio of the two principal component variances
    if n >= 15
        mx = mean(xs); my = mean(ys)
        cxx = sum((x - mx)^2 for x in xs) / n
        cyy = sum((y - my)^2 for y in ys) / n
        cxy = sum((xs[i] - mx) * (ys[i] - my) for i in 1:n) / n
        tr = cxx + cyy
        det_cov = cxx * cyy - cxy^2
        disc = max(tr^2 - 4*det_cov, 0.0)
        l1 = (tr + sqrt(disc)) / 2
        l2 = (tr - sqrt(disc)) / 2
        ratio = l2 / max(l1, eps())
        # Curve-like if l2 << l1 (ratio < ~0.05)
        if ratio < 0.05
            return :torus
        elseif ratio < 0.15 && fill < 0.15
            return :torus
        end
    end
    return fill < 0.08 ? :torus : :chaotic
end

# ============================================================================
# RHOMBUS BASIN: Section definition
# Section: y-component of r_plus (= x2-x1) crosses zero from below.
# This triggers ~once per orbital period for rotating configurations.
# Record (r34, dr34/dt) on the section.
# ============================================================================
function run_rhombus_poincare(; a_vals=[1.3, 1.4, 1.5, 1.6, 1.7],
                                b_vals=[1.2, 1.3, 1.4, 1.45, 1.5, 1.6],
                                eta_vals=[0.6, 0.7, 0.75, 0.8],
                                tmax=200.0, dt=5e-4)
    results = NamedTuple[]
    section_dir = joinpath(HERE, "section_data")

    for a in a_vals, b in b_vals, eta in eta_vals
        label = @sprintf("rhombus_a%.2f_b%.2f_eta%.2f", a, b, eta)
        @info "Running $label"
        q0, p0, masses, charges = rhombus(a=a, b=b, energy_fraction=eta,
                                           velocity_mode=:rotating)

        sol = nothing
        try
            sol = SharedSurvey.run(q0, p0, masses, charges;
                                    tmax=tmax, dt=dt, c=1.0, dims=D,
                                    bounce_r=0.02)
        catch err
            @warn "Failed" label err
            push!(results, (config=:rhombus, a=a, b=b, eta=eta,
                           retcode=:Exception, t_final=0.0,
                           n_crossings=0, classification=:escape,
                           drift_pct=NaN))
            continue
        end

        en = compute_energy_timeseries(sol)
        drift = en.statistics.global_error_percent_max

        # Section: y-component of (x2-x1) crosses zero from below
        g_func(q, p) = pair_angle_y(q, 1, 2)
        obs_func(q, p) = (pair_dist(q, 3, 4), pair_rdot(q, p, 3, 4))

        crossings = collect_crossings(sol, g_func, obs_func)
        cls = classify_crossings(crossings, sol.retcode; drift_pct=drift)

        push!(results, (config=:rhombus, a=a, b=b, eta=eta,
                       retcode=sol.retcode, t_final=sol.t[end],
                       n_crossings=length(crossings), classification=cls,
                       drift_pct=drift))

        # Save crossing data
        if !isempty(crossings)
            open(joinpath(section_dir, "$label.csv"), "w") do io
                println(io, "t,r34,dr34_dt")
                for c in crossings
                    @printf(io, "%.8e,%.8e,%.8e\n", c.t, c.x, c.y)
                end
            end
        end

        @info "  => $(sol.retcode) t=$(sol.t[end]) crossings=$(length(crossings)) class=$cls drift=$(drift)%"
    end
    return results
end

# ============================================================================
# DOUBLE-ORBITER: IC generator
# Symmetric double-orbiter: ++ pair at separation r_pp on x-axis,
# -- pair orbiting at radius R with speed orb.
# ============================================================================
function double_orbiter_ic(; r_pp=4.0, R=4.0, orb=1.3, m=1.0, q=1.0)
    # Positive pair: particles 1,2 on x-axis
    # Negative pair: particles 3,4 on y-axis (symmetric)
    X = zeros(4, 2)
    X[1, :] = [-r_pp/2, 0.0]
    X[2, :] = [r_pp/2, 0.0]
    X[3, :] = [0.0, R]
    X[4, :] = [0.0, -R]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    # Orbital velocities for the negative pair (tangential)
    P = zeros(4, 2)
    P[3, :] = m * orb * [1.0, 0.0]
    P[4, :] = m * orb * [-1.0, 0.0]
    # Zero COM and total momentum
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

function run_double_orbiter_poincare(; orb_vals=[1.0, 1.1, 1.2, 1.3, 1.4, 1.5],
                                       rpp_vals=[3.5, 4.0, 4.5],
                                       tmax=200.0, dt=5e-4)
    results = NamedTuple[]
    section_dir = joinpath(HERE, "section_data")

    for orb in orb_vals, rpp in rpp_vals
        label = @sprintf("double_orb_rpp%.1f_orb%.1f", rpp, orb)
        @info "Running $label"
        q0, p0, masses, charges = double_orbiter_ic(r_pp=rpp, R=4.0, orb=orb)

        # Section: y-component of (x4-x3) crosses zero from below.
        # The -- pair orbits, so this triggers once per orbit.
        sol = nothing
        try
            sol = SharedSurvey.run(q0, p0, masses, charges;
                                    tmax=tmax, dt=dt, c=1.0, dims=D,
                                    bounce_r=0.02)
        catch err
            @warn "Failed" label err
            push!(results, (config=:double_orbiter, r_pp=rpp, orb=orb,
                           retcode=:Exception, t_final=0.0,
                           n_crossings=0, classification=:escape,
                           drift_pct=NaN))
            continue
        end

        en = compute_energy_timeseries(sol)
        drift = en.statistics.global_error_percent_max

        g_func(q, p) = pair_angle_y(q, 3, 4)
        obs_func(q, p) = (pair_dist(q, 2, 4), pair_rdot(q, p, 2, 4))

        crossings = collect_crossings(sol, g_func, obs_func)
        cls = classify_crossings(crossings, sol.retcode; drift_pct=drift)

        push!(results, (config=:double_orbiter, r_pp=rpp, orb=orb,
                       retcode=sol.retcode, t_final=sol.t[end],
                       n_crossings=length(crossings), classification=cls,
                       drift_pct=drift))

        if !isempty(crossings)
            open(joinpath(section_dir, "$label.csv"), "w") do io
                println(io, "t,r24,dr24_dt")
                for c in crossings
                    @printf(io, "%.8e,%.8e,%.8e\n", c.t, c.x, c.y)
                end
            end
        end

        @info "  => $(sol.retcode) t=$(sol.t[end]) crossings=$(length(crossings)) class=$cls drift=$(drift)%"
    end
    return results
end

# ============================================================================
# Write combined summary
# ============================================================================
function write_summary(rhombus_results, orbiter_results)
    open(joinpath(HERE, "poincare_summary.csv"), "w") do io
        println(io, "config,param1,param2,param3,retcode,t_final,n_crossings,classification,drift_pct")
        for r in rhombus_results
            @printf(io, "%s,%.2f,%.2f,%.2f,%s,%.4f,%d,%s,%.6e\n",
                    r.config, r.a, r.b, r.eta, r.retcode, r.t_final,
                    r.n_crossings, r.classification, r.drift_pct)
        end
        for r in orbiter_results
            @printf(io, "%s,%.2f,%.2f,NaN,%s,%.4f,%d,%s,%.6e\n",
                    r.config, r.r_pp, r.orb, r.retcode, r.t_final,
                    r.n_crossings, r.classification, r.drift_pct)
        end
    end
end

# ============================================================================
# Main
# ============================================================================
function main()
    println("=" ^ 70)
    println("Agent 12 — High-resolution Poincaré section analysis")
    println("=" ^ 70)

    println("\n--- Phase 1: Rhombus basin (5x6x4 = 120 ICs) ---")
    rhombus_results = run_rhombus_poincare()

    println("\n--- Phase 2: Double-orbiter (6x3 = 18 ICs) ---")
    orbiter_results = run_double_orbiter_poincare()

    write_summary(rhombus_results, orbiter_results)

    # Print summary statistics
    println("\n" * "=" ^ 70)
    println("RHOMBUS BASIN SUMMARY")
    println("=" ^ 70)
    n_success = count(r -> r.retcode == :Success, rhombus_results)
    n_total = length(rhombus_results)
    println("  Success: $n_success / $n_total")
    for cls in (:torus, :chaotic, :sparse, :escape)
        nc = count(r -> r.classification == cls, rhombus_results)
        println("  $cls: $nc")
    end

    # Best candidates (most crossings)
    good = filter(r -> r.n_crossings >= 10, rhombus_results)
    sort!(good, by=r->r.n_crossings, rev=true)
    println("\nTop rhombus candidates (>= 10 crossings):")
    for r in good[1:min(20, length(good))]
        @printf("  a=%.2f b=%.2f eta=%.2f: %d crossings, class=%s, drift=%.2e%%\n",
                r.a, r.b, r.eta, r.n_crossings, r.classification, r.drift_pct)
    end

    println("\n" * "=" ^ 70)
    println("DOUBLE-ORBITER SUMMARY")
    println("=" ^ 70)
    n_success = count(r -> r.retcode == :Success, orbiter_results)
    n_total = length(orbiter_results)
    println("  Success: $n_success / $n_total")
    for cls in (:torus, :chaotic, :sparse, :escape)
        nc = count(r -> r.classification == cls, orbiter_results)
        println("  $cls: $nc")
    end

    good_orb = filter(r -> r.n_crossings >= 5, orbiter_results)
    sort!(good_orb, by=r->r.n_crossings, rev=true)
    println("\nTop orbiter candidates (>= 5 crossings):")
    for r in good_orb[1:min(10, length(good_orb))]
        @printf("  rpp=%.1f orb=%.1f: %d crossings, class=%s, drift=%.2e%%\n",
                r.r_pp, r.orb, r.n_crossings, r.classification, r.drift_pct)
    end

    return (rhombus=rhombus_results, orbiter=orbiter_results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
