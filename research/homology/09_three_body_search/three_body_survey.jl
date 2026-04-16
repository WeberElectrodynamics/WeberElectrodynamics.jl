"""
Systematic 3-body Weber bound orbit search.
Agent 09 — first-ever systematic survey of 3-body configurations.

Runs ~200 integrations across:
  1. 2+/1- (equilateral, collinear, planetary atom, isosceles)
  2. 1+/2- (same geometries, opposite charge assignment)
  3. 3 same-sign (+++ or ---)
  4. Asymmetric masses/charges (helium-like)
  5. Perturbations of bound orbits

Output: survey_results.csv
"""

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using LinearAlgebra
using Random
using Printf
using Dates

include(joinpath(@__DIR__, "ic_generators_3body.jl"))

# ==========================================================================
# Solver wrapper
# ==========================================================================
const _SYS_CACHE = Dict{Tuple{Int,Int},Any}()
function _get_sys(n, d)
    key = (n, d)
    haskey(_SYS_CACHE, key) && return _SYS_CACHE[key]
    sys = WeberSystem(n, d)
    _SYS_CACHE[key] = sys
    return sys
end

function run_3body(q0, p0, masses, charges;
    tmax=50.0, dt=1e-3, c=1.0, dims=2, bounce_r=0.0,
    zollner_enabled=false, zollner_a=0.0)
    sys = _get_sys(3, dims)
    reg = RegularizationOptions(collision_bounce_radius=bounce_r)
    zol = ZollnerOptions(enabled=zollner_enabled, a=zollner_a)
    prob = WeberProblem(sys, (0.0, tmax), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt,
        regularization=reg, zollner=zol)
    sol = solve(prob, SymmetricProjectionIntegrator())
    return sol
end

# ==========================================================================
# Diagnostics
# ==========================================================================
function unflatten(v, n, d)
    reshape(copy(v), d, n)'
end

function pair_distances_all(sol)
    n = 3; d = length(sol.q[1]) ÷ n
    dists = Vector{Vector{Float64}}()
    for k in 1:length(sol.t)
        X = unflatten(sol.q[k], n, d)
        ds = Float64[]
        for i = 1:n, j = (i+1):n
            push!(ds, norm(view(X, i, :) .- view(X, j, :)))
        end
        push!(dists, ds)
    end
    return dists
end

function bound_check(sol; escape_factor=10.0)
    sol.retcode == :Success || return false
    dists = pair_distances_all(sol)
    isempty(dists) && return false
    d0_max = maximum(dists[1])
    threshold = escape_factor * d0_max
    for ds in dists
        maximum(ds) > threshold && return false
    end
    return true
end

function max_pair_distance(sol)
    dists = pair_distances_all(sol)
    isempty(dists) && return Inf
    return maximum(maximum(ds) for ds in dists)
end

function min_pair_distance(sol)
    dists = pair_distances_all(sol)
    isempty(dists) && return 0.0
    return minimum(minimum(ds) for ds in dists)
end

function energy_drift(sol)
    try
        en = compute_energy_timeseries(sol)
        return en.statistics.global_error_percent_max
    catch
        return NaN
    end
end

function brake_return_3body(sol; ptol=5e-2, skip=50)
    n = length(sol.t)
    n < skip + 5 && return (0, NaN)
    for i = (skip+1):(n-1)
        v = maximum(abs, sol.p[i])
        if v < ptol && v < maximum(abs, sol.p[i-1]) && v < maximum(abs, sol.p[i+1])
            return (i, 2 * (sol.t[i] - sol.t[1]))
        end
    end
    return (0, NaN)
end

# ==========================================================================
# Result storage
# ==========================================================================
mutable struct RunResult
    id::Int
    category::String
    config::String
    charges::String
    masses::String
    c::Float64
    eta::Float64
    side_or_sep::Float64
    bounce_r::Float64
    tmax::Float64
    dt::Float64
    retcode::Symbol
    n_steps::Int
    t_final::Float64
    E_drift_pct::Float64
    bound::Bool
    max_r::Float64
    min_r::Float64
    brake_idx::Int
    brake_period::Float64
    notes::String
end

function result_header()
    return "id,category,config,charges,masses,c,eta,side_or_sep,bounce_r,tmax,dt,retcode,n_steps,t_final,E_drift_pct,bound,max_r,min_r,brake_idx,brake_period,notes"
end

function result_csv(r::RunResult)
    return @sprintf("%d,%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.1f,%.1e,%s,%d,%.4f,%.6f,%s,%.6f,%.6f,%d,%.6f,%s",
        r.id, r.category, r.config, r.charges, r.masses,
        r.c, r.eta, r.side_or_sep, r.bounce_r, r.tmax, r.dt,
        r.retcode, r.n_steps, r.t_final, r.E_drift_pct,
        r.bound, r.max_r, r.min_r, r.brake_idx, r.brake_period, r.notes)
end

# ==========================================================================
# Main survey
# ==========================================================================
function run_survey()
    results = RunResult[]
    run_id = 0
    csv_path = joinpath(@__DIR__, "survey_results.csv")

    function record!(sol, category, config, charges_str, masses_str, c, eta, side_or_sep, bounce_r, tmax, dt; notes="")
        run_id += 1
        rc = sol.retcode
        ns = length(sol.t)
        tf = sol.t[end]
        ed = energy_drift(sol)
        bd = bound_check(sol)
        mr = max_pair_distance(sol)
        mnr = min_pair_distance(sol)
        bi, bp = brake_return_3body(sol)
        r = RunResult(run_id, category, config, charges_str, masses_str,
            c, eta, side_or_sep, bounce_r, tmax, dt,
            rc, ns, tf, isnan(ed) ? -1.0 : ed, bd, mr, mnr, bi, bp, notes)
        push!(results, r)
        status = bd ? "BOUND" : (rc == :Success ? "unbound" : string(rc))
        @printf("  [%3d] %-20s %-25s ret=%-8s E_drift=%.4f%% bound=%-5s max_r=%.3f\n",
            run_id, category, config, rc, isnan(ed) ? -1.0 : ed, bd, mr)
        return r
    end

    println("="^80)
    println("THREE-BODY WEBER BOUND ORBIT SURVEY")
    println("Started: ", now())
    println("="^80)

    # ======================================================================
    # CATEGORY 1: 2+/1- configurations (~60 runs)
    # ======================================================================
    println("\n--- CATEGORY 1: 2+/1- ---")

    # 1a. Equilateral triangle, rotating
    for eta in [0.1, 0.25, 0.5, 0.75]
        for side in [1.0, 2.0, 3.0]
            for c in [1.0, 4.0]
                q0, p0, m, ch = equilateral_triangle(
                    side=side, charges=[1.0, 1.0, -1.0], masses=[1.0, 1.0, 1.0],
                    energy_fraction=eta, velocity_mode=:rotating)
                sol = run_3body(q0, p0, m, ch; tmax=50.0, dt=1e-3, c=c)
                record!(sol, "2+1-", "equi_rot_s$(side)_eta$(eta)", "[+1+1-1]", "[1 1 1]",
                    c, eta, side, 0.0, 50.0, 1e-3)
            end
        end
    end

    # 1b. Collinear +−+ (Euler-like)
    for spacing in [1.0, 2.0, 3.0]
        for kick in [0.0, 0.1, 0.3]
            q0, p0, m, ch = collinear_symmetric(
                spacing=spacing, charges=[1.0, -1.0, 1.0], masses=[1.0, 1.0, 1.0],
                transverse_kick=kick)
            sol = run_3body(q0, p0, m, ch; tmax=50.0, dt=1e-3, c=1.0)
            record!(sol, "2+1-", "collin_+-+_d$(spacing)_k$(kick)", "[+1-1+1]", "[1 1 1]",
                1.0, 0.0, spacing, 0.0, 50.0, 1e-3)
        end
    end

    # 1c. Planetary atom (++- nucleus + orbiter) — most promising
    for eta in [0.3, 0.5, 0.6, 0.7, 0.8, 0.9]
        for R in [0.5, 1.0, 2.0]
            for c in [1.0, 4.0]
                q0, p0, m, ch = planetary_atom_3(
                    nucleus_sep=0.05, orbit_radius=R,
                    charges=[1.0, 1.0, -1.0], masses=[1.0, 1.0, 1.0],
                    eta=eta, c=c)
                sol = run_3body(q0, p0, m, ch; tmax=100.0, dt=1e-4, c=c, bounce_r=0.02)
                record!(sol, "2+1-", "planet_R$(R)_eta$(eta)", "[+1+1-1]", "[1 1 1]",
                    c, eta, R, 0.02, 100.0, 1e-4)
            end
        end
    end

    # ======================================================================
    # CATEGORY 2: 1+/2- configurations (~40 runs)
    # ======================================================================
    println("\n--- CATEGORY 2: 1+/2- ---")

    # 2a. Equilateral triangle, rotating
    for eta in [0.1, 0.25, 0.5, 0.75]
        for side in [1.0, 2.0]
            for c in [1.0, 4.0]
                q0, p0, m, ch = equilateral_triangle(
                    side=side, charges=[-1.0, -1.0, 1.0], masses=[1.0, 1.0, 1.0],
                    energy_fraction=eta, velocity_mode=:rotating)
                sol = run_3body(q0, p0, m, ch; tmax=50.0, dt=1e-3, c=c)
                record!(sol, "1+2-", "equi_rot_s$(side)_eta$(eta)", "[-1-1+1]", "[1 1 1]",
                    c, eta, side, 0.0, 50.0, 1e-3)
            end
        end
    end

    # 2b. Planetary atom with (−−) nucleus, (+) orbiter
    for eta in [0.5, 0.7, 0.8, 0.9]
        for R in [0.5, 1.0, 2.0]
            q0, p0, m, ch = planetary_atom_3(
                nucleus_sep=0.05, orbit_radius=R,
                charges=[-1.0, -1.0, 1.0], masses=[1.0, 1.0, 1.0],
                eta=eta, c=4.0)
            sol = run_3body(q0, p0, m, ch; tmax=100.0, dt=1e-4, c=4.0, bounce_r=0.02)
            record!(sol, "1+2-", "planet_R$(R)_eta$(eta)", "[-1-1+1]", "[1 1 1]",
                4.0, eta, R, 0.02, 100.0, 1e-4)
        end
    end

    # ======================================================================
    # CATEGORY 3: Same-sign (+++ or ---) (~30 runs)
    # ======================================================================
    println("\n--- CATEGORY 3: Same-sign ---")

    # 3a. Sub-critical equilateral triangle with rotation
    for side in [0.05, 0.08, 0.10]
        for eta in [0.1, 0.25, 0.5]
            for c in [1.0, 4.0]
                q0, p0, m, ch = same_sign_triangle(
                    side=side, charge=1.0, mass=1.0, c=c,
                    energy_fraction=eta, velocity_mode=:rotating)
                sol = run_3body(q0, p0, m, ch; tmax=5.0, dt=1e-4, c=c, bounce_r=0.02)
                record!(sol, "+++", "equi_rot_s$(side)_eta$(eta)", "[+1+1+1]", "[1 1 1]",
                    c, eta, side, 0.02, 5.0, 1e-4)
            end
        end
    end

    # 3b. Collinear breathing at sub-critical spacing
    for side in [0.05, 0.10]
        for c in [1.0, 4.0]
            q0, p0, m, ch = collinear_symmetric(
                spacing=side, charges=[1.0, 1.0, 1.0], masses=[1.0, 1.0, 1.0],
                transverse_kick=0.0, velocity_mode=:static)
            sol = run_3body(q0, p0, m, ch; tmax=2.0, dt=1e-5, c=c, bounce_r=0.02)
            record!(sol, "+++", "collin_breath_d$(side)", "[+1+1+1]", "[1 1 1]",
                c, 0.0, side, 0.02, 2.0, 1e-5)
        end
    end

    # 3c. Large triangle with collision bounce (super-critical distances)
    for side in [1.0, 2.0]
        q0, p0, m, ch = same_sign_triangle(
            side=side, charge=1.0, mass=1.0, c=1.0,
            energy_fraction=0.5, velocity_mode=:rotating)
        sol = run_3body(q0, p0, m, ch; tmax=20.0, dt=1e-3, c=1.0, bounce_r=0.02)
        record!(sol, "+++", "equi_rot_super_s$(side)", "[+1+1+1]", "[1 1 1]",
            1.0, 0.5, side, 0.02, 20.0, 1e-3)
    end

    # ======================================================================
    # CATEGORY 4: Asymmetric (helium-like) (~40 runs)
    # ======================================================================
    println("\n--- CATEGORY 4: Asymmetric ---")

    # 4a. Heavy nucleus (+2) + two light electrons (-1,-1)
    for eta in [0.1, 0.25, 0.5, 0.75]
        for R in [1.0, 2.0]
            for angle in [Float64(π), Float64(2π/3)]
                q0, p0, m, ch = helium_like(
                    orbit_radius=R, angle_sep=angle,
                    charges=[2.0, -1.0, -1.0], masses=[10.0, 1.0, 1.0],
                    energy_fraction=eta)
                sol = run_3body(q0, p0, m, ch; tmax=50.0, dt=1e-3, c=1.0)
                angle_str = angle ≈ π ? "pi" : "2pi3"
                record!(sol, "asym", "He_q2m10_R$(R)_a$(angle_str)_eta$(eta)",
                    "[+2-1-1]", "[10 1 1]", 1.0, eta, R, 0.0, 50.0, 1e-3)
            end
        end
    end

    # 4b. (+1,+1,-2) "inverted helium" with light nucleus
    for eta in [0.25, 0.5, 0.75]
        for R in [1.0, 2.0]
            q0, p0, m, ch = helium_like(
                orbit_radius=R, angle_sep=Float64(π),
                charges=[-2.0, 1.0, 1.0], masses=[1.0, 1.0, 1.0],
                energy_fraction=eta)
            sol = run_3body(q0, p0, m, ch; tmax=50.0, dt=1e-3, c=1.0)
            record!(sol, "asym", "inv_He_q-2_R$(R)_eta$(eta)",
                "[-2+1+1]", "[1 1 1]", 1.0, eta, R, 0.0, 50.0, 1e-3)
        end
    end

    # 4c. Heavy nucleus (+2) with c=4 and smaller orbits
    for eta in [0.5, 0.7, 0.8]
        for R in [0.5, 1.0]
            q0, p0, m, ch = helium_like(
                orbit_radius=R, angle_sep=Float64(π),
                charges=[2.0, -1.0, -1.0], masses=[10.0, 1.0, 1.0],
                energy_fraction=eta)
            sol = run_3body(q0, p0, m, ch; tmax=50.0, dt=1e-3, c=4.0)
            record!(sol, "asym", "He_c4_R$(R)_eta$(eta)",
                "[+2-1-1]", "[10 1 1]", 4.0, eta, R, 0.0, 50.0, 1e-3)
        end
    end

    # ======================================================================
    # CATEGORY 5: Perturbations of bound orbits (~30 runs)
    # Will be filled after identifying bound configs above
    # ======================================================================
    println("\n--- CATEGORY 5: Perturbation tests ---")

    # Find bound results to perturb
    bound_results = filter(r -> r.bound && r.E_drift_pct < 1.0, results)
    println("  Found $(length(bound_results)) bound orbits to perturb")

    # Take up to 10 best bound orbits (lowest energy drift)
    sort!(bound_results, by = r -> r.E_drift_pct)
    n_perturb = min(10, length(bound_results))

    for idx in 1:n_perturb
        br = bound_results[idx]
        # Re-generate the IC for this run
        # We'll use a generic approach: re-run the same config with perturbations
        # Parse config to regenerate — we store enough info
        for (delta, label) in [(0.01, "d001"), (0.05, "d005"), (0.1, "d01")]
            # Re-generate base IC from the config name
            base_q0, base_p0, base_m, base_ch = _regenerate_ic(br)
            if isnothing(base_q0)
                continue
            end
            q_pert, p_pert, m_pert, ch_pert = perturb_ic(
                base_q0, base_p0, base_m, base_ch;
                delta_q=delta, delta_p=delta, seed=42+idx)
            sol = run_3body(q_pert, p_pert, m_pert, ch_pert;
                tmax=br.tmax, dt=br.dt, c=br.c, bounce_r=br.bounce_r)
            record!(sol, "perturb", "$(br.config)_$(label)",
                br.charges, br.masses, br.c, br.eta, br.side_or_sep,
                br.bounce_r, br.tmax, br.dt;
                notes="delta=$(delta)_base=$(br.id)")
        end
    end

    # ======================================================================
    # Write CSV
    # ======================================================================
    println("\n" * "="^80)
    println("WRITING RESULTS")
    println("="^80)

    open(csv_path, "w") do f
        println(f, result_header())
        for r in results
            println(f, result_csv(r))
        end
    end

    # Summary statistics
    total = length(results)
    success = count(r -> r.retcode == :Success, results)
    bound = count(r -> r.bound, results)
    good_bound = count(r -> r.bound && r.E_drift_pct < 1.0, results)

    println("\nSUMMARY:")
    println("  Total runs: $total")
    println("  Success: $success / $total ($(round(100*success/total, digits=1))%)")
    println("  Bound: $bound / $total ($(round(100*bound/total, digits=1))%)")
    println("  Good bound (E_drift < 1%): $good_bound / $total")

    # By category
    println("\nBy category:")
    for cat in unique(r.category for r in results)
        cat_results = filter(r -> r.category == cat, results)
        ct = length(cat_results)
        cs = count(r -> r.retcode == :Success, cat_results)
        cb = count(r -> r.bound, cat_results)
        cg = count(r -> r.bound && r.E_drift_pct < 1.0, cat_results)
        println("  $cat: $ct runs, $cs success, $cb bound, $cg good-bound")
    end

    println("\nResults written to: $csv_path")
    println("Finished: ", now())

    return results
end

# Helper to regenerate ICs for perturbation tests
function _regenerate_ic(r::RunResult)
    try
        cfg = r.config
        if startswith(cfg, "planet_")
            # Parse R and eta from config name
            m_R = match(r"R([\d.]+)", cfg)
            m_eta = match(r"eta([\d.]+)", cfg)
            isnothing(m_R) && return (nothing, nothing, nothing, nothing)
            R = parse(Float64, m_R.captures[1])
            eta = parse(Float64, m_eta.captures[1])
            if r.charges == "[+1+1-1]"
                return planetary_atom_3(nucleus_sep=0.05, orbit_radius=R,
                    charges=[1.0, 1.0, -1.0], masses=[1.0, 1.0, 1.0],
                    eta=eta, c=r.c)
            elseif r.charges == "[-1-1+1]"
                return planetary_atom_3(nucleus_sep=0.05, orbit_radius=R,
                    charges=[-1.0, -1.0, 1.0], masses=[1.0, 1.0, 1.0],
                    eta=eta, c=r.c)
            end
        elseif startswith(cfg, "equi_rot_")
            m_s = match(r"s([\d.]+)", cfg)
            m_eta = match(r"eta([\d.]+)", cfg)
            isnothing(m_s) && return (nothing, nothing, nothing, nothing)
            s = parse(Float64, m_s.captures[1])
            eta = parse(Float64, m_eta.captures[1])
            ch = r.charges == "[+1+1-1]" ? [1.0, 1.0, -1.0] :
                 r.charges == "[-1-1+1]" ? [-1.0, -1.0, 1.0] :
                 [1.0, 1.0, 1.0]
            return equilateral_triangle(side=s, charges=ch,
                masses=[1.0, 1.0, 1.0], energy_fraction=eta,
                velocity_mode=:rotating)
        elseif contains(cfg, "He_")
            m_R = match(r"R([\d.]+)", cfg)
            m_eta = match(r"eta([\d.]+)", cfg)
            isnothing(m_R) && return (nothing, nothing, nothing, nothing)
            R = parse(Float64, m_R.captures[1])
            eta = parse(Float64, m_eta.captures[1])
            return helium_like(orbit_radius=R, angle_sep=Float64(π),
                charges=[2.0, -1.0, -1.0], masses=[10.0, 1.0, 1.0],
                energy_fraction=eta)
        end
    catch e
        @warn "Failed to regenerate IC for run $(r.id): $e"
    end
    return (nothing, nothing, nothing, nothing)
end

# Run!
results = run_survey()
