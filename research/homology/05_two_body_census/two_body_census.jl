"""
Exhaustive 2-body Weber orbit census.

Systematically surveys bound orbits across charge signs, mass ratios,
speed-of-light values, angular momenta, and energies.

Outputs census_results.csv with one row per integration.
"""

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using LinearAlgebra
using Printf
using Dates

# Cache compiled systems
const SYS2D = WeberSystem(2, 2)

# ---------------------------------------------------------------------------
# Initial condition constructors
# ---------------------------------------------------------------------------

"""
Build 2-body ICs for unlike charges (q1*q2 < 0) at apoapsis.
Returns (q0, p0, r_apoapsis) or nothing if infeasible.
"""
function unlike_charge_ics(q1, q2, m1, m2, E, L)
    mu = m1 * m2 / (m1 + m2)
    k = q1 * q2  # negative
    M = m1 + m2

    if L == 0.0
        # Head-on: E = k/r0 => r0 = k/E
        r0 = k / E
        r0 <= 0 && return nothing
        q0 = [-(m2 / M) * r0, 0.0, (m1 / M) * r0, 0.0]
        p0 = [0.0, 0.0, 0.0, 0.0]
        return (q0, p0, r0)
    end

    # Quadratic for u = 1/r at turning points
    disc = k^2 + 2 * E * L^2 / mu
    disc < 0 && return nothing
    sqrt_disc = sqrt(disc)
    L2_over_mu = L^2 / mu
    u_a = (-k - sqrt_disc) / L2_over_mu
    u_a <= 0 && return nothing
    r_a = 1.0 / u_a

    # Tangential speed at apoapsis
    p_perp = L / r_a  # = mu * v_perp

    q0 = [-(m2 / M) * r_a, 0.0, (m1 / M) * r_a, 0.0]
    p0 = [0.0, p_perp, 0.0, -p_perp]
    return (q0, p0, r_a)
end

"""
Build 2-body ICs for like charges in the sub-critical regime.
r0 < rho = q1*q2/(mu*c^2) required.
"""
function like_charge_ics(q1, q2, m1, m2, c, E, L)
    mu = m1 * m2 / (m1 + m2)
    k = q1 * q2  # positive
    M = m1 + m2
    rho = k / (mu * c^2)

    if L == 0.0
        E <= 0 && return nothing
        r0 = k / E
        r0 >= rho && return nothing
        q0 = [-(m2 / M) * r0, 0.0, (m1 / M) * r0, 0.0]
        p0 = [0.0, 0.0, 0.0, 0.0]
        return (q0, p0, r0, rho)
    end

    # E = L^2/(2mu r0^2) + k/r0 => E r0^2 - k r0 - L^2/(2mu) = 0
    disc = k^2 + 2 * E * L^2 / mu
    disc < 0 && return nothing
    sqrt_disc = sqrt(disc)
    r0_plus = (k + sqrt_disc) / (2 * E)
    r0_minus = (k - sqrt_disc) / (2 * E)

    # Pick valid sub-critical turning point
    r0 = nothing
    for candidate in [r0_plus, r0_minus]
        if candidate > 0 && candidate < rho
            if r0 === nothing || candidate > r0
                r0 = candidate
            end
        end
    end
    r0 === nothing && return nothing

    p_perp = L / r0
    q0 = [-(m2 / M) * r0, 0.0, (m1 / M) * r0, 0.0]
    p0 = [0.0, p_perp, 0.0, -p_perp]
    return (q0, p0, r0, rho)
end

# ---------------------------------------------------------------------------
# Integration and analysis
# ---------------------------------------------------------------------------

function run_and_analyze(q0, p0, q1_charge, q2_charge, m1, m2, c;
    tmax=100.0, dt=1e-3, bounce_r=0.0, use_reg=false)

    reg = if use_reg
        RegularizationOptions(enabled=true, collision_bounce_radius=bounce_r)
    elseif bounce_r > 0
        RegularizationOptions(collision_bounce_radius=bounce_r)
    else
        RegularizationOptions()
    end

    prob = WeberProblem(
        SYS2D, (0.0, tmax), q0, p0;
        masses=[m1, m2], charges=[q1_charge, q2_charge], c=c, dt=dt,
        regularization=reg
    )
    sol = solve(prob, SymmetricProjectionIntegrator())

    success = sol.retcode == :Success

    r0_init = norm([q0[3] - q0[1], q0[4] - q0[2]])
    max_sep = 0.0
    min_sep = Inf
    stride = max(1, length(sol.t) ÷ 1000)

    crossings = Float64[]
    prev_r = r0_init

    for k in 1:stride:length(sol.t)
        dx = sol.q[k][3] - sol.q[k][1]
        dy = sol.q[k][4] - sol.q[k][2]
        r = sqrt(dx^2 + dy^2)
        max_sep = max(max_sep, r)
        min_sep = min(min_sep, r)
        if prev_r < r0_init && r >= r0_init && k > stride
            push!(crossings, sol.t[k])
        end
        prev_r = r
    end

    drift_pct = NaN
    if success
        try
            en = compute_energy_timeseries(sol)
            drift_pct = en.statistics.global_error_percent_max
        catch
            drift_pct = NaN
        end
    end

    bound = success && max_sep < 10.0 * r0_init && (isnan(drift_pct) || drift_pct < 1.0)

    period = NaN
    if length(crossings) >= 2
        periods = diff(crossings)
        sort!(periods)
        n = length(periods)
        period = n % 2 == 1 ? periods[(n+1) ÷ 2] : (periods[n ÷ 2] + periods[n ÷ 2 + 1]) / 2
    end

    ecc = (max_sep - min_sep) / (max_sep + min_sep + 1e-30)
    orbit_type = classify_orbit(bound, success, ecc, q1_charge * q2_charge)

    return (bound=bound, period=period, orbit_type=orbit_type,
            drift_pct=drift_pct, max_sep=max_sep, min_sep=min_sep,
            retcode=sol.retcode)
end

function classify_orbit(bound, success, ecc, k_sign)
    !success && return "failed"
    !bound && return "unbound"
    if k_sign > 0
        return ecc < 0.05 ? "like_quasi_circular" : "like_oscillation"
    else
        ecc < 0.02 && return "circular"
        ecc < 0.3 && return "low_ecc_elliptic"
        ecc < 0.7 && return "moderate_ecc_elliptic"
        return "high_ecc_elliptic"
    end
end

function write_row(io, q1, q2, m1, m2, c, r0, L, E, res, label)
    @printf(io, "%.4f,%.4f,%.2f,%.2f,%.2f,%.8f,%.4f,%.6f,%s,%.6f,%s,%.6f,%.6f,%.6f,%s,%s\n",
            q1, q2, m1, m2, c, r0, L, E,
            res.bound, res.period, res.orbit_type,
            res.drift_pct, res.max_sep, res.min_sep, res.retcode, label)
end

function write_error_row(io, q1, q2, m1, m2, c, r0, L, E, label)
    @printf(io, "%.4f,%.4f,%.2f,%.2f,%.2f,%.8f,%.4f,%.6f,false,NaN,error,NaN,NaN,NaN,Error,%s\n",
            q1, q2, m1, m2, c, r0, L, E, label)
end

# ---------------------------------------------------------------------------
# Main survey
# ---------------------------------------------------------------------------

function run_census(outfile::String)
    io = open(outfile, "w")
    println(io, "q1,q2,m1,m2,c,r0,L,E,bound,period,orbit_type,drift_pct,max_sep,min_sep,retcode,case_label")

    run_count = 0
    t_start = time()

    # ===================================================================
    # SECTION 1: Unlike charges — hydrogen-like orbits
    # ===================================================================
    println("=== Section 1: Unlike charges (q1*q2 < 0) ===")

    for (m1, m2) in [(1.0, 1.0), (1.0, 2.0), (1.0, 10.0), (1.0, 100.0)]
        for c in [1.0, 2.0, 4.0, 10.0, 100.0]
            for L in [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]
                for E in [-0.01, -0.025, -0.05, -0.1, -0.25, -0.5, -1.0, -2.0]
                    q1, q2 = 1.0, -1.0
                    ics = unlike_charge_ics(q1, q2, m1, m2, E, L)
                    ics === nothing && continue
                    q0, p0, r0 = ics

                    # Reject orbits that are too tight for the integrator
                    r0 < 0.05 && continue

                    mu = m1 * m2 / (m1 + m2)
                    k_abs = abs(q1 * q2)
                    a_semi = k_abs / (2 * abs(E))
                    T_kepler = 2 * pi * a_semi^1.5 * sqrt(mu / k_abs)

                    # Scale integration time: enough periods but not too long
                    n_periods = 15
                    tmax = min(n_periods * T_kepler, 500.0)
                    tmax = max(tmax, 20.0)
                    # dt: at least 100 steps per period, but not too small
                    dt = min(T_kepler / 200, 5e-3)
                    dt = max(dt, 5e-5)

                    bounce_r = L == 0.0 ? max(0.02 * r0, 0.01) : 0.0

                    try
                        res = run_and_analyze(q0, p0, q1, q2, m1, m2, c;
                            tmax=tmax, dt=dt, bounce_r=bounce_r)
                        write_row(io, q1, q2, m1, m2, c, r0, L, E, res, "unlike")
                    catch ex
                        write_error_row(io, q1, q2, m1, m2, c, r0, L, E, "unlike")
                        println("  ERR unlike: m=($m1,$m2) c=$c L=$L E=$E : $ex")
                    end
                    run_count += 1
                    run_count % 50 == 0 && println("  [$run_count runs, $(round(time()-t_start,digits=1))s]")
                end
            end
        end
    end
    sec1 = run_count
    println("Section 1 done: $sec1 runs")

    # ===================================================================
    # SECTION 2: Like charges — sub-critical bound oscillations
    # ===================================================================
    println("\n=== Section 2: Like charges (q1*q2 > 0) ===")

    for (m1, m2) in [(1.0, 1.0), (1.0, 2.0), (1.0, 10.0), (1.0, 100.0)]
        for c in [1.0, 2.0, 4.0, 10.0, 100.0]
            mu = m1 * m2 / (m1 + m2)
            rho = 1.0 / (mu * c^2)
            E_min = mu * c^2

            # For like-charge sub-critical: r0 = 1/E, need r0 < rho
            # So E > E_min. Try a range of r0/rho fractions instead.
            r0_fractions = [0.9, 0.7, 0.5, 0.3, 0.1, 0.05, 0.01]

            for L in [0.0, 0.1, 0.5]
                for frac in r0_fractions
                    q1, q2 = 1.0, 1.0
                    if L == 0.0
                        # Radial: r0 = frac * rho, E = k/r0 = 1/(frac*rho)
                        r0 = frac * rho
                        E = 1.0 / r0  # k = q1*q2 = 1
                    else
                        # With angular momentum: specify r0, compute E
                        r0 = frac * rho
                        E = L^2 / (2 * mu * r0^2) + 1.0 / r0
                    end

                    # Skip if r0 is extremely tiny (numerical issues)
                    r0 < 1e-6 && continue

                    # Period estimate
                    T_est = 2 * sqrt(2) * r0 / c
                    T_est < 1e-8 && continue

                    n_osc = 20
                    tmax = min(n_osc * T_est, 100.0)
                    tmax = max(tmax, 2.0)
                    # Need very small dt for sub-critical oscillation
                    dt = min(T_est / 100, 1e-3)
                    dt = max(dt, 1e-6)

                    bounce_r = L == 0.0 ? max(0.02 * r0, 1e-6) : 0.0
                    use_reg = L > 0.0

                    try
                        res = run_and_analyze(q0_p0(q1, q2, m1, m2, r0, L, mu)...,
                            q1, q2, m1, m2, c;
                            tmax=tmax, dt=dt, bounce_r=bounce_r, use_reg=use_reg)
                        write_row(io, q1, q2, m1, m2, c, r0, L, E, res, "like")
                    catch ex
                        write_error_row(io, q1, q2, m1, m2, c, r0, L, E, "like")
                        println("  ERR like: m=($m1,$m2) c=$c L=$L frac=$frac : $ex")
                    end
                    run_count += 1
                    run_count % 50 == 0 && println("  [$run_count runs, $(round(time()-t_start,digits=1))s]")
                end
            end
        end
    end
    sec2 = run_count - sec1
    println("Section 2 done: $sec2 runs (total $run_count)")

    # ===================================================================
    # SECTION 3: Asymmetric charge magnitudes
    # ===================================================================
    println("\n=== Section 3: Asymmetric charges ===")

    asym_configs = [
        (1.0, -2.0, "asym_1_-2"),
        (2.0, -1.0, "asym_2_-1"),
        (1.0, -0.5, "asym_1_-0.5"),
        (0.5, -1.0, "asym_0.5_-1"),
    ]

    for (q1, q2, label) in asym_configs
        for (m1, m2) in [(1.0, 1.0), (1.0, 2.0), (1.0, 10.0)]
            for c in [2.0, 4.0, 10.0]
                for L in [0.0, 0.25, 0.5, 1.0, 2.0]
                    for E in [-0.01, -0.05, -0.1, -0.25, -0.5, -1.0]
                        ics = unlike_charge_ics(q1, q2, m1, m2, E, L)
                        ics === nothing && continue
                        q0, p0, r0 = ics
                        r0 < 0.05 && continue

                        mu = m1 * m2 / (m1 + m2)
                        k_abs = abs(q1 * q2)
                        a_semi = k_abs / (2 * abs(E))
                        T_kepler = 2 * pi * a_semi^1.5 * sqrt(mu / k_abs)
                        tmax = min(15 * T_kepler, 500.0)
                        tmax = max(tmax, 20.0)
                        dt = min(T_kepler / 200, 5e-3)
                        dt = max(dt, 5e-5)

                        bounce_r = L == 0.0 ? max(0.02 * r0, 0.01) : 0.0

                        try
                            res = run_and_analyze(q0, p0, q1, q2, m1, m2, c;
                                tmax=tmax, dt=dt, bounce_r=bounce_r)
                            write_row(io, q1, q2, m1, m2, c, r0, L, E, res, label)
                        catch ex
                            write_error_row(io, q1, q2, m1, m2, c, r0, L, E, label)
                            println("  ERR asym: $label m=($m1,$m2) c=$c L=$L E=$E : $ex")
                        end
                        run_count += 1
                        run_count % 50 == 0 && println("  [$run_count runs, $(round(time()-t_start,digits=1))s]")
                    end
                end
            end
        end
    end

    total_elapsed = time() - t_start
    close(io)
    println("\n=== Census complete: $run_count total runs in $(round(total_elapsed, digits=1))s ===")
    println("Results written to $outfile")
    return run_count
end

"""
Build q0, p0 for like-charge pair at turning point r0 with angular momentum L.
"""
function q0_p0(q1, q2, m1, m2, r0, L, mu)
    M = m1 + m2
    q0 = [-(m2 / M) * r0, 0.0, (m1 / M) * r0, 0.0]
    if L == 0.0
        p0 = [0.0, 0.0, 0.0, 0.0]
    else
        p_perp = L / r0
        p0 = [0.0, p_perp, 0.0, -p_perp]
    end
    return (q0, p0)
end

# Run the census
basedir = @__DIR__
run_census(joinpath(basedir, "census_results.csv"))
