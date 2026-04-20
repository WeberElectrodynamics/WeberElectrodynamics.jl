#!/usr/bin/env julia
# Agent 14 — Bound orbit verification pipeline for Weber electrodynamics.
#
# Provides:
#   verify_orbit()       — full 5-step verification pipeline for a single IC
#   detect_period()      — autocorrelation-based period detection
#   classify_orbit()     — periodic / quasi-periodic / chaotic-bound / escape
#
# Usage:
#     julia --project=. research/homology/14_verification/verify_bound.jl
#
# Writes: verification_results.csv alongside this script.

using LinearAlgebra
using Printf
using Statistics
using DelimitedFiles

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator

const HERE = @__DIR__
const REPO = normpath(joinpath(HERE, "..", "..", ".."))

# ============================================================================
# Utility: unflatten q/p vectors into n×d matrices
# ============================================================================
function unflatten(v::AbstractVector, n::Int, d::Int)
    return reshape(copy(v), d, n)'  # n×d
end

function pair_distances(q::AbstractVector, n::Int, d::Int)
    X = unflatten(q, n, d)
    out = Float64[]
    for i = 1:n, j = (i+1):n
        push!(out, norm(view(X, i, :) .- view(X, j, :)))
    end
    return out
end

# ============================================================================
# Step 0: Integration helper
# ============================================================================
function integrate(sys, q0, p0, masses, charges, c, tmax, dt;
                   bounce_r=0.0, dims=nothing)
    prob = HamiltonianProblem(
        sys, (0.0, tmax), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt,
    )
    return bounce_r > 0 ?
        solve(prob, SymmetricProjectionIntegrator(); callbacks=CollisionBounce(bounce_r)) :
        solve(prob, SymmetricProjectionIntegrator())
end

# ============================================================================
# Step 1: Period detection via autocorrelation of pair-distance signal
# ============================================================================
"""
    detect_period(sol, n, d; min_periods=1.5, skip_fraction=0.05) -> (T, quality)

Detect the dominant period from autocorrelation of the pair-distance signal.
Returns (NaN, 0.0) if no period is found.
`quality` ∈ [0,1] is the normalized autocorrelation peak height.
"""
function detect_period(sol, n::Int, d::Int;
                       min_periods::Float64=1.5, skip_fraction::Float64=0.05)
    npts = length(sol.t)
    npts < 100 && return (NaN, 0.0)

    # Build pair-distance signal (use first pair for 2-body, all pairs for n-body)
    r_signal = Vector{Float64}(undef, npts)
    for k = 1:npts
        dists = pair_distances(sol.q[k], n, d)
        r_signal[k] = dists[1]  # primary pair distance
    end

    dt_avg = (sol.t[end] - sol.t[1]) / (npts - 1)

    # Subtract mean for autocorrelation
    r_mean = mean(r_signal)
    r_centered = r_signal .- r_mean
    var_r = dot(r_centered, r_centered)
    var_r < 1e-30 && return (NaN, 0.0)  # constant signal (circular orbit)

    # For circular orbits, check phase-space return directly
    z0 = vcat(sol.q[1], sol.p[1])
    z_norm = norm(z0)
    if z_norm > 0 && var_r / (r_mean^2 * npts) < 1e-10
        # Nearly constant r — likely circular. Period from angular velocity.
        # Check full phase-space return at small intervals
        skip = max(10, round(Int, skip_fraction * npts))
        # Find FIRST local minimum of return distance that's good enough
        prev_ret = Inf
        for i = skip:npts
            zi = vcat(sol.q[i], sol.p[i])
            ret = norm(zi .- z0) / max(z_norm, 1e-15)
            # Detect first local minimum: distance was decreasing and starts increasing
            if ret > prev_ret && prev_ret < 0.01
                return (sol.t[i-1] - sol.t[1], 1.0 - prev_ret)
            end
            prev_ret = ret
        end
        # Check if the last point is the best
        if prev_ret < 0.01
            return (sol.t[end] - sol.t[1], 1.0 - prev_ret)
        end
    end

    # Compute autocorrelation for lags up to npts/2
    max_lag = npts ÷ 2
    skip_lag = max(1, round(Int, skip_fraction * npts))
    acf = Vector{Float64}(undef, max_lag)
    for lag = 1:max_lag
        s = 0.0
        @inbounds for i = 1:(npts - lag)
            s += r_centered[i] * r_centered[i + lag]
        end
        acf[lag] = s / var_r
    end

    # Find first significant peak after the skip window
    best_peak = -Inf
    best_lag = 0
    for lag = (skip_lag+1):(max_lag-1)
        if acf[lag] > acf[lag-1] && acf[lag] > acf[lag+1] && acf[lag] > 0.3
            if acf[lag] > best_peak
                best_peak = acf[lag]
                best_lag = lag
                break  # take the first (shortest period) significant peak
            end
        end
    end

    best_lag == 0 && return (NaN, 0.0)

    T_candidate = best_lag * dt_avg
    tmax = sol.t[end] - sol.t[1]
    # Require at least min_periods periods visible
    if T_candidate * min_periods > tmax
        return (NaN, 0.0)
    end

    return (T_candidate, best_peak)
end

# ============================================================================
# Phase-space return distance at time T
# ============================================================================
function phase_return_distance(sol, T::Float64)
    # Find index closest to T
    i_T = 1
    best_dt = Inf
    for i = 1:length(sol.t)
        d = abs(sol.t[i] - (sol.t[1] + T))
        if d < best_dt
            best_dt = d
            i_T = i
        end
    end
    z0 = vcat(sol.q[1], sol.p[1])
    zT = vcat(sol.q[i_T], sol.p[i_T])
    return norm(zT .- z0)
end

# ============================================================================
# Energy computation (Weber Hamiltonian)
# ============================================================================
function weber_energy(q, p, masses, charges, c, n, d)
    X = unflatten(q, n, d)
    P = unflatten(p, n, d)
    V = P ./ masses
    # Kinetic energy
    T_kin = 0.0
    for i = 1:n
        T_kin += 0.5 * masses[i] * dot(view(V, i, :), view(V, i, :))
    end
    # Weber potential
    U = 0.0
    for i = 1:n, j = (i+1):n
        dx = view(X, i, :) .- view(X, j, :)
        dv = view(V, i, :) .- view(V, j, :)
        r = norm(dx)
        rdot = dot(dx, dv) / r
        U += charges[i] * charges[j] / r * (1.0 - rdot^2 / (2.0 * c^2))
    end
    return T_kin + U
end

# ============================================================================
# Step 2-5: Full verification pipeline
# ============================================================================
"""
    VerificationResult

Complete classification of a single trajectory.

Fields:
- `id`: identifier string
- `classification`: :Periodic, :QuasiPeriodic, :ChaoticBound, :Escape, :Failure
- `period`: detected period T (NaN if none)
- `period_quality`: autocorrelation peak height
- `phase_return`: ||z(T) - z(0)|| at the detected period
- `energy_drift_base`: max |ΔH/H₀| at base resolution
- `energy_drift_extended`: max |ΔH/H₀| over 5T
- `convergence_discrepancy`: max ||z_base - z_half|| at T
- `max_separation_ratio`: max r_ij / min r_ij over integration
- `lyapunov_separation`: ||z(2T) - z_pert(2T)|| for perturbation test
- `n_particles`, `dims`, `c`: problem parameters
- `notes`: free-form string
"""
struct VerificationResult
    id::String
    classification::Symbol
    period::Float64
    period_quality::Float64
    phase_return::Float64
    energy_drift_base::Float64
    energy_drift_extended::Float64
    convergence_discrepancy::Float64
    max_separation_ratio::Float64
    lyapunov_separation::Float64
    n_particles::Int
    dims::Int
    c::Float64
    notes::String
end

"""
    verify_orbit(id, q0, p0, masses, charges, c;
                 n_particles, dims, tmax=100.0, dt=1e-3,
                 n_periods_check=5, energy_tol=1e-4,
                 periodic_tol=1e-6, bounce_r=0.0,
                 escape_ratio=10.0, perturbation_eps=1e-4,
                 verbose=true) -> VerificationResult

Run the full 5-step verification pipeline on a single orbit.

Steps:
1. Base integration + period detection
2. Convergence check (dt/2 re-integration)
3. Extended integration (5T)
4. Perturbation test (Lyapunov-type)
5. Classification
"""
function verify_orbit(id::String, q0::Vector{Float64}, p0::Vector{Float64},
                      masses::Vector{Float64}, charges::Vector{Float64}, c::Float64;
                      n_particles::Int, dims::Int,
                      tmax::Float64=100.0, dt::Float64=1e-3,
                      n_periods_check::Int=5,
                      energy_tol::Float64=1e-4,
                      periodic_tol::Float64=1e-6,
                      bounce_r::Float64=0.0,
                      escape_ratio::Float64=10.0,
                      perturbation_eps::Float64=1e-4,
                      verbose::Bool=true)

    sys = HamiltonianSystem(n_particles, dims)
    n = n_particles
    d = dims

    # ------------------------------------------------------------------
    # Step 1: Base integration
    # ------------------------------------------------------------------
    verbose && @printf("  [%s] Step 1: base integration (tmax=%.1f, dt=%.1e) ...\n", id, tmax, dt)
    sol = integrate(sys, q0, p0, masses, charges, c, tmax, dt; bounce_r=bounce_r)
    if sol.retcode != :Success
        verbose && @printf("  [%s] FAILED: retcode=%s\n", id, sol.retcode)
        return VerificationResult(id, :Failure, NaN, 0.0, NaN, NaN, NaN, NaN, NaN, NaN,
                                   n, d, c, "retcode=$(sol.retcode)")
    end

    # Energy drift at base resolution
    E0 = weber_energy(sol.q[1], sol.p[1], masses, charges, c, n, d)
    E_max_drift = 0.0
    for k = 1:length(sol.t)
        Ek = weber_energy(sol.q[k], sol.p[k], masses, charges, c, n, d)
        drift = abs(Ek - E0) / max(abs(E0), 1e-15)
        E_max_drift = max(E_max_drift, drift)
    end
    verbose && @printf("  [%s]   energy drift (base): %.2e\n", id, E_max_drift)

    # Check escape
    r0_max = maximum(pair_distances(q0, n, d))
    r_max_seen = 0.0
    r_min_seen = Inf
    for k = 1:length(sol.t)
        dists = pair_distances(sol.q[k], n, d)
        r_max_seen = max(r_max_seen, maximum(dists))
        r_min_seen = min(r_min_seen, minimum(dists))
    end
    separation_ratio = r_max_seen / max(r_min_seen, 1e-15)

    if r_max_seen > escape_ratio * max(r0_max, 1.0)
        verbose && @printf("  [%s] ESCAPE: r_max=%.2f > %.2f × r0_max=%.2f\n",
                           id, r_max_seen, escape_ratio, r0_max)
        return VerificationResult(id, :Escape, NaN, 0.0, NaN, E_max_drift, NaN, NaN,
                                   separation_ratio, NaN, n, d, c,
                                   "max_r=$(round(r_max_seen, digits=3))")
    end

    # Period detection
    T_detected, quality = detect_period(sol, n, d)
    verbose && @printf("  [%s]   period: T=%.4f  quality=%.3f\n", id,
                       isnan(T_detected) ? -1.0 : T_detected, quality)

    # If no period detected from autocorrelation, try direct phase-space scan
    if isnan(T_detected)
        # Scan for first local minimum of return distance
        z0 = vcat(q0, p0)
        z_norm = max(norm(z0), 1e-15)
        skip = max(50, length(sol.t) ÷ 20)
        prev_ret = Inf
        for i = skip:length(sol.t)
            zi = vcat(sol.q[i], sol.p[i])
            ret = norm(zi .- z0) / z_norm
            if ret > prev_ret && prev_ret < 0.05
                T_detected = sol.t[i-1] - sol.t[1]
                quality = 1.0 - prev_ret
                verbose && @printf("  [%s]   phase-scan period: T=%.4f  return=%.4e\n",
                                   id, T_detected, prev_ret)
                break
            end
            prev_ret = ret
        end
    end

    phase_ret = isnan(T_detected) ? NaN : phase_return_distance(sol, T_detected)
    if !isnan(phase_ret) && verbose
        @printf("  [%s]   phase return ||z(T)-z(0)||=%.4e\n", id, phase_ret)
    end

    # ------------------------------------------------------------------
    # Step 2: Convergence check (dt/2)
    # ------------------------------------------------------------------
    verbose && @printf("  [%s] Step 2: convergence check (dt/2) ...\n", id)
    tmax_conv = isnan(T_detected) ? tmax : min(2.0 * T_detected, tmax)
    sol_half = integrate(sys, q0, p0, masses, charges, c, tmax_conv, dt/2; bounce_r=bounce_r)
    convergence_disc = NaN
    if sol_half.retcode == :Success
        # Compare at the end of the shorter integration
        t_compare = min(sol.t[end], sol_half.t[end])
        # Find matching indices
        i_base = 1
        for i = 1:length(sol.t)
            if abs(sol.t[i] - t_compare) < abs(sol.t[i_base] - t_compare)
                i_base = i
            end
        end
        i_half = 1
        for i = 1:length(sol_half.t)
            if abs(sol_half.t[i] - t_compare) < abs(sol_half.t[i_half] - t_compare)
                i_half = i
            end
        end
        z_base = vcat(sol.q[i_base], sol.p[i_base])
        z_half = vcat(sol_half.q[i_half], sol_half.p[i_half])
        convergence_disc = norm(z_base .- z_half)
        verbose && @printf("  [%s]   convergence discrepancy: %.4e\n", id, convergence_disc)
    else
        verbose && @printf("  [%s]   convergence check: half-dt solve failed\n", id)
    end

    # ------------------------------------------------------------------
    # Step 3: Extended integration (n_periods_check × T)
    # ------------------------------------------------------------------
    E_max_drift_ext = E_max_drift
    if !isnan(T_detected)
        tmax_ext = n_periods_check * T_detected
        if tmax_ext > tmax
            verbose && @printf("  [%s] Step 3: extended integration (%.1f = %d×T) ...\n",
                               id, tmax_ext, n_periods_check)
            sol_ext = integrate(sys, q0, p0, masses, charges, c, tmax_ext, dt; bounce_r=bounce_r)
            if sol_ext.retcode == :Success
                for k = 1:length(sol_ext.t)
                    Ek = weber_energy(sol_ext.q[k], sol_ext.p[k], masses, charges, c, n, d)
                    drift = abs(Ek - E0) / max(abs(E0), 1e-15)
                    E_max_drift_ext = max(E_max_drift_ext, drift)
                end
                # Update escape check
                for k = 1:length(sol_ext.t)
                    dists = pair_distances(sol_ext.q[k], n, d)
                    if maximum(dists) > escape_ratio * max(r0_max, 1.0)
                        verbose && @printf("  [%s] ESCAPE in extended integration\n", id)
                        return VerificationResult(id, :Escape, T_detected, quality,
                                                   phase_ret, E_max_drift, E_max_drift_ext,
                                                   convergence_disc, separation_ratio, NaN,
                                                   n, d, c, "escape in extended integration")
                    end
                end
                verbose && @printf("  [%s]   extended energy drift: %.2e\n", id, E_max_drift_ext)
            else
                verbose && @printf("  [%s]   extended integration: FAILED\n", id)
            end
        else
            verbose && @printf("  [%s] Step 3: skipped (tmax >= %d×T)\n", id, n_periods_check)
        end
    else
        verbose && @printf("  [%s] Step 3: skipped (no period detected)\n", id)
    end

    # ------------------------------------------------------------------
    # Step 4: Perturbation test (Lyapunov-type separation)
    # ------------------------------------------------------------------
    verbose && @printf("  [%s] Step 4: perturbation test (eps=%.1e) ...\n", id, perturbation_eps)
    lyapunov_sep = NaN
    tmax_pert = isnan(T_detected) ? tmax : min(2.0 * T_detected, tmax)
    # Create perturbation in a random but reproducible direction
    z0 = vcat(q0, p0)
    # Use a deterministic perturbation direction based on id hash
    rng_seed = abs(hash(id)) % 10000
    delta = zeros(length(z0))
    for i = 1:length(delta)
        # Simple deterministic pseudo-random
        delta[i] = sin(rng_seed * 0.618033988749 * i)
    end
    delta .*= perturbation_eps / max(norm(delta), 1e-15)

    q0_pert = q0 .+ delta[1:n*d]
    p0_pert = p0 .+ delta[n*d+1:end]

    sol_pert = integrate(sys, q0_pert, p0_pert, masses, charges, c, tmax_pert, dt; bounce_r=bounce_r)
    if sol_pert.retcode == :Success
        # Compare at end
        z_ref = vcat(sol.q[end], sol.p[end])
        t_end_ref = sol.t[end]
        # Find matching time in perturbed solution
        i_match = 1
        for i = 1:length(sol_pert.t)
            if abs(sol_pert.t[i] - t_end_ref) < abs(sol_pert.t[i_match] - t_end_ref)
                i_match = i
            end
        end
        z_pert_end = vcat(sol_pert.q[i_match], sol_pert.p[i_match])

        # Also find separation at tmax_pert
        z_end = vcat(sol.q[end], sol.p[end])
        # If tmax_pert > sol.t[end], use the extended solution
        if tmax_pert <= sol.t[end] + dt
            lyapunov_sep = norm(z_pert_end .- z_ref)
        else
            # Re-integrate reference to tmax_pert
            sol_ref2 = integrate(sys, q0, p0, masses, charges, c, tmax_pert, dt; bounce_r=bounce_r)
            if sol_ref2.retcode == :Success
                i_end = length(sol_ref2.t)
                i_pert_end = 1
                for i = 1:length(sol_pert.t)
                    if abs(sol_pert.t[i] - sol_ref2.t[i_end]) < abs(sol_pert.t[i_pert_end] - sol_ref2.t[i_end])
                        i_pert_end = i
                    end
                end
                z_ref_end = vcat(sol_ref2.q[i_end], sol_ref2.p[i_end])
                z_pert_final = vcat(sol_pert.q[i_pert_end], sol_pert.p[i_pert_end])
                lyapunov_sep = norm(z_pert_final .- z_ref_end)
            end
        end
        if !isnan(lyapunov_sep)
            verbose && @printf("  [%s]   Lyapunov separation: %.4e  (amplification: %.1f×)\n",
                               id, lyapunov_sep, lyapunov_sep / perturbation_eps)
        end
    else
        verbose && @printf("  [%s]   perturbation test: perturbed solve failed\n", id)
    end

    # ------------------------------------------------------------------
    # Step 5: Classification
    # ------------------------------------------------------------------
    classification = classify(phase_ret, E_max_drift, E_max_drift_ext,
                              lyapunov_sep, perturbation_eps, T_detected,
                              periodic_tol, energy_tol, r_max_seen, escape_ratio,
                              r0_max, z0)

    verbose && @printf("  [%s] CLASSIFICATION: %s\n\n", id, classification)

    return VerificationResult(id, classification, T_detected, quality, phase_ret,
                               E_max_drift, E_max_drift_ext, convergence_disc,
                               separation_ratio, lyapunov_sep, n, d, c, "")
end

function classify(phase_ret, E_drift_base, E_drift_ext, lyapunov_sep, pert_eps,
                  T, periodic_tol, energy_tol, r_max, escape_ratio, r0_max, z0)
    z_norm = max(norm(z0), 1e-15)

    # Escape
    if r_max > escape_ratio * max(r0_max, 1.0)
        return :Escape
    end

    # Periodic: tight phase-space return AND good energy conservation
    if !isnan(phase_ret) && !isnan(T)
        relative_return = phase_ret / z_norm
        if relative_return < periodic_tol && E_drift_ext < 1e-8
            return :Periodic
        end
        # Relaxed periodic: still close return with good energy
        if relative_return < 0.01 && E_drift_base < energy_tol
            return :Periodic
        end
    end

    # Chaotic but bound: large Lyapunov amplification (exponential growth)
    if !isnan(lyapunov_sep) && pert_eps > 0
        amplification = lyapunov_sep / pert_eps
        if amplification > 1e4
            return :ChaoticBound
        end
    end

    # Quasi-periodic: bounded, non-periodic return, stable under perturbation
    if !isnan(phase_ret) && !isnan(T)
        relative_return = phase_ret / z_norm
        if relative_return < 0.5 && E_drift_base < energy_tol
            return :QuasiPeriodic
        end
    end

    # Bounded but no period detected — still classify based on confinement
    if E_drift_base < energy_tol && r_max < escape_ratio * max(r0_max, 1.0)
        if isnan(T)
            return :QuasiPeriodic  # bounded but no clear period
        end
    end

    return :Unclassified
end

# ============================================================================
# Test case definitions
# ============================================================================

function test_case_2body_circular()
    # 2-body circular orbit: q₁=+1, q₂=-1, c=4, r₀=2
    # Circular orbit: v = sqrt(|q1*q2| / (μ * r0)) with μ = m1*m2/(m1+m2)
    n = 2; d = 2; c = 4.0
    m = [1.0, 1.0]; q_charges = [1.0, -1.0]
    r0 = 2.0
    mu = 0.5
    v_circ = sqrt(abs(q_charges[1] * q_charges[2]) / (mu * r0))
    # Place particles at ±r0/2 along x, circular velocity along y
    q0 = [r0/2, 0.0, -r0/2, 0.0]
    p0 = [0.0, m[1]*v_circ/2, 0.0, -m[2]*v_circ/2]
    # Analytical period: T = 2πr0 / v_circ (COM frame, each at r0/2)
    # Actually in COM: each particle at r0/2, speed v_circ/2, so T = 2π(r0/2)/(v_circ/2) = 2πr0/v_circ
    T_analytical = 2π * r0 / v_circ
    return (id="2body_circular", q0=q0, p0=p0, masses=m, charges=q_charges,
            c=c, n=n, d=d, tmax=5*T_analytical, dt=1e-3, bounce_r=0.0, escape_ratio=10.0,
            expected=:Periodic, T_expected=T_analytical)
end

function test_case_2body_elliptical()
    # 2-body elliptical orbit: q₁=+1, q₂=-1, c=4, eccentricity ≈ 0.3
    n = 2; d = 2; c = 4.0
    m = [1.0, 1.0]; q_charges = [1.0, -1.0]
    mu = 0.5
    k = abs(q_charges[1] * q_charges[2])  # = 1.0
    # For e=0.3, start at apoapsis: r_apo = a(1+e), v_apo from vis-viva
    a = 2.0  # semi-major axis
    e = 0.3
    r_apo = a * (1 + e)
    # Energy: E = -k/(2a), angular momentum: L = μ*sqrt(k*a*(1-e²)/μ)
    # At apoapsis: v = L/(μ*r_apo)
    L = mu * sqrt(k * a * (1 - e^2) / mu)
    v_apo = L / (mu * r_apo)
    q0 = [r_apo/2, 0.0, -r_apo/2, 0.0]
    p0 = [0.0, m[1]*v_apo/2, 0.0, -m[2]*v_apo/2]
    T_kepler = 2π * sqrt(mu * a^3 / k)
    # Note: Weber correction causes orbital precession, so the orbit is
    # quasi-periodic (pair distance repeats but full phase space does not close
    # at the Kepler period). The radial period is close to but not equal to T_kepler.
    return (id="2body_elliptical_e03", q0=q0, p0=p0, masses=m, charges=q_charges,
            c=c, n=n, d=d, tmax=8*T_kepler, dt=1e-3, bounce_r=0.0, escape_ratio=10.0,
            expected=:QuasiPeriodic, T_expected=T_kepler)
end

function test_case_2body_circular_large_c()
    # 2-body circular orbit with large c (near-Coulomb limit)
    # This tests that near-Coulomb orbits are cleanly periodic.
    n = 2; d = 2; c = 100.0
    m = [1.0, 1.0]; q_charges = [1.0, -1.0]
    r0 = 2.0
    mu = 0.5
    v_circ = sqrt(abs(q_charges[1] * q_charges[2]) / (mu * r0))
    q0 = [r0/2, 0.0, -r0/2, 0.0]
    p0 = [0.0, m[1]*v_circ/2, 0.0, -m[2]*v_circ/2]
    T_analytical = 2π * r0 / v_circ
    return (id="2body_circular_large_c", q0=q0, p0=p0, masses=m, charges=q_charges,
            c=c, n=n, d=d, tmax=5*T_analytical, dt=1e-3, bounce_r=0.0, escape_ratio=10.0,
            expected=:Periodic, T_expected=T_analytical)
end

function test_case_4body_rhombus_chaotic()
    # 4-body rhombus: a=1.5, b=1.45, η=0.75, rotating, c=1
    n = 4; d = 2; c = 1.0
    a = 1.5; b = 1.45
    X = [a 0.0; -a 0.0; 0.0 b; 0.0 -b]
    masses = fill(1.0, 4)
    charges = [1.0, 1.0, -1.0, -1.0]
    # Energy fraction
    eta = 0.5
    U = 0.0
    for i = 1:4, j = (i+1):4
        r = norm(X[i, :] .- X[j, :])
        U += charges[i] * charges[j] / r
    end
    T_target = eta * abs(U)
    speed = sqrt(2 * T_target / (4 * 1.0))
    P = zeros(4, 2)
    for i = 1:4
        r = X[i, :]
        t = [-r[2], r[1]] / max(norm(r), eps())
        P[i, :] = 1.0 * speed * t
    end
    # Zero COM and total P
    com = vec(sum(masses .* X, dims=1) / sum(masses))
    X .-= com'
    P .-= sum(P, dims=1) ./ 4
    q0 = Float64[]
    for i = 1:4, k = 1:2
        push!(q0, X[i, k])
    end
    p0 = Float64[]
    for i = 1:4, k = 1:2
        push!(p0, P[i, k])
    end
    return (id="4body_rhombus_chaotic", q0=q0, p0=p0, masses=masses, charges=charges,
            c=c, n=n, d=d, tmax=200.0, dt=1e-4, bounce_r=0.02, escape_ratio=50.0,
            expected=:ChaoticBound, T_expected=NaN)
end

function test_case_2body_escape()
    # Clearly unbound: same-sign charges, positive energy
    n = 2; d = 2; c = 4.0
    m = [1.0, 1.0]; q_charges = [1.0, 1.0]
    q0 = [1.0, 0.0, -1.0, 0.0]
    p0 = [0.5, 0.0, -0.5, 0.0]  # moving apart
    return (id="2body_escape", q0=q0, p0=p0, masses=m, charges=q_charges,
            c=c, n=n, d=d, tmax=50.0, dt=1e-3, bounce_r=0.0, escape_ratio=10.0,
            expected=:Escape, T_expected=NaN)
end

# ============================================================================
# Main: run all test cases
# ============================================================================
function run_verification()
    test_cases = [
        test_case_2body_circular(),
        test_case_2body_elliptical(),
        test_case_2body_circular_large_c(),
        test_case_4body_rhombus_chaotic(),
        test_case_2body_escape(),
    ]

    results = VerificationResult[]

    println("=" ^ 80)
    println("Agent 14 — Bound Orbit Verification Pipeline")
    println("=" ^ 80)

    for tc in test_cases
        println("\n--- Test case: $(tc.id) (expected: $(tc.expected)) ---")
        result = verify_orbit(
            tc.id, Float64.(tc.q0), Float64.(tc.p0),
            Float64.(tc.masses), Float64.(tc.charges), tc.c;
            n_particles=tc.n, dims=tc.d,
            tmax=tc.tmax, dt=tc.dt,
            bounce_r=tc.bounce_r,
            escape_ratio=tc.escape_ratio,
            verbose=true,
        )
        push!(results, result)

        # Check against expected
        if result.classification == tc.expected
            @printf("  >>> PASS: classified as %s (expected %s)\n", result.classification, tc.expected)
        else
            @printf("  >>> MISMATCH: classified as %s (expected %s)\n", result.classification, tc.expected)
        end
        if !isnan(tc.T_expected) && !isnan(result.period)
            T_err = abs(result.period - tc.T_expected) / tc.T_expected
            @printf("  >>> Period: detected=%.4f  expected=%.4f  error=%.2e\n",
                    result.period, tc.T_expected, T_err)
        end
    end

    # Write CSV (compact format with requested columns)
    csv_path = joinpath(HERE, "verification_results.csv")
    open(csv_path, "w") do io
        println(io, "orbit_id,classification,period,energy_drift,lyapunov,n_periods_verified")
        for (r, tc) in zip(results, test_cases)
            n_periods = (!isnan(r.period) && r.period > 0) ? round(Int, tc.tmax / r.period) : 0
            cls_str = uppercase(replace(string(r.classification), r"([a-z])([A-Z])" => s"\1-\2"))
            @printf(io, "%s,%s,%s,%s,%s,%d\n",
                    r.id, cls_str,
                    isnan(r.period) ? "" : @sprintf("%.4f", r.period),
                    isnan(r.energy_drift_base) ? "" : @sprintf("%.2e", r.energy_drift_base),
                    isnan(r.lyapunov_separation) ? "" : @sprintf("%.2e", r.lyapunov_separation),
                    n_periods)
        end
    end
    @printf("\nWrote %d results to %s\n", length(results), csv_path)

    # Write detailed CSV with all fields
    csv_detail_path = joinpath(HERE, "verification_results_detailed.csv")
    open(csv_detail_path, "w") do io
        println(io, "id,classification,period,period_quality,phase_return,energy_drift_base,energy_drift_extended,convergence_discrepancy,max_separation_ratio,lyapunov_separation,n_particles,dims,c,notes")
        for r in results
            @printf(io, "%s,%s,%.8f,%.6f,%.6e,%.6e,%.6e,%.6e,%.4f,%.6e,%d,%d,%.2f,%s\n",
                    r.id, r.classification,
                    isnan(r.period) ? -1.0 : r.period,
                    r.period_quality,
                    isnan(r.phase_return) ? -1.0 : r.phase_return,
                    isnan(r.energy_drift_base) ? -1.0 : r.energy_drift_base,
                    isnan(r.energy_drift_extended) ? -1.0 : r.energy_drift_extended,
                    isnan(r.convergence_discrepancy) ? -1.0 : r.convergence_discrepancy,
                    isnan(r.max_separation_ratio) ? -1.0 : r.max_separation_ratio,
                    isnan(r.lyapunov_separation) ? -1.0 : r.lyapunov_separation,
                    r.n_particles, r.dims, r.c, r.notes)
        end
    end

    # Summary table
    println("\n" * "=" ^ 80)
    println("SUMMARY")
    println("=" ^ 80)
    @printf("%-30s  %-15s  %-10s  %-12s  %-12s\n",
            "ID", "Classification", "Period", "E_drift", "Lyap_sep")
    println("-" ^ 80)
    for r in results
        @printf("%-30s  %-15s  %-10s  %-12s  %-12s\n",
                r.id, r.classification,
                isnan(r.period) ? "   -" : @sprintf("%.4f", r.period),
                isnan(r.energy_drift_base) ? "   -" : @sprintf("%.2e", r.energy_drift_base),
                isnan(r.lyapunov_separation) ? "   -" : @sprintf("%.2e", r.lyapunov_separation))
    end

    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_verification()
end
