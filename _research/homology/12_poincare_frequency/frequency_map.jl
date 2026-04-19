"""
Agent 12 -- Laskar frequency map analysis for the 4-body 2+/2- Weber Hamiltonian.

Extracts fundamental frequencies from pair-distance time series using windowed
DFT (no FFTW dependency -- uses manual DFT for peak finding).

Frequency map: plot (omega1, omega2) for a grid of ICs.
Frequency drift: track omega(t) over successive windows to measure Arnold diffusion.
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
# Pair distance time series extraction
# ============================================================================
function pair_dist_series(sol, i, j)
    n = length(sol.t)
    r = Vector{Float64}(undef, n)
    for k in 1:n
        X = reshape(sol.q[k], D, N)
        r[k] = norm(X[:, i] .- X[:, j])
    end
    return r
end

# ============================================================================
# Windowed DFT: find dominant frequency in a signal segment
# Uses zero-padded DFT with Hanning window for accuracy.
# ============================================================================
"""
Dominant frequency via manual DFT with Hanning window.
Signals are downsampled to at most `max_pts` before DFT to keep the O(n^2)
manual DFT tractable (no FFTW dependency).
"""
function dominant_frequency(signal::AbstractVector{Float64}, dt::Float64;
                             n_pad_factor::Int=4, max_pts::Int=4000)
    n = length(signal)
    n < 4 && return NaN

    # Downsample if necessary to keep manual DFT tractable
    stride = max(1, n ÷ max_pts)
    if stride > 1
        sig_ds = signal[1:stride:end]
        dt_ds = dt * stride
    else
        sig_ds = signal
        dt_ds = dt
    end
    n_ds = length(sig_ds)

    # Remove mean
    mu = mean(sig_ds)
    s = sig_ds .- mu

    # Hanning window
    for i in 1:n_ds
        s[i] *= 0.5 * (1 - cos(2pi * (i-1) / (n_ds-1)))
    end

    # Zero-pad for frequency resolution
    n_fft = n_ds * n_pad_factor
    padded = zeros(n_fft)
    padded[1:n_ds] .= s

    # Manual DFT (real input, only need positive frequencies)
    n_freq = n_fft ÷ 2
    power = zeros(n_freq)
    freq = zeros(n_freq)
    for k in 1:n_freq
        freq[k] = (k - 1) / (n_fft * dt_ds)
        re = 0.0
        im = 0.0
        phase_step = -2pi * (k-1) / n_fft
        for j in 1:n_ds  # only non-zero entries
            phase = phase_step * (j - 1)
            re += padded[j] * cos(phase)
            im += padded[j] * sin(phase)
        end
        power[k] = re^2 + im^2
    end

    # Skip DC (k=1) and very low frequencies
    min_k = max(2, round(Int, 0.01 / (1.0/(n_fft*dt_ds))) + 1)
    max_k = n_freq

    if min_k > max_k
        return NaN
    end

    # Find peak
    peak_k = min_k
    peak_pow = power[min_k]
    for k in (min_k+1):max_k
        if power[k] > peak_pow
            peak_pow = power[k]
            peak_k = k
        end
    end

    # Parabolic interpolation around peak for sub-bin accuracy
    if peak_k > 1 && peak_k < n_freq
        alpha = log(max(power[peak_k-1], 1e-300))
        beta = log(max(power[peak_k], 1e-300))
        gamma = log(max(power[peak_k+1], 1e-300))
        delta = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma)
        if abs(delta) < 1.0
            return (peak_k - 1 + delta) / (n_fft * dt_ds)
        end
    end

    return freq[peak_k]
end

# ============================================================================
# Double-orbiter IC (same as poincare_sections.jl)
# ============================================================================
function double_orbiter_ic(; r_pp=4.0, R=4.0, orb=1.3, m=1.0, q=1.0)
    X = zeros(4, 2)
    X[1, :] = [-r_pp/2, 0.0]
    X[2, :] = [r_pp/2, 0.0]
    X[3, :] = [0.0, R]
    X[4, :] = [0.0, -R]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    P = zeros(4, 2)
    P[3, :] = m * orb * [1.0, 0.0]
    P[4, :] = m * orb * [-1.0, 0.0]
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# ============================================================================
# Frequency map for a single IC
# ============================================================================
function compute_frequencies(sol; window_frac=0.5)
    dt = sol.t[2] - sol.t[1]
    nsteps = length(sol.t)
    window_len = round(Int, window_frac * nsteps)
    window_len = min(window_len, nsteps)

    r12 = pair_dist_series(sol, 1, 2)
    r34 = pair_dist_series(sol, 3, 4)

    # Use the first window_frac of the trajectory
    seg12 = r12[1:window_len]
    seg34 = r34[1:window_len]

    omega1 = dominant_frequency(seg12, dt)
    omega2 = dominant_frequency(seg34, dt)
    return (omega1=omega1, omega2=omega2)
end

# ============================================================================
# Frequency drift: track omega over successive windows
# ============================================================================
function frequency_drift(sol; n_windows=10, pair1=(1,2), pair2=(3,4))
    dt = sol.t[2] - sol.t[1]
    nsteps = length(sol.t)
    window_len = nsteps ÷ n_windows

    if window_len < 20
        return (t_centers=Float64[], omega1=Float64[], omega2=Float64[])
    end

    r1 = pair_dist_series(sol, pair1[1], pair1[2])
    r2 = pair_dist_series(sol, pair2[1], pair2[2])

    t_centers = Float64[]
    omega1_series = Float64[]
    omega2_series = Float64[]

    for w in 1:n_windows
        i_start = (w-1) * window_len + 1
        i_end = min(w * window_len, nsteps)
        i_end - i_start < 10 && continue

        seg1 = r1[i_start:i_end]
        seg2 = r2[i_start:i_end]

        o1 = dominant_frequency(seg1, dt)
        o2 = dominant_frequency(seg2, dt)

        push!(t_centers, (sol.t[i_start] + sol.t[i_end]) / 2)
        push!(omega1_series, o1)
        push!(omega2_series, o2)
    end

    return (t_centers=t_centers, omega1=omega1_series, omega2=omega2_series)
end

# ============================================================================
# Compute diffusion coefficient from frequency drift
# ============================================================================
function diffusion_coefficient(omega_series, t_centers)
    n = length(omega_series)
    n < 3 && return NaN
    # Filter NaN
    valid = findall(!isnan, omega_series)
    length(valid) < 3 && return NaN

    omegas = omega_series[valid]
    times = t_centers[valid]

    # Compute |Delta omega|^2 / Delta t over successive pairs
    diffs = Float64[]
    for i in 2:length(omegas)
        dt = times[i] - times[i-1]
        dt > 0 || continue
        push!(diffs, (omegas[i] - omegas[i-1])^2 / dt)
    end
    isempty(diffs) && return NaN
    return mean(diffs)
end

# ============================================================================
# Run frequency map for rhombus grid
# ============================================================================
function run_rhombus_frequency_map(; a_vals=[1.4, 1.45, 1.5, 1.55, 1.6],
                                     b_vals=[1.3, 1.4, 1.45, 1.5],
                                     eta_vals=[0.7, 0.75, 0.8],
                                     tmax=200.0, dt=5e-4)
    results = NamedTuple[]
    for a in a_vals, b in b_vals, eta in eta_vals
        label = @sprintf("rhombus_a%.2f_b%.2f_eta%.2f", a, b, eta)
        @info "Frequency: $label"
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
                           omega1=NaN, omega2=NaN, D1=NaN, D2=NaN))
            continue
        end

        if sol.t[end] < 20.0
            push!(results, (config=:rhombus, a=a, b=b, eta=eta,
                           retcode=sol.retcode, t_final=sol.t[end],
                           omega1=NaN, omega2=NaN, D1=NaN, D2=NaN))
            continue
        end

        freqs = compute_frequencies(sol)
        drift = frequency_drift(sol; n_windows=10)
        D1 = diffusion_coefficient(drift.omega1, drift.t_centers)
        D2 = diffusion_coefficient(drift.omega2, drift.t_centers)

        push!(results, (config=:rhombus, a=a, b=b, eta=eta,
                       retcode=sol.retcode, t_final=sol.t[end],
                       omega1=freqs.omega1, omega2=freqs.omega2,
                       D1=D1, D2=D2))

        @info "  => omega1=$(freqs.omega1), omega2=$(freqs.omega2), D1=$D1, D2=$D2"
    end
    return results
end

# ============================================================================
# Run frequency map for double-orbiter grid
# ============================================================================
function run_orbiter_frequency_map(; orb_vals=[1.0, 1.1, 1.2, 1.3, 1.4, 1.5],
                                     rpp_vals=[3.5, 4.0, 4.5],
                                     tmax=200.0, dt=5e-4)
    results = NamedTuple[]
    for orb in orb_vals, rpp in rpp_vals
        label = @sprintf("orbiter_rpp%.1f_orb%.1f", rpp, orb)
        @info "Frequency: $label"
        q0, p0, masses, charges = double_orbiter_ic(r_pp=rpp, R=4.0, orb=orb)
        sol = nothing
        try
            sol = SharedSurvey.run(q0, p0, masses, charges;
                                    tmax=tmax, dt=dt, c=1.0, dims=D,
                                    bounce_r=0.02)
        catch err
            @warn "Failed" label err
            push!(results, (config=:double_orbiter, r_pp=rpp, orb=orb,
                           retcode=:Exception, t_final=0.0,
                           omega1=NaN, omega2=NaN, D1=NaN, D2=NaN))
            continue
        end

        if sol.t[end] < 20.0
            push!(results, (config=:double_orbiter, r_pp=rpp, orb=orb,
                           retcode=sol.retcode, t_final=sol.t[end],
                           omega1=NaN, omega2=NaN, D1=NaN, D2=NaN))
            continue
        end

        freqs = compute_frequencies(sol)
        drift = frequency_drift(sol; n_windows=10)
        D1 = diffusion_coefficient(drift.omega1, drift.t_centers)
        D2 = diffusion_coefficient(drift.omega2, drift.t_centers)

        push!(results, (config=:double_orbiter, r_pp=rpp, orb=orb,
                       retcode=sol.retcode, t_final=sol.t[end],
                       omega1=freqs.omega1, omega2=freqs.omega2,
                       D1=D1, D2=D2))

        @info "  => omega1=$(freqs.omega1), omega2=$(freqs.omega2), D1=$D1, D2=$D2"
    end
    return results
end

# ============================================================================
# Detailed frequency drift for best candidates
# ============================================================================
function detailed_drift_analysis(; candidates=nothing, tmax=300.0, dt=5e-4, n_windows=15)
    if candidates === nothing
        # Default best candidates from prior studies
        candidates = [
            (type=:rhombus, params=(a=1.5, b=1.45, eta=0.75)),
            (type=:rhombus, params=(a=1.5, b=1.45, eta=0.70)),
            (type=:rhombus, params=(a=1.5, b=1.45, eta=0.80)),
            (type=:rhombus, params=(a=1.5, b=1.40, eta=0.75)),
            (type=:rhombus, params=(a=1.5, b=1.50, eta=0.75)),
        ]
    end

    results = NamedTuple[]
    for cand in candidates
        if cand.type == :rhombus
            p = cand.params
            label = @sprintf("drift_rhombus_a%.2f_b%.2f_eta%.2f", p.a, p.b, p.eta)
            @info "Detailed drift: $label"
            q0, p0, masses, charges = rhombus(a=p.a, b=p.b,
                                               energy_fraction=p.eta,
                                               velocity_mode=:rotating)
        else
            continue
        end

        sol = nothing
        try
            sol = SharedSurvey.run(q0, p0, masses, charges;
                                    tmax=tmax, dt=dt, c=1.0, dims=D,
                                    bounce_r=0.02)
        catch err
            @warn "Failed" label err
            continue
        end

        if sol.t[end] < 50.0
            @info "  Too short: t=$(sol.t[end])"
            continue
        end

        drift = frequency_drift(sol; n_windows=n_windows)
        D1 = diffusion_coefficient(drift.omega1, drift.t_centers)
        D2 = diffusion_coefficient(drift.omega2, drift.t_centers)

        # Compute omega variance as a function of time lag
        if length(drift.omega1) >= 3
            valid_o1 = filter(!isnan, drift.omega1)
            valid_o2 = filter(!isnan, drift.omega2)
            omega1_std = length(valid_o1) >= 2 ? std(valid_o1) : NaN
            omega2_std = length(valid_o2) >= 2 ? std(valid_o2) : NaN
            omega1_mean = length(valid_o1) >= 1 ? mean(valid_o1) : NaN
            omega2_mean = length(valid_o2) >= 1 ? mean(valid_o2) : NaN
        else
            omega1_std = omega2_std = omega1_mean = omega2_mean = NaN
        end

        push!(results, (label=label, t_final=sol.t[end],
                       omega1_mean=omega1_mean, omega2_mean=omega2_mean,
                       omega1_std=omega1_std, omega2_std=omega2_std,
                       D1=D1, D2=D2,
                       drift_data=drift))

        # Save drift data
        open(joinpath(HERE, "section_data", "$label.csv"), "w") do io
            println(io, "t_center,omega1,omega2")
            for i in 1:length(drift.t_centers)
                @printf(io, "%.8e,%.8e,%.8e\n",
                        drift.t_centers[i], drift.omega1[i], drift.omega2[i])
            end
        end

        @info "  => omega1=$(omega1_mean)+-$(omega1_std), D1=$D1"
        @info "  => omega2=$(omega2_mean)+-$(omega2_std), D2=$D2"
    end
    return results
end

# ============================================================================
# Write frequency results CSV
# ============================================================================
function write_frequency_results(rhombus_freq, orbiter_freq, drift_results)
    open(joinpath(HERE, "frequency_results.csv"), "w") do io
        println(io, "config,param1,param2,param3,retcode,t_final,omega1,omega2,D1,D2")
        for r in rhombus_freq
            @printf(io, "%s,%.3f,%.3f,%.3f,%s,%.4f,%.8e,%.8e,%.8e,%.8e\n",
                    r.config, r.a, r.b, r.eta, r.retcode, r.t_final,
                    r.omega1, r.omega2, r.D1, r.D2)
        end
        for r in orbiter_freq
            @printf(io, "%s,%.3f,%.3f,NaN,%s,%.4f,%.8e,%.8e,%.8e,%.8e\n",
                    r.config, r.r_pp, r.orb, r.retcode, r.t_final,
                    r.omega1, r.omega2, r.D1, r.D2)
        end
    end

    # Drift summary
    open(joinpath(HERE, "drift_summary.csv"), "w") do io
        println(io, "label,t_final,omega1_mean,omega2_mean,omega1_std,omega2_std,D1,D2")
        for r in drift_results
            @printf(io, "%s,%.4f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    r.label, r.t_final, r.omega1_mean, r.omega2_mean,
                    r.omega1_std, r.omega2_std, r.D1, r.D2)
        end
    end
end

# ============================================================================
# Main
# ============================================================================
function main()
    println("=" ^ 70)
    println("Agent 12 — Laskar frequency map analysis")
    println("=" ^ 70)

    println("\n--- Phase 1: Rhombus frequency map ---")
    rhombus_freq = run_rhombus_frequency_map()

    println("\n--- Phase 2: Double-orbiter frequency map ---")
    orbiter_freq = run_orbiter_frequency_map()

    println("\n--- Phase 3: Detailed drift analysis (top 5 candidates) ---")
    drift_results = detailed_drift_analysis()

    write_frequency_results(rhombus_freq, orbiter_freq, drift_results)

    # Print frequency map summary
    println("\n" * "=" ^ 70)
    println("FREQUENCY MAP RESULTS")
    println("=" ^ 70)

    println("\nRhombus ICs with valid frequencies:")
    valid_r = filter(r -> !isnan(r.omega1) && !isnan(r.omega2), rhombus_freq)
    for r in valid_r
        @printf("  a=%.2f b=%.2f eta=%.2f: omega1=%.4f omega2=%.4f D1=%.2e D2=%.2e\n",
                r.a, r.b, r.eta, r.omega1, r.omega2, r.D1, r.D2)
    end

    println("\nDouble-orbiter ICs with valid frequencies:")
    valid_o = filter(r -> !isnan(r.omega1) && !isnan(r.omega2), orbiter_freq)
    for r in valid_o
        @printf("  rpp=%.1f orb=%.1f: omega1=%.4f omega2=%.4f D1=%.2e D2=%.2e\n",
                r.r_pp, r.orb, r.omega1, r.omega2, r.D1, r.D2)
    end

    println("\n--- Drift analysis ---")
    for r in drift_results
        @printf("  %s: omega1=%.4f+-%.4f  omega2=%.4f+-%.4f  D1=%.2e D2=%.2e\n",
                r.label, r.omega1_mean, r.omega1_std,
                r.omega2_mean, r.omega2_std, r.D1, r.D2)
        # Relative drift
        if !isnan(r.omega1_std) && !isnan(r.omega1_mean) && r.omega1_mean != 0
            @printf("    omega1 relative drift: %.2f%%\n", 100*r.omega1_std/abs(r.omega1_mean))
        end
        if !isnan(r.omega2_std) && !isnan(r.omega2_mean) && r.omega2_mean != 0
            @printf("    omega2 relative drift: %.2f%%\n", 100*r.omega2_std/abs(r.omega2_mean))
        end
    end

    return (rhombus_freq=rhombus_freq, orbiter_freq=orbiter_freq, drift=drift_results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
