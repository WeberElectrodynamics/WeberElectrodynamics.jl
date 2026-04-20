# Agent 13 — c-continuation: track orbits from Coulomb (c→∞) to strong Weber (c=1)
#
# Run:  julia --project=. research/homology/13_c_continuation/continuation.jl
#
# Outputs:
#   bifurcation_diagram.csv — c vs period, energy, max Floquet, drift for each family
#   NOTES.md — auto-appended summary (written at end)

using LinearAlgebra
using Printf
using DelimitedFiles

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator

const OUTDIR = @__DIR__

# ============================================================================
# Core helpers
# ============================================================================

"""2-body circular orbit ICs in center-of-mass frame (2D)."""
function circular_ic_2d(r0; masses=[1.0,1.0], charges=[1.0,-1.0])
    mu = masses[1]*masses[2]/(masses[1]+masses[2])
    k = abs(charges[1]*charges[2])
    v = sqrt(k/(mu*r0))
    T = 2pi * sqrt(mu * r0^3 / k)
    # COM frame: particle 1 at -r0/2, particle 2 at +r0/2 along x
    q0 = [-r0/2, 0.0, r0/2, 0.0]
    v1 = -masses[2]/(masses[1]+masses[2]) * v
    v2 =  masses[1]/(masses[1]+masses[2]) * v
    p0 = [0.0, masses[1]*v1, 0.0, masses[2]*v2]
    E = -k/(2*r0)
    return q0, p0, T, E
end

"""2-body elliptical orbit ICs (2D), starting at periapsis along x-axis."""
function elliptical_ic_2d(E, ecc; masses=[1.0,1.0], charges=[1.0,-1.0])
    mu = masses[1]*masses[2]/(masses[1]+masses[2])
    k = abs(charges[1]*charges[2])
    a = k / (2*abs(E))  # semi-major axis
    r_peri = a * (1 - ecc)
    # At periapsis, all velocity is tangential
    # E = 0.5*mu*v^2 - k/r  =>  v = sqrt(2*(E + k/r)/mu)
    v_peri = sqrt(2*(E + k/r_peri)/mu)
    T_kepler = 2pi * sqrt(mu * a^3 / k)
    q0 = [-r_peri/2, 0.0, r_peri/2, 0.0]
    v1 = -masses[2]/(masses[1]+masses[2]) * v_peri
    v2 =  masses[1]/(masses[1]+masses[2]) * v_peri
    p0 = [0.0, masses[1]*v1, 0.0, masses[2]*v2]
    return q0, p0, T_kepler, E
end

"""Integrate a 2-body system for a given time span."""
function integrate_2body(sys, q0, p0, T, masses, charges, c; dt_hint=1e-3)
    dt = min(dt_hint, T/200)
    prob = HamiltonianProblem(sys, (0.0, T), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt)
    sol = solve(prob, SymmetricProjectionIntegrator())
    return sol
end

"""Compute energy drift as max fractional error."""
function energy_drift(sol)
    sol.retcode == :Success || return NaN
    en = compute_energy_timeseries(sol)
    return en.statistics.global_error_percent_max
end

"""Compute the monodromy matrix via centered finite differences."""
function monodromy_2body(sys, q0, p0, T, masses, charges, c; eps_fd=1e-7, dt_hint=1e-3)
    n = length(q0) + length(p0)
    z0 = vcat(q0, p0)
    ndof = length(q0)

    # Reference trajectory
    sol_ref = integrate_2body(sys, q0, p0, T, masses, charges, c; dt_hint=dt_hint)
    sol_ref.retcode == :Success || return nothing

    M = zeros(n, n)
    for k in 1:n
        zp = copy(z0); zp[k] += eps_fd
        zm = copy(z0); zm[k] -= eps_fd
        sp = integrate_2body(sys, zp[1:ndof], zp[ndof+1:end], T, masses, charges, c; dt_hint=dt_hint)
        sm = integrate_2body(sys, zm[1:ndof], zm[ndof+1:end], T, masses, charges, c; dt_hint=dt_hint)
        (sp.retcode != :Success || sm.retcode != :Success) && return nothing
        M[:, k] = (vcat(sp.q[end], sp.p[end]) - vcat(sm.q[end], sm.p[end])) / (2*eps_fd)
    end
    return M
end

"""Classify eigenvalue type for a 2x2 symplectic block."""
function classify_eigenvalue(lam)
    a = abs(lam)
    if abs(a - 1.0) < 1e-6
        return :elliptic
    elseif abs(imag(lam)) < 1e-8 && real(lam) > 0
        return :positive_hyperbolic
    elseif abs(imag(lam)) < 1e-8 && real(lam) < 0
        return :negative_hyperbolic
    else
        return :loxodromic
    end
end

# ============================================================================
# Part 1: 2-Body c-continuation for unlike charges
# ============================================================================

println("="^70)
println("PART 1: 2-Body c-Continuation (Unlike Charges)")
println("="^70)

sys2 = HamiltonianSystem(2, 2)
masses = [1.0, 1.0]
charges = [1.0, -1.0]

c_values = [100.0, 50.0, 20.0, 10.0, 5.0, 4.0, 3.0, 2.0, 1.5, 1.0]

# Store all results
struct ContinuationResult
    family::String
    c::Float64
    T::Float64
    E::Float64
    max_floquet::Float64
    drift_pct::Float64
    status::Symbol
    eigenvalue_types::String
end

all_results = ContinuationResult[]

# --- Circular orbit seeds ---
circular_energies = [-0.5, -1.0, -2.0]

for E_target in circular_energies
    r0 = abs(charges[1]*charges[2]) / (2*abs(E_target))
    q0, p0, T_kepler, E0 = circular_ic_2d(r0)
    family = @sprintf("circular_E%.1f", E_target)

    println("\n--- Family: $family  (r0=$(@sprintf("%.4f",r0)), T_kepler=$(@sprintf("%.4f",T_kepler))) ---")

    T_current = T_kepler
    q_current = copy(q0)
    p_current = copy(p0)

    for c in c_values
        # For circular orbits with ṙ=0, the Weber correction vanishes identically.
        # The orbit remains exact at all c. But the linearized stability changes.
        # We integrate for one period and check closure + monodromy.

        dt_hint = min(1e-3, T_current/500)
        sol = integrate_2body(sys2, q_current, p_current, T_current, masses, charges, c; dt_hint=dt_hint)

        if sol.retcode != :Success
            println("  c=$(@sprintf("%6.1f",c))  FAILED (retcode=$(sol.retcode))")
            push!(all_results, ContinuationResult(family, c, NaN, E_target, NaN, NaN, :Failure, ""))
            continue
        end

        # Measure closure error
        q_end = sol.q[end]
        p_end = sol.p[end]
        closure = norm(vcat(q_end - q_current, p_end - p_current))

        drift = energy_drift(sol)

        # Monodromy
        M = monodromy_2body(sys2, q_current, p_current, T_current, masses, charges, c; dt_hint=dt_hint)
        max_floq = NaN
        eig_types = ""
        if M !== nothing
            ev = eigvals(M)
            max_floq = maximum(abs, ev)
            types = [classify_eigenvalue(e) for e in ev]
            eig_types = join(string.(types), ";")
        end

        status = closure < 1e-4 ? :Closed : :Open
        push!(all_results, ContinuationResult(family, c, T_current, E_target, max_floq, drift, status, eig_types))

        @printf("  c=%6.1f  T=%8.4f  closure=%.2e  drift=%.4f%%  |λ|max=%.4f  types=%s\n",
                c, T_current, closure, isnan(drift) ? 0.0 : drift,
                isnan(max_floq) ? 0.0 : max_floq, eig_types)
    end
end

# --- Elliptical orbit seeds ---
elliptical_cases = [(E=-1.0, e=0.3), (E=-1.0, e=0.5), (E=-1.0, e=0.7)]

for case in elliptical_cases
    E_target = case.E
    ecc = case.e
    q0, p0, T_kepler, E0 = elliptical_ic_2d(E_target, ecc)
    family = @sprintf("elliptic_E%.1f_e%.1f", E_target, ecc)

    println("\n--- Family: $family  (T_kepler=$(@sprintf("%.4f",T_kepler))) ---")

    T_current = T_kepler
    q_current = copy(q0)
    p_current = copy(p0)

    for c in c_values
        dt_hint = min(5e-4, T_current/500)
        sol = integrate_2body(sys2, q_current, p_current, T_current, masses, charges, c; dt_hint=dt_hint)

        if sol.retcode != :Success
            println("  c=$(@sprintf("%6.1f",c))  FAILED")
            push!(all_results, ContinuationResult(family, c, NaN, E_target, NaN, NaN, :Failure, ""))
            continue
        end

        q_end = sol.q[end]
        p_end = sol.p[end]
        closure = norm(vcat(q_end - q_current, p_end - p_current))
        drift = energy_drift(sol)

        # For elliptical orbits, try a simple period correction via bisection on closure
        # The Weber precession means the Kepler period doesn't close the orbit.
        # We search for the actual period near T_kepler.
        best_T = T_current
        best_closure = closure
        if closure > 1e-4 && c < 100
            # Scan around T_kepler for minimum closure
            for dT_frac in range(-0.3, 0.3, length=61)
                T_try = T_kepler * (1.0 + dT_frac)
                T_try < 0.1 && continue
                sol_try = integrate_2body(sys2, q0, p0, T_try, masses, charges, c; dt_hint=min(5e-4, T_try/300))
                sol_try.retcode != :Success && continue
                cl = norm(vcat(sol_try.q[end] - q0, sol_try.p[end] - p0))
                if cl < best_closure
                    best_closure = cl
                    best_T = T_try
                end
            end
            # Refine with finer scan if we found something promising
            if best_closure < 0.5
                T_lo = best_T * 0.99
                T_hi = best_T * 1.01
                for T_try in range(T_lo, T_hi, length=41)
                    sol_try = integrate_2body(sys2, q0, p0, T_try, masses, charges, c; dt_hint=min(5e-4, T_try/300))
                    sol_try.retcode != :Success && continue
                    cl = norm(vcat(sol_try.q[end] - q0, sol_try.p[end] - p0))
                    if cl < best_closure
                        best_closure = cl
                        best_T = T_try
                    end
                end
            end
        end

        # Monodromy at the best period
        M = monodromy_2body(sys2, q0, p0, best_T, masses, charges, c; dt_hint=dt_hint)
        max_floq = NaN
        eig_types = ""
        if M !== nothing
            ev = eigvals(M)
            max_floq = maximum(abs, ev)
            types = [classify_eigenvalue(e) for e in ev]
            eig_types = join(string.(types), ";")
        end

        status = best_closure < 1e-3 ? :Closed : (best_closure < 0.1 ? :NearClosed : :Open)
        push!(all_results, ContinuationResult(family, c, best_T, E_target, max_floq, drift, status, eig_types))

        @printf("  c=%6.1f  T=%8.4f  closure=%.2e  drift=%.4f%%  |λ|max=%.6f  status=%s\n",
                c, best_T, best_closure, isnan(drift) ? 0.0 : drift,
                isnan(max_floq) ? 0.0 : max_floq, status)
    end
end

# ============================================================================
# Part 2: 4-Body c-continuation of known candidates
# ============================================================================

println("\n" * "="^70)
println("PART 2: 4-Body c-Continuation")
println("="^70)

sys4 = HamiltonianSystem(4, 2)
masses4 = fill(1.0, 4)
charges4 = [1.0, 1.0, -1.0, -1.0]

"""Breathing alternating square IC."""
function breathing_square_ic(side, vrad)
    s = side / 2
    X = [s s; -s -s; -s s; s -s]
    q0 = Float64[]
    for i in 1:4
        push!(q0, X[i,1]); push!(q0, X[i,2])
    end
    p0 = zeros(8)
    if vrad != 0.0
        for i in 1:4
            ix = 2*(i-1)+1; iy = ix+1
            nrm = hypot(q0[ix], q0[iy])
            p0[ix] = vrad * q0[ix] / nrm
            p0[iy] = vrad * q0[iy] / nrm
        end
    end
    return q0, p0
end

"""Double orbiter IC: two +/- dimers orbiting each other."""
function double_orbiter_ic(R, v_orbit; dyad_length=0.3)
    a = dyad_length
    q0 = [R/2, a/2, -R/2, a/2, R/2, -a/2, -R/2, -a/2]
    mu_pair = 0.5
    v_intra = sqrt(abs(1.0)/(mu_pair * a)) * 0.4
    p0 = [
         v_orbit,  v_intra,   # particle 1
        -v_orbit, -v_intra,   # particle 2
         v_orbit, -v_intra,   # particle 3
        -v_orbit,  v_intra,   # particle 4
    ]
    return q0, p0
end

# --- Breathing square at various c ---
println("\n--- Breathing alternating square ---")
for vrad in (0.3, 0.5, 0.7)
    for c in [1.0, 2.0, 4.0, 10.0]
        q0, p0 = breathing_square_ic(1.0, vrad)
        dt = 5e-4
        tmax = 15.0
        prob = HamiltonianProblem(sys4, (0.0, tmax), q0, p0;
            masses=masses4, charges=charges4, c=c, dt=dt)
        sol = solve(prob, SymmetricProjectionIntegrator())

        if sol.retcode != :Success
            @printf("  vrad=%.1f  c=%5.1f  FAILED (retcode=%s)\n", vrad, c, sol.retcode)
            push!(all_results, ContinuationResult("breathing_sq_v$(vrad)", c, NaN, NaN, NaN, NaN, :Failure, ""))
            continue
        end

        # Check for brake return (all momenta near zero)
        T_brake = NaN
        for i in max(100, length(sol.t)÷10):length(sol.t)-1
            pmag = maximum(abs, sol.p[i])
            if pmag < 0.05 && maximum(abs, sol.p[i-1]) > pmag && maximum(abs, sol.p[i+1]) > pmag
                T_brake = 2.0 * (sol.t[i] - sol.t[1])
                break
            end
        end

        drift = energy_drift(sol)

        # Check boundedness
        d0 = maximum(hypot(q0[2i-1], q0[2i]) for i in 1:4)
        dmax = maximum(maximum(hypot(sol.q[k][2i-1], sol.q[k][2i]) for i in 1:4) for k in 1:length(sol.t))
        bounded = dmax < 10*d0

        @printf("  vrad=%.1f  c=%5.1f  T_brake=%8.4f  drift=%.4f%%  bounded=%s  dmax/d0=%.2f\n",
                vrad, c, isnan(T_brake) ? -1.0 : T_brake,
                isnan(drift) ? 0.0 : drift, bounded, dmax/d0)

        push!(all_results, ContinuationResult("breathing_sq_v$(vrad)", c,
            isnan(T_brake) ? NaN : T_brake, NaN, NaN,
            isnan(drift) ? NaN : drift, bounded ? :Bounded : :Unbound, ""))
    end
end

# --- Double orbiter survival at various c ---
println("\n--- Double orbiter survival time ---")
for c in [1.0, 2.0, 4.0, 10.0]
    q0, p0 = double_orbiter_ic(3.0, 0.1)
    tmax = min(200.0, 50.0 * c^2)  # test t* ~ c^2 prediction
    dt = min(2e-4, 0.01)
    prob = HamiltonianProblem(sys4, (0.0, tmax), q0, p0;
        masses=masses4, charges=charges4, c=c, dt=dt)
    sol = solve(prob, SymmetricProjectionIntegrator())

    if sol.retcode != :Success
        @printf("  c=%5.1f  FAILED at t=%.2f (retcode=%s)\n", c, sol.t[end], sol.retcode)
        push!(all_results, ContinuationResult("double_orbiter", c, NaN, NaN, NaN, NaN, :Failure, ""))
        continue
    end

    # Check when it becomes unbound (any pair distance > 20)
    t_unbound = NaN
    for k in 1:length(sol.t)
        q = sol.q[k]
        dmax = 0.0
        for i in 1:4, j in i+1:4
            dx = q[2i-1] - q[2j-1]; dy = q[2i] - q[2j]
            dmax = max(dmax, sqrt(dx^2+dy^2))
        end
        if dmax > 20.0
            t_unbound = sol.t[k]
            break
        end
    end

    drift = energy_drift(sol)
    t_survive = isnan(t_unbound) ? sol.t[end] : t_unbound

    @printf("  c=%5.1f  t_survive=%8.2f  t_max=%8.2f  drift=%.4f%%  (t*/c^2=%.2f)\n",
            c, t_survive, tmax, isnan(drift) ? 0.0 : drift, t_survive/c^2)

    push!(all_results, ContinuationResult("double_orbiter", c, t_survive, NaN, NaN,
        isnan(drift) ? NaN : drift, isnan(t_unbound) ? :Bounded : :Unbound, "t*/c^2=$(@sprintf("%.2f",t_survive/c^2))"))
end

# ============================================================================
# Part 3: Detailed bifurcation analysis for 2-body circular orbits
# ============================================================================

println("\n" * "="^70)
println("PART 3: Bifurcation Analysis — Eigenvalue Tracking")
println("="^70)

# Fine c-grid near predicted bifurcation points
# From Agent 04: c_bif = sqrt(k/(2*mu*a*(n-1))) where k=1, mu=0.5
# For E=-1.0 (a=0.5): c_bif(1) = sqrt(1/(2*0.5*0.5*1)) = sqrt(2) ≈ 1.414
# For E=-0.5 (a=1.0): c_bif(1) = sqrt(1/(2*0.5*1.0*1)) = 1.0

fine_c_grids = Dict(
    -0.5 => sort(vcat(collect(range(0.5, 2.0, length=31)), [0.95, 1.0, 1.05])),
    -1.0 => sort(vcat(collect(range(0.8, 3.0, length=45)), [1.35, 1.40, 1.414, 1.45, 1.50])),
    -2.0 => sort(vcat(collect(range(0.5, 5.0, length=46)), [0.95, 1.0, 1.414, 2.0])),
)

struct BifurcationPoint
    E::Float64
    c_approx::Float64
    type::String
    description::String
end

bifurcations_found = BifurcationPoint[]

for E_target in [-0.5, -1.0, -2.0]
    r0 = abs(charges[1]*charges[2]) / (2*abs(E_target))
    q0, p0, T_kepler, _ = circular_ic_2d(r0)
    family = @sprintf("circular_E%.1f_bif", E_target)

    println("\n--- Bifurcation scan: E = $E_target (r0=$(@sprintf("%.4f",r0))) ---")

    c_grid = fine_c_grids[E_target]
    prev_max_floq = NaN
    prev_eigenvalues = nothing

    for c in c_grid
        dt_hint = min(5e-4, T_kepler/500)
        M = monodromy_2body(sys2, q0, p0, T_kepler, masses, charges, c; eps_fd=1e-7, dt_hint=dt_hint)
        if M === nothing
            continue
        end

        ev = eigvals(M)
        max_floq = maximum(abs, ev)

        # Check for eigenvalue crossing +1 or -1
        for e in ev
            if abs(e - 1.0) < 0.05 && abs(imag(e)) < 0.01
                @printf("  c=%6.3f  eigenvalue near +1: %.4f + %.4fi\n", c, real(e), imag(e))
            end
            if abs(e + 1.0) < 0.05 && abs(imag(e)) < 0.01
                @printf("  c=%6.3f  eigenvalue near -1: %.4f + %.4fi\n", c, real(e), imag(e))
            end
        end

        # Detect stability transition
        if !isnan(prev_max_floq) && ((prev_max_floq < 1.01 && max_floq > 1.01) ||
                                      (prev_max_floq > 1.01 && max_floq < 1.01))
            @printf("  *** STABILITY TRANSITION at c ≈ %.3f: |λ|max %.4f -> %.4f\n",
                    c, prev_max_floq, max_floq)
            trans_type = max_floq > 1.01 ? "stable→unstable" : "unstable→stable"
            push!(bifurcations_found, BifurcationPoint(E_target, c, trans_type,
                "max|λ|: $(@sprintf("%.4f",prev_max_floq)) → $(@sprintf("%.4f",max_floq))"))
        end

        # Track eigenvalue collision on unit circle (Krein collision)
        if prev_eigenvalues !== nothing
            # Check if two eigenvalues on the unit circle have merged
            on_circle_prev = [e for e in prev_eigenvalues if abs(abs(e) - 1.0) < 0.05 && abs(imag(e)) > 0.01]
            on_circle_now = [e for e in ev if abs(abs(e) - 1.0) < 0.05 && abs(imag(e)) > 0.01]
            if length(on_circle_prev) >= 4 && length(on_circle_now) < length(on_circle_prev)
                @printf("  *** POSSIBLE KREIN COLLISION at c ≈ %.3f\n", c)
                push!(bifurcations_found, BifurcationPoint(E_target, c, "Krein collision",
                    "circle eigenvalues: $(length(on_circle_prev)) → $(length(on_circle_now))"))
            end
        end

        prev_max_floq = max_floq
        prev_eigenvalues = ev

        push!(all_results, ContinuationResult(family, c, T_kepler, E_target, max_floq, NaN, :Computed, ""))
    end
end

# ============================================================================
# Part 4: New orbit search at each c
# ============================================================================

println("\n" * "="^70)
println("PART 4: New Orbit Search (T-brake shooting)")
println("="^70)

# For 2-body: try T-brake orbits with various breathing amplitudes
# A "brake" orbit starts with p=0 at some configuration and returns to p=0

for c in [1.0, 2.0, 5.0, 10.0]
    println("\n--- c = $c ---")
    for r_start in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
        # Start at rest (brake) along x-axis
        q0 = [-r_start/2, 0.0, r_start/2, 0.0]
        p0 = [0.0, 0.0, 0.0, 0.0]

        # Head-on collision: particles fall toward each other then repel (for unlike charges, they attract)
        # For unlike charges, head-on means they fall together. Need collision bounce.
        # Instead, give a small transverse kick to avoid singularity
        for L_kick in [0.05, 0.2, 0.5]
            p_kick = [0.0, -L_kick/r_start, 0.0, L_kick/r_start]  # angular momentum = L_kick
            tmax = 20.0
            dt = min(5e-4, 0.001)
            sol = integrate_2body(sys2, q0, p_kick, tmax, masses, charges, c; dt_hint=dt)
            sol.retcode != :Success && continue

            # Look for periodic return
            best_closure = Inf
            best_T = NaN
            for k in max(50, length(sol.t)÷20):length(sol.t)
                cl = norm(vcat(sol.q[k] - q0, sol.p[k] - p_kick))
                if cl < best_closure
                    best_closure = cl
                    best_T = sol.t[k]
                end
            end

            if best_closure < 0.1
                @printf("  c=%5.1f  r0=%.1f  L=%.2f  T_approx=%8.4f  closure=%.2e  (CANDIDATE)\n",
                        c, r_start, L_kick, best_T, best_closure)
                push!(all_results, ContinuationResult("new_orbit_search", c, best_T, NaN, NaN, NaN,
                    :Candidate, "r0=$r_start, L=$L_kick"))
            end
        end
    end
end

# ============================================================================
# Write CSV output
# ============================================================================

csv_path = joinpath(OUTDIR, "bifurcation_diagram.csv")
open(csv_path, "w") do io
    println(io, "family,c,T,E,max_floquet,drift_pct,status,eigenvalue_types")
    for r in all_results
        @printf(io, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%s\n",
            r.family, r.c,
            isnan(r.T) ? -1.0 : r.T,
            isnan(r.E) ? -1.0 : r.E,
            isnan(r.max_floquet) ? -1.0 : r.max_floquet,
            isnan(r.drift_pct) ? -1.0 : r.drift_pct,
            r.status,
            r.eigenvalue_types)
    end
end
println("\nWrote $(length(all_results)) rows to $csv_path")

# ============================================================================
# Summary of bifurcations
# ============================================================================

println("\n" * "="^70)
println("BIFURCATION SUMMARY")
println("="^70)

if isempty(bifurcations_found)
    println("No bifurcations detected in the scanned c-range.")
else
    for b in bifurcations_found
        @printf("  E=%.1f  c≈%.3f  type=%s  (%s)\n", b.E, b.c_approx, b.type, b.description)
    end
end

println("\nDone.")
