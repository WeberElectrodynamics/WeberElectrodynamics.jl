"""
Reeb flow integrator for the 2-body Weber energy surface.

The Reeb field on a contact-type energy surface Σ_E is R_α = X_H / (Y·H)
where Y·H = T = kinetic energy = E - V(q). This is the Hamiltonian flow
with time reparameterized by 1/T.

For 2-body unlike charges at negative energies, the Hill region is an annulus
and the energy surface is contact type (proven in NOTES.md). We integrate the
Reeb flow and detect closed orbits (= periodic Hamiltonian orbits) via
Poincare section crossings.
"""

using LinearAlgebra, Printf
using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator

# ============================================================================
# Hamiltonian integration and Reeb reparameterization
# ============================================================================

"""
Build compiled Weber Hamiltonian for 2-body 2D, returning
(sys, H, V, params, dof) where H(q,p) and V(q) are callables.
"""
function build_2body(charges, masses, c)
    sys = HamiltonianSystem(2, 2)
    kappas = [1.0]
    params = vcat(masses, charges, [c], kappas)
    dof = 4
    zero_p = zeros(dof)

    H = (q, p) -> sys.hamiltonian_compiled(q, p, params)
    V = q -> sys.hamiltonian_compiled(q, zero_p, params)

    return (sys=sys, H=H, V=V, params=params, dof=dof)
end

"""
Compute kinetic energy T = E - V(q) = H - V at a phase-space point.
"""
function kinetic_energy(H_func, V_func, q, p)
    return H_func(q, p) - V_func(q)
end

"""
    integrate_reeb_2body(charges, masses, c, q0, p0, E, reeb_time;
                         ham_dt=1e-4, n_steps=nothing)

Integrate the Reeb flow on Σ_E by integrating the Hamiltonian flow and
reparameterizing time. Returns arrays of (reeb_τ, ham_t, q, p, T_kinetic).

The Reeb time dτ = T(t) dt, so τ = ∫₀ᵗ T(s) ds. We integrate in
Hamiltonian time and record the cumulative Reeb time.
"""
function integrate_reeb_2body(charges, masses, c, q0, p0, E;
                              tmax::Float64=200.0, dt::Float64=1e-4)
    wb = build_2body(charges, masses, c)

    # Verify energy
    H0 = wb.H(q0, p0)
    @assert abs(H0 - E) < 1e-6 "Initial condition not on energy surface: H=$H0, E=$E"

    # Integrate with the Weber integrator
    prob = HamiltonianProblem(wb.sys, (0.0, tmax), q0, p0;
                        masses=masses, charges=charges, c=c, dt=dt)
    sol = solve(prob, SymmetricProjectionIntegrator())

    # Extract trajectory and compute Reeb time
    n_out = length(sol.t)
    reeb_times = zeros(n_out)
    T_kinetic = zeros(n_out)

    for i in 1:n_out
        qi = sol.q[i]
        pi = sol.p[i]
        T_kinetic[i] = wb.H(qi, pi) - wb.V(qi)
    end

    # Cumulative Reeb time: τ = ∫ T dt via trapezoidal rule
    for i in 2:n_out
        Δt = sol.t[i] - sol.t[i-1]
        reeb_times[i] = reeb_times[i-1] + 0.5 * (T_kinetic[i] + T_kinetic[i-1]) * Δt
    end

    return (t=sol.t, τ=reeb_times, q=sol.q, p=sol.p, T=T_kinetic,
            H_err=[wb.H(sol.q[i], sol.p[i]) - E for i in 1:n_out])
end

# ============================================================================
# Closed orbit detection via Poincare section
# ============================================================================

"""
Relative coordinates: r = q1 - q2, p_r = μ(v1 - v2).
For 2D: r = (q[1]-q[3], q[2]-q[4]), etc.
"""
function relative_coords(q, p, masses)
    m1, m2 = masses
    μ = m1 * m2 / (m1 + m2)
    r = [q[1] - q[3], q[2] - q[4]]
    v_rel = [p[1]/m1 - p[3]/m2, p[2]/m1 - p[4]/m2]
    p_r = μ .* v_rel
    return r, p_r
end

"""
    detect_closed_orbits(sol_data, masses; section_angle=0.0, tol=0.02)

Detect closed Reeb orbits by looking for near-recurrences in the Poincare
section defined by angle(r) = section_angle.

Returns list of (period_ham, period_reeb, energy_error, recurrence_dist).
"""
function detect_closed_orbits(sol_data, masses; section_tol=0.05, recurrence_tol=0.05)
    t, τ, qs, ps = sol_data.t, sol_data.τ, sol_data.q, sol_data.p

    # Collect Poincare section crossings: y_rel = 0, x_rel > 0, vy_rel > 0
    crossings = Int[]
    for i in 2:length(t)
        r_prev, _ = relative_coords(qs[i-1], ps[i-1], masses)
        r_curr, _ = relative_coords(qs[i], ps[i], masses)
        # y = 0 crossing with x > 0
        if r_prev[2] * r_curr[2] < 0 && r_curr[1] > 0
            push!(crossings, i)
        end
    end

    if length(crossings) < 2
        return []
    end

    # For each pair of crossings, check recurrence in (r, p_r) space
    orbits = NamedTuple{(:ham_period, :reeb_period, :H_err_max, :recurrence_dist), NTuple{4,Float64}}[]

    # First crossing as reference
    ref_idx = crossings[1]
    r_ref, p_ref = relative_coords(qs[ref_idx], ps[ref_idx], masses)
    state_ref = vcat(r_ref, p_ref)

    for k in 2:length(crossings)
        idx = crossings[k]
        r_k, p_k = relative_coords(qs[idx], ps[idx], masses)
        state_k = vcat(r_k, p_k)
        dist = norm(state_k - state_ref)

        if dist < recurrence_tol
            T_ham = t[idx] - t[ref_idx]
            T_reeb = τ[idx] - τ[ref_idx]
            H_err = maximum(abs.(sol_data.H_err[ref_idx:idx]))
            push!(orbits, (ham_period=T_ham, reeb_period=T_reeb,
                           H_err_max=H_err, recurrence_dist=dist))
            break  # Take the first recurrence (shortest period)
        end
    end

    return orbits
end

# ============================================================================
# Stability analysis via monodromy approximation
# ============================================================================

"""
Estimate stability of a periodic orbit by checking if nearby trajectories
remain close (stable) or diverge (unstable) over one period.
"""
function estimate_stability(charges, masses, c, q0, p0, E, period;
                           dt::Float64=1e-4, ε::Float64=1e-5)
    wb = build_2body(charges, masses, c)

    # Reference orbit
    prob0 = HamiltonianProblem(wb.sys, (0.0, period * 1.2), q0, p0;
                         masses=masses, charges=charges, c=c, dt=dt)
    sol0 = solve(prob0, SymmetricProjectionIntegrator())

    # Perturb in relative py direction (transverse to orbit)
    p_pert = copy(p0)
    p_pert[2] += ε  # perturb py1
    p_pert[4] -= ε  # keep total momentum zero

    # Rescale to match energy
    H_pert = wb.H(q0, p_pert)
    # Simple energy correction: scale all momenta
    if abs(H_pert - E) > 1e-10
        # Approximate: T ≈ p²/(2μ), so scale p by sqrt(T_target/T_actual)
        T_target = E - wb.V(q0)
        T_actual = H_pert - wb.V(q0)
        if T_actual > 0 && T_target > 0
            scale = sqrt(T_target / T_actual)
            p_pert .*= scale
        end
    end

    prob1 = HamiltonianProblem(wb.sys, (0.0, period * 1.2), q0, p_pert;
                         masses=masses, charges=charges, c=c, dt=dt)
    sol1 = solve(prob1, SymmetricProjectionIntegrator())

    # Compare at t ≈ period
    idx0 = argmin(abs.(sol0.t .- period))
    idx1 = argmin(abs.(sol1.t .- period))

    δq = norm(sol1.q[idx1] - sol0.q[idx0])
    δp = norm(sol1.p[idx1] - sol0.p[idx0])

    # Lyapunov-like exponent estimate
    δ0 = ε
    δf = sqrt(δq^2 + δp^2)
    λ_est = log(δf / δ0) / period

    return (δq=δq, δp=δp, λ_estimate=λ_est,
            stability = λ_est < 0.1 ? :stable : :unstable)
end

# ============================================================================
# Main: scan energies for 2-body unlike charges
# ============================================================================

function run_reeb_scan()
    charges = [1.0, -1.0]
    masses = [1.0, 1.0]
    c = 1.0
    μ = 0.5

    println("="^70)
    println("REEB FLOW INTEGRATION: 2-body unlike charges (+1,-1)")
    println("="^70)
    println()

    results = []

    for E in [-5.0, -3.0, -2.0, -1.0, -0.5, -0.2, -0.1]
        @printf "\n--- Energy E = %.2f ---\n" E

        # For circular orbit: V(r) = -1/r = E => r_circ = 1/|E|
        # But we need T > 0, so use an orbit with some kinetic energy.
        # For circular Coulomb orbit: E = -q1*q2/(2r) = 1/(2r) (since q1*q2=-1)
        # Wait: V(r) = q1*q2/r = -1/r, so for circular Coulomb:
        # E = T + V = V/2 = -1/(2r) => r = 1/(2|E|)
        r_circ = 1.0 / (2.0 * abs(E))

        # Circular orbit: p_θ = μ * r * ω, with v_circ = sqrt(|q1q2|/(μr)) for Coulomb
        # T = 1/(2r), V = -1/r, E = -1/(2r)
        v_circ = sqrt(1.0 / (μ * r_circ))  # |F|/μ = |q1q2|/r² = 1/r²; v²/r = 1/(μr²)
        p_circ = μ * v_circ

        # Place particles in COM: particle 1 at (r/2, 0), particle 2 at (-r/2, 0)
        # Circular orbit: momenta tangential
        q0 = [r_circ/2, 0.0, -r_circ/2, 0.0]
        p0 = [0.0, p_circ/2, 0.0, -p_circ/2]  # py only, tangential

        # Check energy
        wb = build_2body(charges, masses, c)
        H0 = wb.H(q0, p0)
        @printf "Circular orbit: r=%.4f, v=%.4f, H=%.6f (target E=%.6f)\n" r_circ v_circ H0 E

        # Adjust to hit exact energy by scaling momenta
        V0 = wb.V(q0)
        T_needed = E - V0
        if T_needed <= 0
            @printf "  Cannot achieve E=%.2f at this r (V=%.4f > E), trying wider orbit\n" E V0
            # Use r such that V(r) < E
            r_try = 1.0 / abs(E) + 0.1
            q0 = [r_try/2, 0.0, -r_try/2, 0.0]
            V0 = wb.V(q0)
            T_needed = E - V0
            if T_needed <= 0
                @printf "  Still cannot achieve energy, skipping\n"
                continue
            end
        end

        # Scale momenta to match energy exactly
        T_current = H0 - V0
        if T_current > 0
            scale = sqrt(T_needed / T_current)
            p0 = p0 .* scale
        end

        H_check = wb.H(q0, p0)
        @printf "  After scaling: H=%.8f, |H-E|=%.2e\n" H_check abs(H_check - E)

        # Estimate orbital period (Kepler: T = 2π a^(3/2) / sqrt(GM))
        # For Coulomb with μ: a = |q1q2|/(2|E|) = 1/(2|E|), T_kepler = 2π μ a² / L
        # For circular: T = 2πr/v
        T_est = 2π * r_circ / v_circ
        @printf "  Estimated period: %.4f\n" T_est

        # Integrate for several periods
        tmax = min(T_est * 8, 500.0)
        dt = min(T_est / 500, 1e-3)

        sol_data = integrate_reeb_2body(charges, masses, c, q0, p0, E;
                                        tmax=tmax, dt=dt)

        @printf "  Integration: %d steps, t_final=%.2f, τ_final=%.4f\n" length(sol_data.t) sol_data.t[end] sol_data.τ[end]
        @printf "  Max |ΔH|: %.2e\n" maximum(abs.(sol_data.H_err))
        @printf "  T_kinetic range: [%.4f, %.4f]\n" minimum(sol_data.T) maximum(sol_data.T)

        # Detect closed orbits
        orbits = detect_closed_orbits(sol_data, masses; recurrence_tol=0.1)

        if isempty(orbits)
            @printf "  No closed orbit detected (recurrence_tol=0.1)\n"
            # Report first few Poincare section crossings for diagnostics
        else
            orb = orbits[1]
            @printf "  CLOSED ORBIT FOUND:\n"
            @printf "    Hamiltonian period: %.6f\n" orb.ham_period
            @printf "    Reeb period:        %.6f\n" orb.reeb_period
            @printf "    Max energy error:   %.2e\n" orb.H_err_max
            @printf "    Recurrence dist:    %.2e\n" orb.recurrence_dist

            # Stability
            stab = estimate_stability(charges, masses, c, q0, p0, E, orb.ham_period; dt=dt)
            @printf "    Stability: %s (λ ≈ %.4f)\n" stab.stability stab.λ_estimate

            push!(results, (E=E, T_ham=orb.ham_period, T_reeb=orb.reeb_period,
                           H_err=orb.H_err_max, stability=stab.stability,
                           λ=stab.λ_estimate))
        end
    end

    println()
    println("="^70)
    println("SUMMARY OF CLOSED REEB ORBITS")
    println("="^70)
    @printf "%-8s %-12s %-12s %-12s %-10s %-10s\n" "E" "T_ham" "T_reeb" "H_err" "stable?" "λ"
    for r in results
        @printf "%-8.2f %-12.4f %-12.4f %-12.2e %-10s %-10.4f\n" r.E r.T_ham r.T_reeb r.H_err r.stability r.λ
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_reeb_scan()
end
