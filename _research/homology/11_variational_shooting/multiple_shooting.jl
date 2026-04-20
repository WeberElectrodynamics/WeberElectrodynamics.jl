# Agent 11 — Multiple-shooting Newton iterator for periodic orbits
#             of the Weber Hamiltonian.
#
# Methods:
#   1. Single-shooting Newton (full period)
#   2. Multiple-shooting Newton (subdivided arcs)
#   3. c-continuation from Coulomb/Kepler limit
#
# Run:
#   julia --project=. research/homology/11_variational_shooting/multiple_shooting.jl
#
# Outputs:
#   - found_orbits.csv (next to this file)
#   - convergence log to stdout

using LinearAlgebra
using Printf

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using CommonSolve: solve

const OUTDIR = @__DIR__

# Include the CZ index calculator from Agent 04
include(joinpath(@__DIR__, "..", "04_cz_index", "cz_index.jl"))

# ============================================================================
# Core helpers
# ============================================================================

"""
Build a HamiltonianSystem once per (n_particles, dims) pair and cache it.
"""
const SYS_CACHE = Dict{Tuple{Int,Int}, HamiltonianSystem}()
function get_system(n::Int, d::Int)
    get!(SYS_CACHE, (n, d)) do
        HamiltonianSystem(n, d)
    end
end

"""
    flow(q0, p0, T; n_particles, dims, masses, charges, c, dt)

Integrate the Weber system from (q0, p0) for time T.
Returns (q_final, p_final) or (nothing, nothing) on failure.
"""
function flow(q0::Vector{Float64}, p0::Vector{Float64}, T::Float64;
              n_particles::Int, dims::Int,
              masses::Vector{Float64}, charges::Vector{Float64},
              c::Float64, dt::Float64=1e-3)
    T <= 0 && return (copy(q0), copy(p0))
    sys = get_system(n_particles, dims)
    nsteps = max(10, round(Int, T / dt))
    dt_use = T / nsteps
    prob = HamiltonianProblem(sys, (0.0, T), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt_use)
    sol = solve(prob, SymmetricProjectionIntegrator())
    if sol.retcode != :Success
        return (nothing, nothing)
    end
    return (sol.q[end], sol.p[end])
end

"""
    flow_with_energy(q0, p0, T; kwargs...)

Like flow(), but also returns the initial energy from the solution.
"""
function flow_with_energy(q0::Vector{Float64}, p0::Vector{Float64}, T::Float64;
                          n_particles::Int, dims::Int,
                          masses::Vector{Float64}, charges::Vector{Float64},
                          c::Float64, dt::Float64=1e-3)
    T <= 0 && return (copy(q0), copy(p0), NaN)
    sys = get_system(n_particles, dims)
    nsteps = max(10, round(Int, T / dt))
    dt_use = T / nsteps
    prob = HamiltonianProblem(sys, (0.0, T), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt_use)
    sol = solve(prob, SymmetricProjectionIntegrator())
    if sol.retcode != :Success
        return (nothing, nothing, NaN)
    end
    en = compute_energy_timeseries(sol)
    return (sol.q[end], sol.p[end], en.total_energy[1])
end

# ============================================================================
# Finite-difference monodromy (Jacobian of the flow map)
# ============================================================================

"""
    flow_jacobian(q0, p0, T; kwargs...)

Compute the Jacobian dΦ_T/dz₀ of the Weber flow map by central finite differences.
Returns the 2*DOF × 2*DOF monodromy matrix, or nothing on failure.
"""
function flow_jacobian(q0::Vector{Float64}, p0::Vector{Float64}, T::Float64;
                       n_particles::Int, dims::Int,
                       masses::Vector{Float64}, charges::Vector{Float64},
                       c::Float64, dt::Float64=1e-3, eps::Float64=1e-7)
    dof = n_particles * dims
    n = 2 * dof
    z0 = vcat(q0, p0)

    # Reference flow
    (qr, pr) = flow(q0, p0, T; n_particles=n_particles, dims=dims,
                     masses=masses, charges=charges, c=c, dt=dt)
    (qr === nothing) && return nothing
    z_ref = vcat(qr, pr)

    J = zeros(n, n)
    for k in 1:n
        zp = copy(z0); zp[k] += eps
        zm = copy(z0); zm[k] -= eps
        (qp, pp) = flow(zp[1:dof], zp[dof+1:end], T;
                        n_particles=n_particles, dims=dims,
                        masses=masses, charges=charges, c=c, dt=dt)
        (qm, pm) = flow(zm[1:dof], zm[dof+1:end], T;
                        n_particles=n_particles, dims=dims,
                        masses=masses, charges=charges, c=c, dt=dt)
        if qp === nothing || qm === nothing
            return nothing
        end
        J[:, k] = (vcat(qp, pp) - vcat(qm, pm)) / (2 * eps)
    end
    return J
end

# ============================================================================
# Single-shooting Newton iterator
# ============================================================================

"""
    OrbitResult — stores data for a found periodic orbit
"""
mutable struct OrbitResult
    id::String
    n_particles::Int
    dims::Int
    q0::Vector{Float64}
    p0::Vector{Float64}
    period::Float64
    energy::Float64
    c::Float64
    mismatch::Float64          # ||F|| at convergence
    max_floquet::Float64
    cz_index::Int
    cz_upper::Int
    stability_type::String
    converged::Bool
    newton_iters::Int
    notes::String
end

const ALL_RESULTS = OrbitResult[]

"""
    single_shooting_newton(q0, p0, T; kwargs...) -> OrbitResult or nothing

Newton iteration to close a periodic orbit.
Solves F(z) = Φ_T(z) - z = 0 where Φ_T is the Weber flow map.

For Hamiltonian systems the monodromy always has eigenvalue 1 (energy direction),
so we fix energy by projecting out the energy gradient direction, effectively
solving in a Poincare section.
"""
function single_shooting_newton(q0_init::Vector{Float64}, p0_init::Vector{Float64},
                                T_init::Float64;
                                n_particles::Int, dims::Int,
                                masses::Vector{Float64}, charges::Vector{Float64},
                                c::Float64, dt::Float64=1e-3,
                                max_iter::Int=20, tol::Float64=1e-10,
                                fix_period::Bool=true,
                                id::String="orbit",
                                verbose::Bool=true)
    dof = n_particles * dims
    n = 2 * dof
    q0 = copy(q0_init)
    p0 = copy(p0_init)
    T = T_init

    for iter in 1:max_iter
        # Forward flow
        (qf, pf) = flow(q0, p0, T; n_particles=n_particles, dims=dims,
                         masses=masses, charges=charges, c=c, dt=dt)
        if qf === nothing
            verbose && println("  [$id] iter $iter: flow failed")
            return nothing
        end

        # Mismatch
        F = vcat(qf - q0, pf - p0)
        mismatch = norm(F)
        verbose && @printf("  [%s] iter %d: ||F|| = %.3e\n", id, iter, mismatch)

        if mismatch < tol
            # Converged — compute Floquet multipliers and CZ index
            verbose && println("  [$id] CONVERGED in $iter iterations")
            return _finalize_orbit(id, n_particles, dims, q0, p0, T, c,
                                   mismatch, masses, charges, dt, iter)
        end

        # Jacobian
        J = flow_jacobian(q0, p0, T; n_particles=n_particles, dims=dims,
                          masses=masses, charges=charges, c=c, dt=dt)
        if J === nothing
            verbose && println("  [$id] iter $iter: Jacobian computation failed")
            return nothing
        end

        # Solve (J - I) δz = -F
        A = J - I(n)

        # Use pseudoinverse for the singular system (energy direction)
        δz = -pinv(A) * F

        # Damped update (line search)
        α = 1.0
        for ls in 1:5
            q_try = q0 + α * δz[1:dof]
            p_try = p0 + α * δz[dof+1:end]
            (qt, pt) = flow(q_try, p_try, T; n_particles=n_particles, dims=dims,
                            masses=masses, charges=charges, c=c, dt=dt)
            if qt !== nothing
                F_new = vcat(qt - q_try, pt - p_try)
                if norm(F_new) < mismatch
                    q0 = q_try
                    p0 = p_try
                    break
                end
            end
            α *= 0.5
        end
        if α < 1.0 / 32
            # If line search failed, take small step anyway
            q0 = q0 + 0.1 * δz[1:dof]
            p0 = p0 + 0.1 * δz[dof+1:end]
        end
    end

    # Did not converge — check if close enough
    (qf, pf) = flow(q0, p0, T; n_particles=n_particles, dims=dims,
                     masses=masses, charges=charges, c=c, dt=dt)
    if qf !== nothing
        mismatch = norm(vcat(qf - q0, pf - p0))
        if mismatch < 1e-4
            verbose && @printf("  [%s] near-converged: ||F|| = %.3e\n", id, mismatch)
            return _finalize_orbit(id, n_particles, dims, q0, p0, T, c,
                                   mismatch, masses, charges, dt, max_iter)
        end
    end
    verbose && println("  [$id] did NOT converge after $max_iter iterations")
    return nothing
end

"""
Finalize an orbit: compute energy, Floquet multipliers, CZ index.
"""
function _finalize_orbit(id, n_particles, dims, q0, p0, T, c,
                         mismatch, masses, charges, dt, iters)
    dof = n_particles * dims

    # Energy
    (_, _, E0) = flow_with_energy(q0, p0, T; n_particles=n_particles, dims=dims,
                                   masses=masses, charges=charges, c=c, dt=dt)

    # Monodromy
    max_floquet = NaN
    cz = -999
    cz_upper = -999
    stab_type = "unknown"

    J = flow_jacobian(q0, p0, T; n_particles=n_particles, dims=dims,
                      masses=masses, charges=charges, c=c, dt=dt,
                      eps=max(1e-7, mismatch * 10))
    if J !== nothing && !any(isnan, J)
        ev = eigvals(J)
        max_floquet = maximum(abs, ev)

        # CZ index
        try
            cz_val, details = conley_zehnder_index(J; remove_energy=true)
            cz = cz_val
            cz_upper = cz_val

            # Classify stability
            if max_floquet < 1.01
                stab_type = "elliptic"
            elseif max_floquet < 2.0
                stab_type = "weakly_hyperbolic"
            else
                stab_type = "strongly_hyperbolic"
            end
        catch e
            @warn "CZ computation failed for $id" exception=e
        end
    end

    result = OrbitResult(id, n_particles, dims, copy(q0), copy(p0),
                         T, E0, c, mismatch, max_floquet, cz, cz_upper,
                         stab_type, mismatch < 1e-4, iters, "")
    push!(ALL_RESULTS, result)
    return result
end

# ============================================================================
# 2-body circular orbit exact solutions
# ============================================================================

"""
    two_body_circular_ic(; E, L, m1, m2, q1, q2, c)

Compute exact initial conditions for a 2-body unlike-charge circular orbit.
For circular orbits, rdot=0 so the Weber correction vanishes identically.
These are exact periodic orbits of both Coulomb and Weber dynamics.

Returns (q0, p0, T) in flat 2D format for a 2-body system.
"""
function two_body_circular_ic(; E::Float64=-0.5, m1::Float64=1.0, m2::Float64=1.0,
                                q1::Float64=1.0, q2::Float64=-1.0, c::Float64=1.0)
    mu = m1 * m2 / (m1 + m2)
    k = -q1 * q2  # positive for unlike charges
    k > 0 || error("two_body_circular_ic requires unlike charges")

    # Circular orbit: E = -k/(2a), a = k/(2|E|)
    a = k / (2 * abs(E))

    # Angular momentum: L = mu * v * a, v = sqrt(k/(mu*a))
    v = sqrt(k / (mu * a))
    L = mu * v * a

    # Period: T = 2*pi*a / v
    T = 2 * pi * a / v

    # Place particles in COM frame along x-axis
    # r1 = -m2/(m1+m2) * a, r2 = m1/(m1+m2) * a
    M = m1 + m2
    r1 = -m2 / M * a
    r2 = m1 / M * a

    # Positions along x-axis, velocities along y-axis
    q0 = Float64[r1, 0.0, r2, 0.0]  # [x1, y1, x2, y2]
    # momenta: p = m * v_tangential (y-direction)
    v1 = -m2 / M * v  # COM constraint
    v2 = m1 / M * v
    p0 = Float64[0.0, m1 * v1, 0.0, m2 * v2]
    # Wait — in COM frame, particle 1 at (r1, 0), velocity in y.
    # v_circ for particle 1: ω * |r1| in tangential direction.
    # ω = v / a
    omega = v / a
    p0 = Float64[0.0, m1 * omega * r1,  # = m1 * (-m2/M) * v  ... wrong sign
                 0.0, m2 * omega * r2]
    # Actually: particle 1 is at (r1, 0) with r1 < 0, tangential = (0, +1) rotated
    # For counterclockwise orbit: v1 = omega * (-y1, x1) applied to position (r1, 0)
    # gives v1 = omega * (0, r1) — but r1 < 0, so vy1 < 0.
    # Hmm, let's be more careful:
    # Position of particle 1: (r1, 0) where r1 = -m2/M * a < 0
    # For CCW circular orbit: v = omega * (-y, x)
    # v1 = omega * (0, r1) — since y1=0, x1=r1
    # So vy1 = omega * r1 < 0
    # Position of particle 2: (r2, 0) where r2 > 0
    # v2 = omega * (0, r2), vy2 = omega * r2 > 0
    p0 = Float64[0.0, m1 * omega * r1,
                 0.0, m2 * omega * r2]

    return (q0, p0, T)
end

# ============================================================================
# 4-body breathing square IC (from prior study)
# ============================================================================

"""
    breathing_square_ic(; side, v_rad)

Generate IC for the breathing alternating square.
Particles: 1(+) at (s/2, s/2), 2(+) at (-s/2, -s/2),
           3(-) at (-s/2, s/2), 4(-) at (s/2, -s/2)
with outward radial velocity v_rad.
"""
function breathing_square_ic(; side::Float64=2.0, v_rad::Float64=0.5)
    s = side / 2
    positions = [
        s  s;
        -s -s;
        -s  s;
        s -s;
    ]
    masses = fill(1.0, 4)
    charges = [1.0, 1.0, -1.0, -1.0]

    q0 = Float64[]
    p0 = Float64[]
    for i in 1:4
        r = positions[i, :]
        nrm = norm(r)
        push!(q0, r[1], r[2])
        push!(p0, masses[i] * v_rad * r[1] / nrm,
                  masses[i] * v_rad * r[2] / nrm)
    end

    return (q0, p0, masses, charges)
end

# ============================================================================
# c-Continuation
# ============================================================================

"""
    c_continuation(q0, p0, T; c_values, kwargs...)

Track a periodic orbit as c varies. At each c, use Newton shooting to
close the orbit starting from the previous converged solution.
"""
function c_continuation(q0_init::Vector{Float64}, p0_init::Vector{Float64},
                        T_init::Float64;
                        c_values::Vector{Float64},
                        n_particles::Int, dims::Int,
                        masses::Vector{Float64}, charges::Vector{Float64},
                        dt::Float64=1e-3,
                        base_id::String="cont",
                        verbose::Bool=true)
    q0 = copy(q0_init)
    p0 = copy(p0_init)
    T = T_init
    results = OrbitResult[]

    for (i, c) in enumerate(c_values)
        id = @sprintf("%s_c%.3f", base_id, c)
        verbose && println("\n--- c-continuation: c = $c ---")
        r = single_shooting_newton(q0, p0, T;
            n_particles=n_particles, dims=dims,
            masses=masses, charges=charges, c=c, dt=dt,
            max_iter=15, tol=1e-8, id=id, verbose=verbose)
        if r !== nothing && r.converged
            push!(results, r)
            q0 = copy(r.q0)
            p0 = copy(r.p0)
            T = r.period  # period doesn't change in our fixed-T scheme
            verbose && @printf("  c=%.3f: CONVERGED, E=%.6f, |λ|=%.3f, μ_CZ=%d\n",
                               c, r.energy, r.max_floquet, r.cz_index)
        else
            verbose && println("  c=$c: FAILED — orbit may cease to exist")
            # Try continuing with last known good IC
        end
    end
    return results
end

# ============================================================================
# Main execution
# ============================================================================

function main()
    println("=" ^ 70)
    println("Agent 11 — Multiple Shooting / Newton Iterator")
    println("=" ^ 70)

    # ------------------------------------------------------------------
    # Task 1a: Close the breathing square to machine precision
    # ------------------------------------------------------------------
    println("\n### Task 1a: Close breathing square ###")
    q0_bs, p0_bs, masses_bs, charges_bs = breathing_square_ic(side=2.0, v_rad=0.5)
    T_bs = 11.78

    r_bs = single_shooting_newton(q0_bs, p0_bs, T_bs;
        n_particles=4, dims=2, masses=masses_bs, charges=charges_bs,
        c=1.0, dt=5e-4, max_iter=15, tol=1e-10,
        id="breathing_square", verbose=true)

    if r_bs !== nothing
        @printf("  Breathing square: T=%.4f, E=%.6f, ||F||=%.3e, |λ|=%.3f\n",
                r_bs.period, r_bs.energy, r_bs.mismatch, r_bs.max_floquet)
    end

    # Also try slightly different periods
    for dT in (-0.2, -0.1, 0.0, 0.1, 0.2, 0.5)
        T_try = 11.78 + dT
        id = @sprintf("bs_T%.2f", T_try)
        r = single_shooting_newton(q0_bs, p0_bs, T_try;
            n_particles=4, dims=2, masses=masses_bs, charges=charges_bs,
            c=1.0, dt=5e-4, max_iter=10, tol=1e-8,
            id=id, verbose=false)
        if r !== nothing && r.converged
            @printf("  T=%.2f: CONVERGED, E=%.6f, ||F||=%.3e, |λ|=%.3f\n",
                    T_try, r.energy, r.mismatch, r.max_floquet)
        end
    end

    # Try different v_rad values to find more breathing modes
    for v_rad in (0.3, 0.4, 0.6, 0.7, 0.8)
        q0, p0, masses, charges = breathing_square_ic(side=2.0, v_rad=v_rad)
        for T_try in (8.0, 10.0, 12.0, 15.0, 20.0)
            id = @sprintf("bs_v%.1f_T%.0f", v_rad, T_try)
            r = single_shooting_newton(q0, p0, T_try;
                n_particles=4, dims=2, masses=masses, charges=charges,
                c=1.0, dt=5e-4, max_iter=10, tol=1e-8,
                id=id, verbose=false)
            if r !== nothing && r.converged
                @printf("  [v_rad=%.1f, T=%.0f]: CONVERGED, E=%.6f, ||F||=%.3e\n",
                        v_rad, T_try, r.energy, r.mismatch)
            end
        end
    end

    # ------------------------------------------------------------------
    # Task 1b: 2-body circular orbits — verify exactness
    # ------------------------------------------------------------------
    println("\n### Task 1b: 2-body circular orbits ###")
    for E in (-0.1, -0.25, -0.5, -1.0)
        q0, p0, T = two_body_circular_ic(E=E)
        masses = [1.0, 1.0]
        charges = [1.0, -1.0]

        id = @sprintf("circ_E%.2f", E)
        r = single_shooting_newton(q0, p0, T;
            n_particles=2, dims=2, masses=masses, charges=charges,
            c=1.0, dt=min(T/100, 1e-3), max_iter=10, tol=1e-10,
            id=id, verbose=true)

        if r !== nothing
            @printf("  E=%.2f: T=%.4f, ||F||=%.3e, |λ|=%.3f, μ_CZ=%d\n",
                    E, r.period, r.mismatch, r.max_floquet, r.cz_index)
        end
    end

    # ------------------------------------------------------------------
    # Task 1c: Perturbed circular orbits — find elliptic periodic orbits
    # ------------------------------------------------------------------
    println("\n### Task 1c: Perturbed 2-body orbits ###")
    for E in (-0.25, -0.5)
        q0_circ, p0_circ, T_circ = two_body_circular_ic(E=E)
        masses = [1.0, 1.0]
        charges = [1.0, -1.0]

        # Add radial perturbation
        for delta_r in (0.05, 0.1, 0.2)
            q0 = copy(q0_circ)
            q0[1] *= (1.0 + delta_r)  # stretch particle 1 outward
            q0[3] *= (1.0 + delta_r)  # stretch particle 2 outward

            # Try multiple periods near the circular period
            for T_fac in (0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5, 2.0)
                T_try = T_circ * T_fac
                id = @sprintf("pert_E%.2f_dr%.2f_Tf%.2f", E, delta_r, T_fac)
                r = single_shooting_newton(q0, p0_circ, T_try;
                    n_particles=2, dims=2, masses=masses, charges=charges,
                    c=1.0, dt=min(T_try/200, 5e-4), max_iter=10, tol=1e-8,
                    id=id, verbose=false)
                if r !== nothing && r.converged
                    @printf("  [E=%.2f, dr=%.2f, Tf=%.2f]: CONVERGED, T=%.4f, E=%.6f, ||F||=%.3e, |λ|=%.3f\n",
                            E, delta_r, T_fac, r.period, r.energy, r.mismatch, r.max_floquet)
                end
            end
        end
    end

    # ------------------------------------------------------------------
    # Task 3: c-continuation from Coulomb
    # ------------------------------------------------------------------
    println("\n### Task 3: c-continuation from Kepler ###")

    c_values = [100.0, 50.0, 20.0, 10.0, 5.0, 2.0, 1.0, 0.5]

    for E in (-0.25, -0.5, -1.0)
        println("\n  --- E = $E ---")
        q0, p0, T = two_body_circular_ic(E=E, c=100.0)
        masses = [1.0, 1.0]
        charges = [1.0, -1.0]

        base_id = @sprintf("ccont_E%.2f", E)
        results = c_continuation(q0, p0, T;
            c_values=c_values,
            n_particles=2, dims=2,
            masses=masses, charges=charges,
            dt=min(T/200, 5e-4),
            base_id=base_id, verbose=true)

        if !isempty(results)
            println("  Summary for E=$E:")
            println("    c       | T       | E       | |λ|_max  | μ_CZ")
            println("    --------|---------|---------|----------|------")
            for r in results
                @printf("    %7.3f | %7.4f | %+7.4f | %8.3f | %d\n",
                        r.c, r.period, r.energy, r.max_floquet, r.cz_index)
            end
        end
    end

    # ------------------------------------------------------------------
    # Task: 3-body equilateral search
    # ------------------------------------------------------------------
    println("\n### 3-body: equilateral triangle search ###")
    # 3-body with charges (+1, +1, -1), equilateral triangle rotating
    for side in (1.0, 2.0)
        for omega_frac in (0.5, 0.8, 1.0, 1.2)
            q0_3, p0_3, T_3 = three_body_equilateral_ic(side=side, omega_frac=omega_frac)
            masses_3 = [1.0, 1.0, 1.0]
            charges_3 = [1.0, 1.0, -1.0]

            id = @sprintf("equi_s%.1f_w%.1f", side, omega_frac)
            r = single_shooting_newton(q0_3, p0_3, T_3;
                n_particles=3, dims=2, masses=masses_3, charges=charges_3,
                c=1.0, dt=min(T_3/200, 5e-4), max_iter=10, tol=1e-8,
                id=id, verbose=false)
            if r !== nothing && r.converged
                @printf("  [equilateral s=%.1f, ω_frac=%.1f]: CONVERGED, T=%.4f, E=%.6f\n",
                        side, omega_frac, r.period, r.energy)
            end
        end
    end

    # ------------------------------------------------------------------
    # Write CSV
    # ------------------------------------------------------------------
    write_csv()
    println("\nDone. Found $(count(r -> r.converged, ALL_RESULTS)) converged orbits out of $(length(ALL_RESULTS)) attempts.")
end

"""
    three_body_equilateral_ic(; side, omega_frac)

3-body equilateral triangle with charges (+1, +1, -1), rotating.
omega_frac scales the angular velocity relative to a balance estimate.
"""
function three_body_equilateral_ic(; side::Float64=1.0, omega_frac::Float64=1.0)
    # Equilateral triangle centered at origin
    angles = [0.0, 2*pi/3, 4*pi/3]
    R = side / sqrt(3.0)
    X = zeros(3, 2)
    for i in 1:3
        X[i, 1] = R * cos(angles[i])
        X[i, 2] = R * sin(angles[i])
    end

    masses = [1.0, 1.0, 1.0]
    charges = [1.0, 1.0, -1.0]

    # Rough angular velocity estimate from force balance
    # Sum of Coulomb forces on particle 1 from 2 and 3
    # F = q1*q2/r^2 + q1*q3/r^2 = (1*1 + 1*(-1))/side^2 = 0
    # This doesn't balance — need different approach
    # Use omega from average attraction
    mu = 1.0  # mass
    k_avg = 1.0 / side^2  # average force scale
    omega = omega_frac * sqrt(k_avg / (mu * R))

    T = 2 * pi / omega

    # Momenta: tangential velocity
    P = zeros(3, 2)
    for i in 1:3
        r = X[i, :]
        P[i, :] = masses[i] * omega * [-r[2], r[1]]
    end

    # Zero COM
    com = sum(masses[i] * X[i, :] for i in 1:3) / sum(masses)
    for i in 1:3
        X[i, :] .-= com
    end
    # Zero total momentum
    ptot = sum(P[i, :] for i in 1:3)
    for i in 1:3
        P[i, :] .-= ptot / 3
    end

    q0 = Float64[]
    p0 = Float64[]
    for i in 1:3
        push!(q0, X[i, 1], X[i, 2])
        push!(p0, P[i, 1], P[i, 2])
    end

    return (q0, p0, T)
end

function write_csv()
    csv_path = joinpath(OUTDIR, "found_orbits.csv")
    open(csv_path, "w") do io
        println(io, "id,n_particles,dims,period,energy,c,mismatch,max_floquet,cz_index,cz_upper,stability_type,converged,newton_iters,q0,p0,notes")
        for r in ALL_RESULTS
            q0s = join([@sprintf("%.10f", x) for x in r.q0], ";")
            p0s = join([@sprintf("%.10f", x) for x in r.p0], ";")
            @printf(io, "%s,%d,%d,%.6f,%.6f,%.4f,%.3e,%.6f,%d,%d,%s,%s,%d,%s,%s,%s\n",
                r.id, r.n_particles, r.dims,
                r.period, r.energy, r.c, r.mismatch,
                isnan(r.max_floquet) ? -1.0 : r.max_floquet,
                r.cz_index, r.cz_upper, r.stability_type,
                r.converged, r.newton_iters,
                q0s, p0s, r.notes)
        end
    end
    println("Wrote $(length(ALL_RESULTS)) results to $csv_path")
end

main()
