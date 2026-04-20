# Agent 11 — Direct action minimization for periodic orbits
#             of the Weber Lagrangian.
#
# The Weber Lagrangian: L = T - V_Weber
#   T = Σ_i |qdot_i|²/(2m_i)
#   V_Weber = Σ_{i<j} q_i q_j / r_{ij} * (1 + rdot_{ij}²/(2c²))
#
# Discretize a loop q(t) on N time points, minimize the action
#   A = Σ_k L(q_k, (q_{k+1}-q_k)/Δt) * Δt
# over {q_k} using gradient descent.
#
# Run:
#   julia --project=. research/homology/11_variational_shooting/action_minimizer.jl
#
# Outputs:
#   - Results appended to found_orbits.csv (if multiple_shooting.jl ran first)
#   - Convergence info to stdout

using LinearAlgebra
using Printf

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using CommonSolve: solve

const OUTDIR = @__DIR__

# ============================================================================
# Discretized Weber action
# ============================================================================

"""
    weber_lagrangian(q_k, q_kp1, dt, masses, charges, c, n_particles, dims)

Evaluate the Weber Lagrangian at the midpoint of a time step.
q_k, q_kp1 are flat position vectors at times t_k and t_{k+1}.
"""
function weber_lagrangian(q_k::Vector{Float64}, q_kp1::Vector{Float64},
                          dt::Float64, masses::Vector{Float64},
                          charges::Vector{Float64}, c::Float64,
                          n_particles::Int, dims::Int)
    dof = n_particles * dims
    # Velocities
    qdot = (q_kp1 - q_k) / dt

    # Kinetic energy
    T = 0.0
    for i in 1:n_particles
        idx = (i-1)*dims
        v2 = 0.0
        for d in 1:dims
            v2 += qdot[idx+d]^2
        end
        T += 0.5 * masses[i] * v2
    end

    # Weber potential: V = Σ_{i<j} q_i q_j / r * (1 + rdot²/(2c²))
    # Use midpoint for positions
    q_mid = 0.5 * (q_k + q_kp1)
    V = 0.0
    for i in 1:n_particles
        idx_i = (i-1)*dims
        for j in (i+1):n_particles
            idx_j = (j-1)*dims
            # Relative position at midpoint
            r2 = 0.0
            for d in 1:dims
                dr = q_mid[idx_i+d] - q_mid[idx_j+d]
                r2 += dr^2
            end
            r = sqrt(r2)
            r < 1e-15 && continue

            # Radial velocity: rdot = (r_vec . v_vec) / r
            rdot_ij = 0.0
            for d in 1:dims
                dr = q_mid[idx_i+d] - q_mid[idx_j+d]
                dv = qdot[idx_i+d] - qdot[idx_j+d]
                rdot_ij += dr * dv
            end
            rdot_ij /= r

            V += charges[i] * charges[j] / r * (1.0 + rdot_ij^2 / (2.0 * c^2))
        end
    end

    return T - V
end

"""
    compute_action(Q, T, masses, charges, c, n_particles, dims)

Compute the discretized action for a loop Q = [q_1, q_2, ..., q_N]
with period T. The loop is closed: q_{N+1} = q_1.
"""
function compute_action(Q::Vector{Vector{Float64}}, T::Float64,
                        masses::Vector{Float64}, charges::Vector{Float64},
                        c::Float64, n_particles::Int, dims::Int)
    N = length(Q)
    dt = T / N
    action = 0.0
    for k in 1:N
        kp1 = k == N ? 1 : k + 1
        action += weber_lagrangian(Q[k], Q[kp1], dt, masses, charges, c,
                                    n_particles, dims) * dt
    end
    return action
end

"""
    action_gradient(Q, T, masses, charges, c, n_particles, dims; eps=1e-7)

Compute the gradient of the action with respect to Q (all node positions)
by finite differences.
"""
function action_gradient(Q::Vector{Vector{Float64}}, T::Float64,
                         masses::Vector{Float64}, charges::Vector{Float64},
                         c::Float64, n_particles::Int, dims::Int;
                         eps::Float64=1e-7)
    N = length(Q)
    dof = n_particles * dims
    A0 = compute_action(Q, T, masses, charges, c, n_particles, dims)

    grad = [zeros(dof) for _ in 1:N]
    for k in 1:N
        for d in 1:dof
            Q_plus = deepcopy(Q)
            Q_plus[k][d] += eps
            A_plus = compute_action(Q_plus, T, masses, charges, c, n_particles, dims)
            grad[k][d] = (A_plus - A0) / eps
        end
    end
    return grad
end

"""
    action_gradient_analytical(Q, T, masses, charges, c, n_particles, dims)

Compute the action gradient analytically (Euler-Lagrange discrete equations).
The gradient dA/dq_k involves the Lagrangian terms from steps k-1 and k.
"""
function action_gradient_analytical(Q::Vector{Vector{Float64}}, T::Float64,
                                    masses::Vector{Float64}, charges::Vector{Float64},
                                    c::Float64, n_particles::Int, dims::Int)
    N = length(Q)
    dof = n_particles * dims
    dt = T / N
    eps_fd = 1e-8

    grad = [zeros(dof) for _ in 1:N]

    for k in 1:N
        km1 = k == 1 ? N : k - 1
        kp1 = k == N ? 1 : k + 1

        for d in 1:dof
            # Perturb q_k[d] and measure change in L_{k-1} and L_k
            # L_{k-1} = L(q_{k-1}, q_k, dt) — q_k appears as "next"
            # L_k     = L(q_k, q_{k+1}, dt) — q_k appears as "current"

            Q_p = copy(Q[k]); Q_p[d] += eps_fd
            Q_m = copy(Q[k]); Q_m[d] -= eps_fd

            # dL_{k-1}/dq_k
            L_km1_p = weber_lagrangian(Q[km1], Q_p, dt, masses, charges, c, n_particles, dims)
            L_km1_m = weber_lagrangian(Q[km1], Q_m, dt, masses, charges, c, n_particles, dims)

            # dL_k/dq_k
            L_k_p = weber_lagrangian(Q_p, Q[kp1], dt, masses, charges, c, n_particles, dims)
            L_k_m = weber_lagrangian(Q_m, Q[kp1], dt, masses, charges, c, n_particles, dims)

            grad[k][d] = ((L_km1_p - L_km1_m) + (L_k_p - L_k_m)) / (2 * eps_fd) * dt
        end
    end

    return grad
end

# ============================================================================
# Action minimization with gradient descent
# ============================================================================

"""
    minimize_action(Q_init, T; kwargs...) -> (Q_final, action, converged)

Minimize the discretized Weber action over a loop by gradient descent
with momentum (heavy ball method).
"""
function minimize_action(Q_init::Vector{Vector{Float64}}, T::Float64;
                         masses::Vector{Float64}, charges::Vector{Float64},
                         c::Float64, n_particles::Int, dims::Int,
                         max_iter::Int=200, lr::Float64=1e-4,
                         tol::Float64=1e-8, momentum::Float64=0.9,
                         verbose::Bool=true, id::String="action")
    N = length(Q_init)
    dof = n_particles * dims
    Q = deepcopy(Q_init)
    velocity = [zeros(dof) for _ in 1:N]

    A_prev = compute_action(Q, T, masses, charges, c, n_particles, dims)
    verbose && @printf("  [%s] iter 0: action = %.8f\n", id, A_prev)

    best_Q = deepcopy(Q)
    best_A = A_prev

    for iter in 1:max_iter
        grad = action_gradient_analytical(Q, T, masses, charges, c, n_particles, dims)

        # Gradient norm
        gnorm = sqrt(sum(sum(g.^2) for g in grad))

        # Heavy ball update
        for k in 1:N
            velocity[k] = momentum * velocity[k] - lr * grad[k]
            Q[k] += velocity[k]
        end

        # Enforce COM = 0 at each node (to remove translation drift)
        for k in 1:N
            com = zeros(dims)
            M_total = sum(masses)
            for i in 1:n_particles
                idx = (i-1)*dims
                for d in 1:dims
                    com[d] += masses[i] * Q[k][idx+d]
                end
            end
            com ./= M_total
            for i in 1:n_particles
                idx = (i-1)*dims
                for d in 1:dims
                    Q[k][idx+d] -= com[d]
                end
            end
        end

        A = compute_action(Q, T, masses, charges, c, n_particles, dims)

        if A < best_A
            best_A = A
            best_Q = deepcopy(Q)
        end

        if iter % 20 == 0 || iter <= 5
            verbose && @printf("  [%s] iter %d: action = %.8f, |grad| = %.3e\n",
                               id, iter, A, gnorm)
        end

        if gnorm < tol
            verbose && @printf("  [%s] CONVERGED at iter %d: |grad| = %.3e\n", id, iter, gnorm)
            return (Q, A, true, iter)
        end

        # Adaptive learning rate
        if A > A_prev + 0.1 * abs(A_prev)
            lr *= 0.5
            velocity = [zeros(dof) for _ in 1:N]
            Q = deepcopy(best_Q)
            verbose && @printf("  [%s] iter %d: action increased, reducing lr to %.2e\n",
                               id, iter, lr)
        end

        A_prev = A
    end

    verbose && println("  [$id] did not converge after $max_iter iterations")
    return (best_Q, best_A, false, max_iter)
end

# ============================================================================
# Seed generators
# ============================================================================

"""
    circular_orbit_loop(; E, n_nodes, n_particles=2, dims=2, kwargs...)

Generate a discretized circular orbit loop for a 2-body unlike-charge system.
"""
function circular_orbit_loop(; E::Float64=-0.5, n_nodes::Int=64,
                               m1::Float64=1.0, m2::Float64=1.0,
                               q1::Float64=1.0, q2::Float64=-1.0,
                               c::Float64=1.0)
    mu = m1 * m2 / (m1 + m2)
    k = -q1 * q2
    a = k / (2 * abs(E))
    v = sqrt(k / (mu * a))
    omega = v / a
    T = 2 * pi / omega

    M = m1 + m2
    r1 = m2 / M * a  # distance from COM for particle 1
    r2 = m1 / M * a

    Q = Vector{Float64}[]
    for i in 0:(n_nodes-1)
        theta = 2 * pi * i / n_nodes
        # Particle 1 (charge q1) orbits at radius r1
        x1 = r1 * cos(theta)
        y1 = r1 * sin(theta)
        # Particle 2 (charge q2) opposite
        x2 = -r2 * cos(theta)
        y2 = -r2 * sin(theta)
        push!(Q, Float64[x1, y1, x2, y2])
    end

    return (Q, T)
end

"""
    breathing_square_loop(; n_nodes, side, v_rad, T_approx)

Generate a discretized breathing square loop by integrating forward.
"""
function breathing_square_loop(; n_nodes::Int=64, side::Float64=2.0,
                                 v_rad::Float64=0.5, T_approx::Float64=11.78,
                                 c::Float64=1.0)
    s = side / 2
    positions = [s s; -s -s; -s s; s -s]
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

    # Integrate to get the trajectory
    sys = HamiltonianSystem(4, 2)
    dt_int = T_approx / (n_nodes * 4)
    prob = HamiltonianProblem(sys, (0.0, T_approx), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt_int)
    sol = solve(prob, SymmetricProjectionIntegrator())

    if sol.retcode != :Success
        @warn "Failed to integrate breathing square for loop seed"
        return nothing
    end

    # Sample n_nodes points from the solution
    stride = max(1, length(sol.t) ÷ n_nodes)
    Q = Vector{Float64}[]
    for i in 1:n_nodes
        idx = min(1 + (i-1) * stride, length(sol.q))
        push!(Q, copy(sol.q[idx]))
    end

    return (Q, T_approx, masses, charges)
end

# ============================================================================
# Main execution
# ============================================================================

function main()
    println("=" ^ 70)
    println("Agent 11 — Direct Action Minimization")
    println("=" ^ 70)

    results_file = joinpath(OUTDIR, "action_results.csv")
    results = []

    # ------------------------------------------------------------------
    # Task 2a: 2-body circular orbit — verify action is stationary
    # ------------------------------------------------------------------
    println("\n### Task 2a: 2-body circular orbit — action stationarity ###")

    for E in (-0.25, -0.5, -1.0)
        masses = [1.0, 1.0]
        charges = [1.0, -1.0]

        Q_circ, T_circ = circular_orbit_loop(E=E, n_nodes=64)
        A_circ = compute_action(Q_circ, T_circ, masses, charges, 1.0, 2, 2)
        grad = action_gradient_analytical(Q_circ, T_circ, masses, charges, 1.0, 2, 2)
        gnorm = sqrt(sum(sum(g.^2) for g in grad))

        @printf("  E=%.2f: action=%.6f, |grad|=%.3e, T=%.4f\n",
                E, A_circ, gnorm, T_circ)

        push!(results, (id=@sprintf("circ_E%.2f", E),
                        action=A_circ, grad_norm=gnorm, T=T_circ, E=E,
                        converged=gnorm < 1e-4, method="exact_seed"))
    end

    # ------------------------------------------------------------------
    # Task 2a continued: Perturbed circular, minimize back
    # ------------------------------------------------------------------
    println("\n### Task 2a: Perturbed circular -> minimize ###")

    for E in (-0.25, -0.5)
        masses = [1.0, 1.0]
        charges = [1.0, -1.0]

        Q_circ, T_circ = circular_orbit_loop(E=E, n_nodes=32)

        # Perturb
        Q_pert = deepcopy(Q_circ)
        for k in 1:length(Q_pert)
            for d in 1:4
                Q_pert[k][d] += 0.05 * randn() * abs(Q_pert[k][d] + 0.01)
            end
        end

        A_before = compute_action(Q_pert, T_circ, masses, charges, 1.0, 2, 2)
        @printf("  E=%.2f: action before=%.6f\n", E, A_before)

        Q_opt, A_opt, conv, iters = minimize_action(Q_pert, T_circ;
            masses=masses, charges=charges, c=1.0,
            n_particles=2, dims=2,
            max_iter=100, lr=1e-5, tol=1e-6,
            verbose=true, id=@sprintf("pert_E%.2f", E))

        @printf("  E=%.2f: action after=%.6f, converged=%s, iters=%d\n",
                E, A_opt, conv, iters)

        push!(results, (id=@sprintf("pert_opt_E%.2f", E),
                        action=A_opt, grad_norm=NaN, T=T_circ, E=E,
                        converged=conv, method="perturbed_minimized"))
    end

    # ------------------------------------------------------------------
    # Task 2c: 4-body breathing square action
    # ------------------------------------------------------------------
    println("\n### Task 2c: 4-body breathing square action ###")

    loop_data = breathing_square_loop(n_nodes=32, T_approx=11.78)
    if loop_data !== nothing
        Q_bs, T_bs, masses_bs, charges_bs = loop_data
        A_bs = compute_action(Q_bs, T_bs, masses_bs, charges_bs, 1.0, 4, 2)
        grad_bs = action_gradient_analytical(Q_bs, T_bs, masses_bs, charges_bs, 1.0, 4, 2)
        gnorm_bs = sqrt(sum(sum(g.^2) for g in grad_bs))

        @printf("  Breathing square: action=%.6f, |grad|=%.3e\n", A_bs, gnorm_bs)

        # Try to minimize
        Q_opt, A_opt, conv, iters = minimize_action(Q_bs, T_bs;
            masses=masses_bs, charges=charges_bs, c=1.0,
            n_particles=4, dims=2,
            max_iter=100, lr=1e-6, tol=1e-6,
            verbose=true, id="bs_action")

        @printf("  After optimization: action=%.6f, converged=%s, iters=%d\n",
                A_opt, conv, iters)

        push!(results, (id="bs_action_opt", action=A_opt, grad_norm=gnorm_bs,
                        T=T_bs, E=NaN, converged=conv, method="bs_minimized"))
    end

    # ------------------------------------------------------------------
    # Task 2: Try different periods for 2-body action to find new orbits
    # ------------------------------------------------------------------
    println("\n### Task 2: Period scan for 2-body action minimization ###")
    masses = [1.0, 1.0]
    charges = [1.0, -1.0]

    for E_seed in (-0.25, -0.5)
        Q_seed, T_seed = circular_orbit_loop(E=E_seed, n_nodes=32)

        for T_fac in (0.5, 1.5, 2.0, 3.0)
            T_try = T_seed * T_fac
            # Rescale the loop in time (stretch/compress the positions)
            Q_try = deepcopy(Q_seed)
            if T_fac != 1.0
                # For period-multiplied orbits, repeat the loop
                if T_fac == 2.0
                    Q_try = vcat(Q_seed, Q_seed)
                elseif T_fac == 3.0
                    Q_try = vcat(Q_seed, Q_seed, Q_seed)
                end
            end

            id = @sprintf("Tscan_E%.2f_Tf%.1f", E_seed, T_fac)
            Q_opt, A_opt, conv, iters = minimize_action(Q_try, T_try;
                masses=masses, charges=charges, c=1.0,
                n_particles=2, dims=2,
                max_iter=80, lr=1e-5, tol=1e-6,
                verbose=false, id=id)

            grad_f = action_gradient_analytical(Q_opt, T_try, masses, charges, 1.0, 2, 2)
            gnorm_f = sqrt(sum(sum(g.^2) for g in grad_f))
            @printf("  E_seed=%.2f, T_fac=%.1f: action=%.6f, |grad|=%.3e, conv=%s\n",
                    E_seed, T_fac, A_opt, gnorm_f, conv)

            push!(results, (id=id, action=A_opt, grad_norm=gnorm_f,
                            T=T_try, E=E_seed, converged=conv, method="period_scan"))
        end
    end

    # ------------------------------------------------------------------
    # Write results
    # ------------------------------------------------------------------
    open(results_file, "w") do io
        println(io, "id,action,grad_norm,period,energy_seed,converged,method")
        for r in results
            @printf(io, "%s,%.8f,%.3e,%.6f,%.4f,%s,%s\n",
                    r.id, r.action,
                    isnan(r.grad_norm) ? -1.0 : r.grad_norm,
                    r.T, r.E, r.converged, r.method)
        end
    end
    println("\nWrote $(length(results)) results to $results_file")
end

main()
