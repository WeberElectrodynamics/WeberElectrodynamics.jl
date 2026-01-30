using CommonSolve
using LinearAlgebra: norm

function CommonSolve.init(
    prob::WeberProblem,
    alg::SymmetricProjectionIntegrator = SymmetricProjectionIntegrator(),
)
    degrees_of_freedom = prob.system.degrees_of_freedom

    # Create workspace buffers (no params needed - baked into compiled functions)
    buffers = SymmetricProjectionBuffers(degrees_of_freedom)

    # Pre-allocate solution storage with all inner vectors
    n_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
    t_history = Vector{Float64}(undef, n_steps + 1)
    q_history = [Vector{Float64}(undef, degrees_of_freedom) for _ = 1:(n_steps+1)]
    p_history = [Vector{Float64}(undef, degrees_of_freedom) for _ = 1:(n_steps+1)]

    # Store initial conditions (in-place copy to pre-allocated vectors)
    t_history[1] = prob.tspan[1]
    q_history[1] .= prob.q_initial
    p_history[1] .= prob.p_initial

    WeberIntegrator(
        prob,
        alg,
        prob.tspan[1],
        prob.tspan[2],
        copy(prob.q_initial),
        copy(prob.p_initial),
        0,
        buffers,
        t_history,
        q_history,
        p_history,
    )
end

function CommonSolve.step!(integrator::WeberIntegrator)
    # Check if already done
    if integrator.t >= integrator.t_end - eps(integrator.t_end)
        return false
    end

    prob = integrator.prob
    dt = prob.dt
    convergence_tolerance = prob.convergence_tolerance
    maximum_iterations = prob.maximum_iterations

    # Compiled functions from WeberSystem (params baked in)
    dq_dt_compiled = prob.system.dq_dt_compiled
    dp_dt_compiled = prob.system.dp_dt_compiled

    buffers = integrator.buffers
    d = buffers.d  # degrees of freedom

    # Extended phase space indices: Z = [Q, X, P, Y] where each component has d elements
    # Q = positions, X = auxiliary positions, P = momenta, Y = auxiliary momenta
    idx_Q_start = 1
    idx_Q_end = d
    idx_X_start = d + 1
    idx_X_end = 2d
    idx_P_start = 2d + 1
    idx_P_end = 3d
    idx_Y_start = 3d + 1
    idx_Y_end = 4d

    # Buffer aliases using paper notation (see docs/theory/SemiExplicitIntegrator.md)
    A = buffers.A              # projection matrix
    Z = buffers.Z              # extended state Zₙ
    Ẑ = buffers.Ẑ              # shifted/evolved state
    Z_result = buffers.Z_result
    position_buffer = buffers.position_buffer
    auxiliary_position_buffer = buffers.auxiliary_position_buffer
    momentum_buffer = buffers.momentum_buffer
    auxiliary_momentum_buffer = buffers.auxiliary_momentum_buffer
    ATμ = buffers.ATμ          # constraint shift A^T μ
    μ = buffers.μ              # Lagrange multipliers
    μ_prev = buffers.μ_prev    # previous iteration
    f_μ = buffers.f_μ          # nonlinear residual f(μ)

    q = integrator.q
    p = integrator.p

    # Step 1: Embed to extended space: Zₙ = (qₙ, qₙ, pₙ, pₙ)
    @views begin
        Z[idx_Q_start:idx_Q_end] .= q
        Z[idx_X_start:idx_X_end] .= q
        Z[idx_P_start:idx_P_end] .= p
        Z[idx_Y_start:idx_Y_end] .= p
    end

    # Strang splitting flow map Φ̂ = Φ^A_{Δt/2} ∘ Φ^B_{Δt} ∘ Φ^A_{Δt/2}
    function Φ̂(Z_vec::Vector{Float64})::Nothing
        @views begin
            Q_component = Z_vec[idx_Q_start:idx_Q_end]
            X_component = Z_vec[idx_X_start:idx_X_end]
            P_component = Z_vec[idx_P_start:idx_P_end]
            Y_component = Z_vec[idx_Y_start:idx_Y_end]

            # Flow A (half step): frozen (Q, Y), evolve (X, P) using H_A(q,y)
            dq_dt_compiled(auxiliary_position_buffer, Q_component, Y_component)
            dp_dt_compiled(momentum_buffer, Q_component, Y_component)

            @. X_component = X_component + auxiliary_position_buffer * (dt / 2)
            @. P_component = P_component + momentum_buffer * (dt / 2)

            # Flow B (full step): frozen (X, P), evolve (Q, Y) using H_B(x,p)
            dq_dt_compiled(position_buffer, X_component, P_component)
            dp_dt_compiled(auxiliary_momentum_buffer, X_component, P_component)

            @. Q_component = Q_component + position_buffer * dt
            @. Y_component = Y_component + auxiliary_momentum_buffer * dt

            # Flow A (half step): frozen (Q, Y), evolve (X, P) using H_A(q,y)
            dq_dt_compiled(auxiliary_position_buffer, Q_component, Y_component)
            dp_dt_compiled(momentum_buffer, Q_component, Y_component)

            @. X_component = X_component + auxiliary_position_buffer * (dt / 2)
            @. P_component = P_component + momentum_buffer * (dt / 2)
        end
        return nothing
    end

    # Nonlinear function f(μ) = A(Φ̂(Zₙ + A^T μ) + A^T μ)
    function f!(result::Vector{Float64}, μ_val::Vector{Float64})::Nothing
        mul!(ATμ, A', μ_val)       # ATμ = A^T μ
        @. Ẑ = Z + ATμ             # Ẑₙ = Zₙ + A^T μ
        Φ̂(Ẑ)                       # Ẑₙ₊₁ = Φ̂(Ẑₙ)
        @. Z_result = Ẑ + ATμ      # Zₙ₊₁ = Ẑₙ₊₁ + A^T μ
        mul!(result, A, Z_result)  # f(μ) = A · Zₙ₊₁
        return nothing
    end

    # Step 2: Solve for μ such that f(μ) = 0 using relaxed fixed-point iteration
    relaxation = integrator.alg.relaxation
    for iter = 1:maximum_iterations
        f!(f_μ, μ)
        copyto!(μ_prev, μ)
        @. μ = μ - relaxation * f_μ

        # Check step convergence
        step_norm = norm(μ .- μ_prev)
        if step_norm < convergence_tolerance
            # Verify constraint satisfaction
            f!(f_μ, μ)
            if norm(f_μ) < convergence_tolerance
                @goto converged
            end
        end
    end

    # Failed to converge
    f!(f_μ, μ)
    error(
        "Fixed-point iteration failed to converge after $maximum_iterations iterations (residual=$(norm(f_μ)), tolerance=$convergence_tolerance)",
    )

    @label converged

    # Steps 3-5: Compute final state with converged μ
    # Recompute to ensure consistency (solver's last f! call used previous μ)
    mul!(ATμ, A', μ)           # ATμ = A^T μ
    @. Ẑ = Z + ATμ             # Ẑₙ = Zₙ + A^T μ
    Φ̂(Ẑ)                       # Ẑₙ₊₁ = Φ̂(Ẑₙ)

    # Step 5: Zₙ₊₁ = Ẑₙ₊₁ + A^T μ (symmetric projection)
    @. Ẑ = Ẑ + ATμ

    # Step 6: Extract (qₙ₊₁, pₙ₊₁) from Zₙ₊₁
    @views begin
        integrator.q .= Ẑ[idx_Q_start:idx_Q_end]
        integrator.p .= Ẑ[idx_P_start:idx_P_end]
    end

    integrator.step_count += 1
    integrator.t = prob.tspan[1] + integrator.step_count * dt

    # Store in history (in-place copy to pre-allocated vectors)
    idx = integrator.step_count + 1
    integrator.t_history[idx] = integrator.t
    integrator.q_history[idx] .= integrator.q
    integrator.p_history[idx] .= integrator.p

    return integrator.t < integrator.t_end
end

function CommonSolve.solve!(integrator::WeberIntegrator)
    retcode = :Success
    try
        while CommonSolve.step!(integrator)
            # Continue stepping
        end
    catch e
        if e isa ErrorException && contains(e.msg, "Fixed-point iteration failed")
            retcode = :Failure
        else
            rethrow()
        end
    end

    # Trim history to actual steps taken
    n = integrator.step_count + 1
    WeberSolution(
        integrator.t_history[1:n],
        integrator.q_history[1:n],
        integrator.p_history[1:n],
        integrator.prob,
        retcode,
    )
end

function CommonSolve.solve(
    prob::WeberProblem,
    alg::WeberAlgorithm = SymmetricProjectionIntegrator(),
)
    integrator = CommonSolve.init(prob, alg)
    CommonSolve.solve!(integrator)
end
