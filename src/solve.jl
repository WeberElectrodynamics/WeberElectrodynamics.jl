using CommonSolve
using LinearAlgebra: norm, mul!

@inline function strang_splitting_flow!(
    Z_vec::Vector{Float64},
    dt::Float64,
    dq_dt_compiled,
    dp_dt_compiled,
    auxiliary_position_buffer::Vector{Float64},
    momentum_buffer::Vector{Float64},
    position_buffer::Vector{Float64},
    auxiliary_momentum_buffer::Vector{Float64},
    d::Int,
)
    idx_Q_start = 1
    idx_Q_end = d
    idx_X_start = d + 1
    idx_X_end = 2d
    idx_P_start = 2d + 1
    idx_P_end = 3d
    idx_Y_start = 3d + 1
    idx_Y_end = 4d

    @inbounds @views begin
        Q_component = Z_vec[idx_Q_start:idx_Q_end]
        X_component = Z_vec[idx_X_start:idx_X_end]
        P_component = Z_vec[idx_P_start:idx_P_end]
        Y_component = Z_vec[idx_Y_start:idx_Y_end]

        dq_dt_compiled(auxiliary_position_buffer, Q_component, Y_component)
        dp_dt_compiled(momentum_buffer, Q_component, Y_component)

        @. X_component = X_component + auxiliary_position_buffer * (dt / 2)
        @. P_component = P_component + momentum_buffer * (dt / 2)

        dq_dt_compiled(position_buffer, X_component, P_component)
        dp_dt_compiled(auxiliary_momentum_buffer, X_component, P_component)

        @. Q_component = Q_component + position_buffer * dt
        @. Y_component = Y_component + auxiliary_momentum_buffer * dt

        dq_dt_compiled(auxiliary_position_buffer, Q_component, Y_component)
        dp_dt_compiled(momentum_buffer, Q_component, Y_component)

        @. X_component = X_component + auxiliary_position_buffer * (dt / 2)
        @. P_component = P_component + momentum_buffer * (dt / 2)
    end
    return nothing
end

@inline function compute_constraint_residual!(
    result::Vector{Float64},
    μ_val::Vector{Float64},
    A::Matrix{Float64},
    ATμ::Vector{Float64},
    Z::Vector{Float64},
    Ẑ::Vector{Float64},
    Z_result::Vector{Float64},
    dt::Float64,
    dq_dt_compiled,
    dp_dt_compiled,
    auxiliary_position_buffer::Vector{Float64},
    momentum_buffer::Vector{Float64},
    position_buffer::Vector{Float64},
    auxiliary_momentum_buffer::Vector{Float64},
    d::Int,
)
    @inbounds begin
        mul!(ATμ, A', μ_val)
        @. Ẑ = Z + ATμ
        strang_splitting_flow!(
            Ẑ,
            dt,
            dq_dt_compiled,
            dp_dt_compiled,
            auxiliary_position_buffer,
            momentum_buffer,
            position_buffer,
            auxiliary_momentum_buffer,
            d,
        )
        @. Z_result = Ẑ + ATμ
        mul!(result, A, Z_result)
    end
    return nothing
end

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

    # Buffer aliases using paper notation (see docs/theory/SemiExplicitIntegrator.md)
    A = buffers.A
    Z = buffers.Z
    Ẑ = buffers.Ẑ
    Z_result = buffers.Z_result
    position_buffer = buffers.position_buffer
    auxiliary_position_buffer = buffers.auxiliary_position_buffer
    momentum_buffer = buffers.momentum_buffer
    auxiliary_momentum_buffer = buffers.auxiliary_momentum_buffer
    ATμ = buffers.ATμ
    μ = buffers.μ
    μ_prev = buffers.μ_prev
    f_μ = buffers.f_μ

    q = integrator.q
    p = integrator.p

    # Step 1: Embed to extended space: Zₙ = (qₙ, qₙ, pₙ, pₙ)
    @inbounds @views begin
        Z[1:d] .= q
        Z[(d+1):(2d)] .= q
        Z[(2d+1):(3d)] .= p
        Z[(3d+1):(4d)] .= p
    end

    # Step 2: Solve for μ such that f(μ) = 0 using relaxed fixed-point iteration
    relaxation = integrator.alg.relaxation
    for _ = 1:maximum_iterations
        compute_constraint_residual!(
            f_μ,
            μ,
            A,
            ATμ,
            Z,
            Ẑ,
            Z_result,
            dt,
            dq_dt_compiled,
            dp_dt_compiled,
            auxiliary_position_buffer,
            momentum_buffer,
            position_buffer,
            auxiliary_momentum_buffer,
            d,
        )
        copyto!(μ_prev, μ)
        @inbounds @. μ = μ - relaxation * f_μ

        # Check step convergence
        step_norm = norm(μ .- μ_prev)
        if step_norm < convergence_tolerance
            # Verify constraint satisfaction
            compute_constraint_residual!(
                f_μ,
                μ,
                A,
                ATμ,
                Z,
                Ẑ,
                Z_result,
                dt,
                dq_dt_compiled,
                dp_dt_compiled,
                auxiliary_position_buffer,
                momentum_buffer,
                position_buffer,
                auxiliary_momentum_buffer,
                d,
            )
            if norm(f_μ) < convergence_tolerance
                @goto converged
            end
        end
    end

    # Failed to converge
    compute_constraint_residual!(
        f_μ,
        μ,
        A,
        ATμ,
        Z,
        Ẑ,
        Z_result,
        dt,
        dq_dt_compiled,
        dp_dt_compiled,
        auxiliary_position_buffer,
        momentum_buffer,
        position_buffer,
        auxiliary_momentum_buffer,
        d,
    )
    error(
        "Fixed-point iteration failed to converge after $maximum_iterations iterations (residual=$(norm(f_μ)), tolerance=$convergence_tolerance)",
    )

    @label converged

    # Steps 3-5: Compute final state with converged μ
    # Recompute to ensure consistency (solver's last call used previous μ)
    @inbounds begin
        mul!(ATμ, A', μ)
        @. Ẑ = Z + ATμ
        strang_splitting_flow!(
            Ẑ,
            dt,
            dq_dt_compiled,
            dp_dt_compiled,
            auxiliary_position_buffer,
            momentum_buffer,
            position_buffer,
            auxiliary_momentum_buffer,
            d,
        )
    end

    # Step 5: Zₙ₊₁ = Ẑₙ₊₁ + A^T μ (symmetric projection)
    @inbounds @. Ẑ = Ẑ + ATμ

    # Step 6: Extract (qₙ₊₁, pₙ₊₁) from Zₙ₊₁
    @inbounds @views begin
        integrator.q .= Ẑ[1:d]
        integrator.p .= Ẑ[(2d+1):(3d)]
    end

    integrator.step_count += 1
    integrator.t = prob.tspan[1] + integrator.step_count * dt

    # Store in history (in-place copy to pre-allocated vectors)
    idx = integrator.step_count + 1
    @inbounds begin
        integrator.t_history[idx] = integrator.t
        integrator.q_history[idx] .= integrator.q
        integrator.p_history[idx] .= integrator.p
    end

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
