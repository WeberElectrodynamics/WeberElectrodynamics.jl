using CommonSolve
using LinearAlgebra: norm, mul!, Transpose

@inline function strang_splitting_flow!(
    Z_vec::Vector{Float64},
    dt::Float64,
    dq_dt_compiled,
    dp_dt_compiled,
    params::Vector{Float64},
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

        dq_dt_compiled(auxiliary_position_buffer, Q_component, Y_component, params)
        dp_dt_compiled(momentum_buffer, Q_component, Y_component, params)

        @. X_component = X_component + auxiliary_position_buffer * (dt / 2)
        @. P_component = P_component + momentum_buffer * (dt / 2)

        dq_dt_compiled(position_buffer, X_component, P_component, params)
        dp_dt_compiled(auxiliary_momentum_buffer, X_component, P_component, params)

        @. Q_component = Q_component + position_buffer * dt
        @. Y_component = Y_component + auxiliary_momentum_buffer * dt

        dq_dt_compiled(auxiliary_position_buffer, Q_component, Y_component, params)
        dp_dt_compiled(momentum_buffer, Q_component, Y_component, params)

        @. X_component = X_component + auxiliary_position_buffer * (dt / 2)
        @. P_component = P_component + momentum_buffer * (dt / 2)
    end
    return nothing
end

@inline function compute_constraint_residual!(
    result::Vector{Float64},
    μ_val::Vector{Float64},
    A::Matrix{Float64},
    A_transpose::Transpose{Float64,Matrix{Float64}},
    ATμ::Vector{Float64},
    Z::Vector{Float64},
    Ẑ::Vector{Float64},
    Z_result::Vector{Float64},
    dt::Float64,
    dq_dt_compiled,
    dp_dt_compiled,
    params::Vector{Float64},
    auxiliary_position_buffer::Vector{Float64},
    momentum_buffer::Vector{Float64},
    position_buffer::Vector{Float64},
    auxiliary_momentum_buffer::Vector{Float64},
    d::Int,
)
    @inbounds begin
        mul!(ATμ, A_transpose, μ_val)
        @. Ẑ = Z + ATμ
        strang_splitting_flow!(
            Ẑ,
            dt,
            dq_dt_compiled,
            dp_dt_compiled,
            params,
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

@inline function _projected_cartesian_step!(
    q::Vector{Float64},
    p::Vector{Float64},
    dt_step::Float64,
    prob::WeberProblem,
    alg::SymmetricProjectionIntegrator,
    buffers::SymmetricProjectionBuffers,
)
    convergence_tolerance = prob.convergence_tolerance
    maximum_iterations = prob.maximum_iterations

    dq_dt_compiled = prob.system.dq_dt_compiled
    dp_dt_compiled = prob.system.dp_dt_compiled
    params = prob.params

    d = buffers.d

    A = buffers.A
    A_transpose = buffers.A_transpose
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
    diff_buffer = buffers.diff_buffer

    @inbounds @views begin
        Z[1:d] .= q
        Z[(d+1):(2d)] .= q
        Z[(2d+1):(3d)] .= p
        Z[(3d+1):(4d)] .= p
    end

    relaxation = alg.relaxation
    for _ = 1:maximum_iterations
        compute_constraint_residual!(
            f_μ,
            μ,
            A,
            A_transpose,
            ATμ,
            Z,
            Ẑ,
            Z_result,
            dt_step,
            dq_dt_compiled,
            dp_dt_compiled,
            params,
            auxiliary_position_buffer,
            momentum_buffer,
            position_buffer,
            auxiliary_momentum_buffer,
            d,
        )
        copyto!(μ_prev, μ)
        @inbounds @. μ = μ - relaxation * f_μ

        @inbounds @. diff_buffer = μ - μ_prev
        step_norm = norm(diff_buffer)
        if step_norm < convergence_tolerance
            compute_constraint_residual!(
                f_μ,
                μ,
                A,
                A_transpose,
                ATμ,
                Z,
                Ẑ,
                Z_result,
                dt_step,
                dq_dt_compiled,
                dp_dt_compiled,
                params,
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

    compute_constraint_residual!(
        f_μ,
        μ,
        A,
        A_transpose,
        ATμ,
        Z,
        Ẑ,
        Z_result,
        dt_step,
        dq_dt_compiled,
        dp_dt_compiled,
        params,
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

    @inbounds begin
        mul!(ATμ, A_transpose, μ)
        @. Ẑ = Z + ATμ
        strang_splitting_flow!(
            Ẑ,
            dt_step,
            dq_dt_compiled,
            dp_dt_compiled,
            params,
            auxiliary_position_buffer,
            momentum_buffer,
            position_buffer,
            auxiliary_momentum_buffer,
            d,
        )
    end

    @inbounds @. Ẑ = Ẑ + ATμ

    @inbounds @views begin
        q .= Ẑ[1:d]
        p .= Ẑ[(2d+1):(3d)]
    end

    return nothing
end

@inline function _update_diagnostics!(
    integrator::WeberIntegrator,
    mode::UInt8,
    substeps::Int,
    min_distance::Float64,
    max_constraint_violation::Float64,
)
    diagnostics = integrator.diagnostics
    idx = integrator.step_count
    if idx <= length(diagnostics.mode_history)
        diagnostics.mode_history[idx] = mode
    end

    diagnostics.total_substeps += substeps
    if substeps > diagnostics.max_substeps_used
        diagnostics.max_substeps_used = substeps
    end

    if isfinite(min_distance)
        diagnostics.min_encounter_distance = min(diagnostics.min_encounter_distance, min_distance)
    end
    if max_constraint_violation > diagnostics.max_constraint_violation
        diagnostics.max_constraint_violation = max_constraint_violation
    end

    if mode == REG_MODE_PAIR
        diagnostics.active_steps += 1
        diagnostics.pair_steps += 1
    elseif mode == REG_MODE_CHAIN
        diagnostics.active_steps += 1
        diagnostics.chain_steps += 1
    else
        diagnostics.unregularized_steps += 1
    end

    return nothing
end

@inline function _store_state!(integrator::WeberIntegrator, dt_step::Float64)
    integrator.step_count += 1
    integrator.t += dt_step
    if integrator.t > integrator.t_end
        integrator.t = integrator.t_end
    end

    idx = integrator.step_count + 1
    @inbounds begin
        integrator.t_history[idx] = integrator.t
        integrator.q_history[idx] .= integrator.q
        integrator.p_history[idx] .= integrator.p
    end

    return nothing
end

@inline function _step_unregularized!(
    integrator::WeberIntegrator,
    dt_step::Float64,
)
    _projected_cartesian_step!(
        integrator.q,
        integrator.p,
        dt_step,
        integrator.prob,
        integrator.alg,
        integrator.buffers,
    )

    _store_state!(integrator, dt_step)
    _update_diagnostics!(integrator, REG_MODE_NONE, 1, Inf, 0.0)
    return nothing
end

@inline function _step_regularized_pair!(
    integrator::WeberIntegrator,
    dt_step::Float64,
    min_distance::Float64,
)
    prob = integrator.prob
    rb = integrator.buffers.regularization_buffers
    dims = rb.dims

    i = rb.active_anchor_i
    j = rb.active_anchor_j

    reg = prob.regularization
    substeps = Int(ceil(rb.r_on / max(min_distance, reg.g_floor)))
    substeps = clamp(substeps, 1, reg.max_substeps)

    dt_sub = dt_step / substeps
    max_constraint = 0.0

    @inbounds for _ = 1:substeps
        _extract_pair_relative_state!(rb, integrator.q, integrator.p, prob.masses, i, j)

        if dims == 1
            x = rb.rel_q[1]
            px = rb.rel_p[1]
            s = x >= 0 ? 1.0 : -1.0
            u = sqrt(abs(x))
            rb.lc_u[1] = u
            rb.lc_U[1] = 2 * s * u * px
        elseif dims == 2
            _lc_lift!(rb)
            _lc_project!(rb.rel_q, rb.rel_p, rb.lc_u, rb.lc_U)
        elseif dims == 3
            _ks_lift!(rb)
            c_err = _ks_project_constraint!(rb.ks_U, rb.ks_u, rb.ks_n)
            if c_err > max_constraint
                max_constraint = c_err
            end
        end

        _projected_cartesian_step!(
            integrator.q,
            integrator.p,
            dt_sub,
            prob,
            integrator.alg,
            integrator.buffers,
        )
    end

    _store_state!(integrator, dt_step)
    _update_diagnostics!(integrator, REG_MODE_PAIR, substeps, min_distance, max_constraint)
    return nothing
end

@inline function _step_regularized_chain!(
    integrator::WeberIntegrator,
    dt_step::Float64,
    min_distance::Float64,
)
    prob = integrator.prob
    rb = integrator.buffers.regularization_buffers
    reg = prob.regularization

    omega = _component_omega(rb)
    g = max(1.0 / max(omega, eps(Float64)), reg.g_floor)

    substeps = Int(ceil(dt_step / g))
    substeps = clamp(substeps, 1, reg.max_substeps)

    dt_sub = dt_step / substeps
    max_constraint = 0.0

    dims = rb.dims
    chain_count = rb.active_count

    @inbounds for _ = 1:substeps
        if dims == 2
            for idx = 1:(chain_count-1)
                i = rb.chain_order[idx]
                j = rb.chain_order[idx+1]
                _extract_pair_relative_state!(rb, integrator.q, integrator.p, prob.masses, i, j)
                _lc_lift!(rb)
            end
        elseif dims == 3
            for idx = 1:(chain_count-1)
                i = rb.chain_order[idx]
                j = rb.chain_order[idx+1]
                _extract_pair_relative_state!(rb, integrator.q, integrator.p, prob.masses, i, j)
                _ks_lift!(rb)
                c_err = _ks_project_constraint!(rb.ks_U, rb.ks_u, rb.ks_n)
                if c_err > max_constraint
                    max_constraint = c_err
                end
            end
        end

        _projected_cartesian_step!(
            integrator.q,
            integrator.p,
            dt_sub,
            prob,
            integrator.alg,
            integrator.buffers,
        )
    end

    _store_state!(integrator, dt_step)
    _update_diagnostics!(integrator, REG_MODE_CHAIN, substeps, min_distance, max_constraint)
    return nothing
end

@inline function _step_regularized_dispatch!(
    integrator::WeberIntegrator,
    dt_step::Float64,
)
    prob = integrator.prob
    reg = prob.regularization
    rb = integrator.buffers.regularization_buffers
    diagnostics = integrator.diagnostics

    was_active = rb.is_active
    active, mode, min_distance =
        _detect_regularization_component!(rb, integrator.q, reg.chain_enabled)

    if !was_active && active
        diagnostics.activation_count += 1
    elseif was_active && !active
        diagnostics.deactivation_count += 1
    end

    if !active
        _step_unregularized!(integrator, dt_step)
        return nothing
    end

    if mode == REG_MODE_CHAIN
        _step_regularized_chain!(integrator, dt_step, min_distance)
    else
        _step_regularized_pair!(integrator, dt_step, min_distance)
    end

    return nothing
end

@inline function _clone_diagnostics(
    diagnostics::RegularizationDiagnostics,
    n_steps::Int,
)
    out = RegularizationDiagnostics(diagnostics.enabled, n_steps)
    out.activation_count = diagnostics.activation_count
    out.deactivation_count = diagnostics.deactivation_count
    out.active_steps = diagnostics.active_steps
    out.pair_steps = diagnostics.pair_steps
    out.chain_steps = diagnostics.chain_steps
    out.unregularized_steps = diagnostics.unregularized_steps
    out.total_substeps = diagnostics.total_substeps
    out.max_substeps_used = diagnostics.max_substeps_used
    out.max_constraint_violation = diagnostics.max_constraint_violation
    out.min_encounter_distance = diagnostics.min_encounter_distance

    if n_steps > 0
        @inbounds out.mode_history[1:n_steps] .= diagnostics.mode_history[1:n_steps]
    end

    return out
end

function CommonSolve.init(
    prob::WeberProblem,
    alg::SymmetricProjectionIntegrator = SymmetricProjectionIntegrator(),
)
    degrees_of_freedom = prob.system.degrees_of_freedom

    buffers = SymmetricProjectionBuffers(prob)

    n_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
    t_history = Vector{Float64}(undef, n_steps + 1)
    q_history = [Vector{Float64}(undef, degrees_of_freedom) for _ = 1:(n_steps+1)]
    p_history = [Vector{Float64}(undef, degrees_of_freedom) for _ = 1:(n_steps+1)]

    t_history[1] = prob.tspan[1]
    q_history[1] .= prob.q_initial
    p_history[1] .= prob.p_initial

    diagnostics = RegularizationDiagnostics(prob.regularization.enabled, n_steps)

    WeberIntegrator(
        prob,
        alg,
        prob.tspan[1],
        prob.tspan[2],
        copy(prob.q_initial),
        copy(prob.p_initial),
        0,
        buffers,
        diagnostics,
        t_history,
        q_history,
        p_history,
    )
end

function CommonSolve.step!(integrator::WeberIntegrator)
    max_steps = length(integrator.t_history) - 1

    if integrator.step_count >= max_steps || integrator.t >= integrator.t_end - eps(integrator.t_end)
        integrator.t = integrator.t_end
        return false
    end

    prob = integrator.prob
    dt = prob.dt
    is_final_step = integrator.step_count == max_steps - 1
    dt_step = is_final_step ? (integrator.t_end - integrator.t) : dt
    if dt_step <= 0
        integrator.t = integrator.t_end
        return false
    end

    if prob.regularization.enabled
        _step_regularized_dispatch!(integrator, dt_step)
    else
        _step_unregularized!(integrator, dt_step)
    end

    return integrator.step_count < max_steps
end

function CommonSolve.solve!(integrator::WeberIntegrator)
    retcode = :Success
    try
        while CommonSolve.step!(integrator)
        end
    catch e
        if e isa ErrorException && contains(e.msg, "Fixed-point iteration failed")
            retcode = :Failure
        else
            rethrow()
        end
    end

    n = integrator.step_count + 1
    diagnostics = _clone_diagnostics(integrator.diagnostics, integrator.step_count)

    WeberSolution(
        integrator.t_history[1:n],
        integrator.q_history[1:n],
        integrator.p_history[1:n],
        integrator.prob,
        retcode,
        diagnostics,
    )
end

function CommonSolve.solve(
    prob::WeberProblem,
    alg::WeberAlgorithm = SymmetricProjectionIntegrator(),
)
    integrator = CommonSolve.init(prob, alg)
    CommonSolve.solve!(integrator)
end
