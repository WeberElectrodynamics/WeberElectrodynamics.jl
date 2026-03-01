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

@inline function _merge_used_backend!(
    diagnostics::RegularizationDiagnostics,
    backend::Symbol,
)
    if !diagnostics.enabled || backend == REG_BACKEND_DISABLED
        return nothing
    end

    if diagnostics.used_backend == REG_BACKEND_DISABLED
        diagnostics.used_backend = backend
    elseif diagnostics.used_backend != backend
        diagnostics.used_backend = REG_BACKEND_MIXED
    end

    return nothing
end

@inline function _update_diagnostics!(
    integrator::WeberIntegrator,
    mode::UInt8,
    substeps::Int,
    min_distance::Float64,
    max_constraint_violation::Float64,
    pair_backend::Symbol,
    backend_fallback_step::Bool,
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
        if pair_backend == REG_BACKEND_LIFTED
            diagnostics.lifted_pair_steps += 1
        else
            diagnostics.adaptive_pair_steps += 1
        end
        _merge_used_backend!(diagnostics, pair_backend)
        if backend_fallback_step
            diagnostics.backend_fallback_steps += 1
        end
    elseif mode == REG_MODE_CHAIN
        diagnostics.active_steps += 1
        diagnostics.chain_steps += 1
        _merge_used_backend!(diagnostics, REG_BACKEND_ADAPTIVE)
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

@inline function _set_pair_params!(
    rb::RegularizationBuffers,
    prob::WeberProblem,
    i::Int,
    j::Int,
)
    n = rb.n_particles
    params_pair = rb.params_pair
    masses = prob.masses
    charges = prob.charges
    kappas = prob.kappas

    @inbounds begin
        for k = 1:n
            params_pair[k] = masses[k]
            params_pair[n+k] = 0.0
        end
        params_pair[n+i] = charges[i]
        params_pair[n+j] = charges[j]
        params_pair[2n+1] = prob.c
        # Copy kappas; charges for non-(i,j) pairs are zero so their κ values
        # do not affect the force, but the buffer must have the right length.
        n_pairs = rb.n_pairs
        for k = 1:n_pairs
            params_pair[2n+1+k] = kappas[k]
        end
    end

    return nothing
end

@inline function _compute_full_pair_external_derivatives!(
    rb::RegularizationBuffers,
    q_state::Vector{Float64},
    p_state::Vector{Float64},
    prob::WeberProblem,
)
    system = prob.system
    system.dq_dt_compiled(rb.dq_full, q_state, p_state, prob.params)
    system.dp_dt_compiled(rb.dp_full, q_state, p_state, prob.params)
    system.dq_dt_compiled(rb.dq_pair, q_state, p_state, rb.params_pair)
    system.dp_dt_compiled(rb.dp_pair, q_state, p_state, rb.params_pair)

    @inbounds @. rb.dq_ext = rb.dq_full - rb.dq_pair
    @inbounds @. rb.dp_ext = rb.dp_full - rb.dp_pair

    return nothing
end

@inline function _external_half_step_midpoint!(
    q::Vector{Float64},
    p::Vector{Float64},
    dt_half::Float64,
    prob::WeberProblem,
    rb::RegularizationBuffers,
)
    if dt_half <= 0
        return nothing
    end

    _compute_full_pair_external_derivatives!(rb, q, p, prob)

    midpoint_scale = dt_half * 0.5
    @inbounds @. rb.q_mid = q + midpoint_scale * rb.dq_ext
    @inbounds @. rb.p_mid = p + midpoint_scale * rb.dp_ext

    _compute_full_pair_external_derivatives!(rb, rb.q_mid, rb.p_mid, prob)

    @inbounds @. q = q + dt_half * rb.dq_ext
    @inbounds @. p = p + dt_half * rb.dp_ext

    return nothing
end

@inline function _extract_pair_2d_state!(
    rb::RegularizationBuffers,
    q_state::Vector{Float64},
    p_state::Vector{Float64},
    masses::Vector{Float64},
    i::Int,
    j::Int,
)
    mi = masses[i]
    mj = masses[j]
    M = mi + mj
    mu = mi * mj / M

    i0 = (i - 1) * 2
    j0 = (j - 1) * 2

    @inbounds for d = 1:2
        qi = q_state[i0+d]
        qj = q_state[j0+d]
        pi = p_state[i0+d]
        pj = p_state[j0+d]

        rb.pair_R[d] = (mi * qi + mj * qj) / M
        rb.pair_P[d] = pi + pj

        rb.rel_q[d] = qi - qj
        rb.rel_p[d] = mu * (pi / mi - pj / mj)
    end

    return mi, mj, mu, M
end

@inline function _extract_pair_2d_derivatives!(
    rb::RegularizationBuffers,
    dq_state::Vector{Float64},
    dp_state::Vector{Float64},
    i::Int,
    j::Int,
    mi::Float64,
    mj::Float64,
    mu::Float64,
    M::Float64,
)
    i0 = (i - 1) * 2
    j0 = (j - 1) * 2

    @inbounds for d = 1:2
        dqi = dq_state[i0+d]
        dqj = dq_state[j0+d]
        dpi = dp_state[i0+d]
        dpj = dp_state[j0+d]

        rb.pair_Rdot[d] = (mi * dqi + mj * dqj) / M
        rb.pair_Pdot[d] = dpi + dpj

        rb.rel_qdot[d] = dqi - dqj
        rb.rel_pdot[d] = mu * (dpi / mi - dpj / mj)
    end

    return nothing
end

@inline function _write_pair_state_2d!(
    q_state::Vector{Float64},
    p_state::Vector{Float64},
    i::Int,
    j::Int,
    mi::Float64,
    mj::Float64,
    M::Float64,
    R::Vector{Float64},
    P::Vector{Float64},
    rel_q::Vector{Float64},
    rel_p::Vector{Float64},
)
    i0 = (i - 1) * 2
    j0 = (j - 1) * 2

    αi = mj / M
    αj = mi / M
    βi = mi / M
    βj = mj / M

    @inbounds begin
        q_state[i0+1] = R[1] + αi * rel_q[1]
        q_state[i0+2] = R[2] + αi * rel_q[2]
        q_state[j0+1] = R[1] - αj * rel_q[1]
        q_state[j0+2] = R[2] - αj * rel_q[2]

        p_state[i0+1] = βi * P[1] + rel_p[1]
        p_state[i0+2] = βi * P[2] + rel_p[2]
        p_state[j0+1] = βj * P[1] - rel_p[1]
        p_state[j0+2] = βj * P[2] - rel_p[2]
    end

    return nothing
end

@inline function _current_pair_r_2d(q::Vector{Float64}, i::Int, j::Int)::Float64
    i0 = (i - 1) * 2
    j0 = (j - 1) * 2
    @inbounds begin
        dx = q[i0+1] - q[j0+1]
        dy = q[i0+2] - q[j0+2]
    end
    return sqrt(dx * dx + dy * dy)
end

# Reflect the relative coordinate through the origin while leaving
# momenta unchanged.  Physically this is the two particles passing
# through each other at r = 0 (C⁰-continuation of the ℓ = 0 collision).
# Energy is exactly preserved because H depends only on |q_rel| and p_rel.
@inline function _reflect_pair_2d!(
    q::Vector{Float64},
    p::Vector{Float64},
    masses::Vector{Float64},
    i::Int,
    j::Int,
)
    mi = masses[i]
    mj = masses[j]
    M = mi + mj

    i0 = (i - 1) * 2
    j0 = (j - 1) * 2

    @inbounds begin
        dx = q[i0+1] - q[j0+1]
        dy = q[i0+2] - q[j0+2]

        # q_rel_new = -q_rel, COM unchanged.
        # q_i_new = R + (mj/M)*(-q_rel) = q_i - 2*(mj/M)*q_rel
        # q_j_new = R - (mi/M)*(-q_rel) = q_j + 2*(mi/M)*q_rel
        fi = 2.0 * mj / M
        fj = 2.0 * mi / M

        q[i0+1] -= fi * dx
        q[i0+2] -= fi * dy
        q[j0+1] += fj * dx
        q[j0+2] += fj * dy
    end

    # Momenta unchanged — particles pass through each other.
    return nothing
end

# General-dimension variants for the adaptive_cartesian backend.
@inline function _current_pair_r(
    q::Vector{Float64},
    dims::Int,
    i::Int,
    j::Int,
)::Float64
    i0 = (i - 1) * dims
    j0 = (j - 1) * dims
    r2 = 0.0
    @inbounds for d = 1:dims
        dx = q[i0+d] - q[j0+d]
        r2 += dx * dx
    end
    return sqrt(r2)
end

# Reflect relative coordinate through the origin (general-dimension).
# Preserves COM, negates q_rel, leaves p unchanged.
@inline function _reflect_pair!(
    q::Vector{Float64},
    masses::Vector{Float64},
    dims::Int,
    i::Int,
    j::Int,
)
    mi = masses[i]
    mj = masses[j]
    M = mi + mj
    fi = 2.0 * mj / M
    fj = 2.0 * mi / M

    i0 = (i - 1) * dims
    j0 = (j - 1) * dims

    @inbounds for d = 1:dims
        dx = q[i0+d] - q[j0+d]
        q[i0+d] -= fi * dx
        q[j0+d] += fj * dx
    end

    return nothing
end

@inline function _compute_lc_tau_derivatives!(
    du_tau::Vector{Float64},
    dU_tau::Vector{Float64},
    u::Vector{Float64},
    p_rel::Vector{Float64},
    qdot_rel::Vector{Float64},
    pdot_rel::Vector{Float64},
    g_scale::Float64,
    g_floor::Float64,
)
    u1 = u[1]
    u2 = u[2]

    r = u1 * u1 + u2 * u2
    r_geom = max(r, g_floor)

    inv_2r = 0.5 / r_geom
    du1_dt = inv_2r * (u1 * qdot_rel[1] + u2 * qdot_rel[2])
    du2_dt = inv_2r * (-u2 * qdot_rel[1] + u1 * qdot_rel[2])

    dU1_dt =
        2 * (u1 * pdot_rel[1] + u2 * pdot_rel[2]) +
        2 * (du1_dt * p_rel[1] + du2_dt * p_rel[2])
    dU2_dt =
        2 * (-u2 * pdot_rel[1] + u1 * pdot_rel[2]) +
        2 * (-du2_dt * p_rel[1] + du1_dt * p_rel[2])

    du_tau[1] = g_scale * du1_dt
    du_tau[2] = g_scale * du2_dt
    dU_tau[1] = g_scale * dU1_dt
    dU_tau[2] = g_scale * dU2_dt

    return r_geom
end

@inline function _lifted_pair_substep_2d!(
    q::Vector{Float64},
    p::Vector{Float64},
    dt_sub::Float64,
    prob::WeberProblem,
    rb::RegularizationBuffers,
    i::Int,
    j::Int,
)
    system = prob.system
    masses = prob.masses
    g_floor = prob.regularization.g_floor

    prev_u1 = rb.lc_u[1]
    prev_u2 = rb.lc_u[2]
    prev_u_norm2 = prev_u1 * prev_u1 + prev_u2 * prev_u2

    system.dq_dt_compiled(rb.dq_pair, q, p, rb.params_pair)
    system.dp_dt_compiled(rb.dp_pair, q, p, rb.params_pair)

    mi, mj, mu, M = _extract_pair_2d_state!(rb, q, p, masses, i, j)
    _extract_pair_2d_derivatives!(rb, rb.dq_pair, rb.dp_pair, i, j, mi, mj, mu, M)

    _lc_lift!(rb)
    if prev_u_norm2 > 0
        dot_u = rb.lc_u[1] * prev_u1 + rb.lc_u[2] * prev_u2
        if dot_u < 0
            rb.lc_u[1] = -rb.lc_u[1]
            rb.lc_u[2] = -rb.lc_u[2]
            rb.lc_U[1] = -rb.lc_U[1]
            rb.lc_U[2] = -rb.lc_U[2]
        end
    end
    r_eff = max(rb.lc_u[1] * rb.lc_u[1] + rb.lc_u[2] * rb.lc_u[2], g_floor)
    _compute_lc_tau_derivatives!(
        rb.lc_du_tau,
        rb.lc_dU_tau,
        rb.lc_u,
        rb.rel_p,
        rb.rel_qdot,
        rb.rel_pdot,
        r_eff,
        g_floor,
    )

    dτ = dt_sub / r_eff

    @inbounds for d = 1:2
        rb.lc_u_mid[d] = rb.lc_u[d] + 0.5 * dτ * rb.lc_du_tau[d]
        rb.lc_U_mid[d] = rb.lc_U[d] + 0.5 * dτ * rb.lc_dU_tau[d]

        rb.pair_R_mid[d] = rb.pair_R[d] + 0.5 * dt_sub * rb.pair_Rdot[d]
        rb.pair_P_mid[d] = rb.pair_P[d] + 0.5 * dt_sub * rb.pair_Pdot[d]
    end

    _lc_project!(rb.temp_rel_q, rb.temp_rel_p, rb.lc_u_mid, rb.lc_U_mid)

    copyto!(rb.q_mid, q)
    copyto!(rb.p_mid, p)
    _write_pair_state_2d!(
        rb.q_mid,
        rb.p_mid,
        i,
        j,
        mi,
        mj,
        M,
        rb.pair_R_mid,
        rb.pair_P_mid,
        rb.temp_rel_q,
        rb.temp_rel_p,
    )

    system.dq_dt_compiled(rb.dq_pair, rb.q_mid, rb.p_mid, rb.params_pair)
    system.dp_dt_compiled(rb.dp_pair, rb.q_mid, rb.p_mid, rb.params_pair)
    _extract_pair_2d_derivatives!(rb, rb.dq_pair, rb.dp_pair, i, j, mi, mj, mu, M)

    _compute_lc_tau_derivatives!(
        rb.lc_du_tau_mid,
        rb.lc_dU_tau_mid,
        rb.lc_u_mid,
        rb.temp_rel_p,
        rb.rel_qdot,
        rb.rel_pdot,
        r_eff,
        g_floor,
    )

    @inbounds for d = 1:2
        rb.lc_u[d] = rb.lc_u[d] + dτ * rb.lc_du_tau_mid[d]
        rb.lc_U[d] = rb.lc_U[d] + dτ * rb.lc_dU_tau_mid[d]

        rb.pair_R[d] = rb.pair_R[d] + dt_sub * rb.pair_Rdot[d]
        rb.pair_P[d] = rb.pair_P[d] + dt_sub * rb.pair_Pdot[d]
    end

    _lc_project!(rb.rel_q, rb.rel_p, rb.lc_u, rb.lc_U)
    _write_pair_state_2d!(
        q,
        p,
        i,
        j,
        mi,
        mj,
        M,
        rb.pair_R,
        rb.pair_P,
        rb.rel_q,
        rb.rel_p,
    )

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
    _update_diagnostics!(
        integrator,
        REG_MODE_NONE,
        1,
        Inf,
        0.0,
        REG_BACKEND_DISABLED,
        false,
    )
    return nothing
end

@inline function _step_regularized_pair_adaptive!(
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
    constraint_tolerance = reg.constraint_tolerance
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
            if c_err <= constraint_tolerance
                c_err = 0.0
            end
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
    _update_diagnostics!(
        integrator,
        REG_MODE_PAIR,
        substeps,
        min_distance,
        max_constraint,
        REG_BACKEND_ADAPTIVE,
        rb.backend_fallback,
    )
    return nothing
end

@inline function _step_regularized_pair_lifted_2d!(
    integrator::WeberIntegrator,
    dt_step::Float64,
    min_distance::Float64,
)
    prob = integrator.prob
    rb = integrator.buffers.regularization_buffers

    i = rb.active_anchor_i
    j = rb.active_anchor_j

    _set_pair_params!(rb, prob, i, j)

    reg = prob.regularization
    g_floor = reg.g_floor
    max_sub = reg.max_substeps

    # Bounded fictitious timestep: calibrated so that at r = r_on
    # the physical substep equals dt_step / 1.5 (matching old code).
    dtau_target = dt_step / (1.5 * rb.r_on)

    t_remaining = dt_step
    substeps = 0

    @inbounds while t_remaining > 1e-14 && substeps < max_sub
        r_current = max(_current_pair_r_2d(integrator.q, i, j), g_floor)

        # Physical substep proportional to current r (bounded fictitious step).
        dt_sub = min(r_current * dtau_target, t_remaining)
        # Safety floor: prevent zero-length substeps from stalling.
        dt_sub = max(dt_sub, min(g_floor * dtau_target, t_remaining))

        dt_half = 0.5 * dt_sub

        _external_half_step_midpoint!(integrator.q, integrator.p, dt_half, prob, rb)
        _lifted_pair_substep_2d!(integrator.q, integrator.p, dt_sub, prob, rb, i, j)
        _external_half_step_midpoint!(integrator.q, integrator.p, dt_half, prob, rb)

        t_remaining -= dt_sub
        substeps += 1
    end

    # Fallback: if max_substeps exhausted, complete remaining time in one step.
    if t_remaining > 1e-14
        dt_half = 0.5 * t_remaining
        _external_half_step_midpoint!(integrator.q, integrator.p, dt_half, prob, rb)
        _lifted_pair_substep_2d!(integrator.q, integrator.p, t_remaining, prob, rb, i, j)
        _external_half_step_midpoint!(integrator.q, integrator.p, dt_half, prob, rb)
        substeps += 1
    end

    _store_state!(integrator, dt_step)
    _update_diagnostics!(
        integrator,
        REG_MODE_PAIR,
        substeps,
        min_distance,
        0.0,
        REG_BACKEND_LIFTED,
        rb.backend_fallback,
    )
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
    constraint_tolerance = reg.constraint_tolerance

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
                if c_err <= constraint_tolerance
                    c_err = 0.0
                end
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
    _update_diagnostics!(
        integrator,
        REG_MODE_CHAIN,
        substeps,
        min_distance,
        max_constraint,
        REG_BACKEND_ADAPTIVE,
        false,
    )
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

    if !active || mode == REG_MODE_NONE
        _step_unregularized!(integrator, dt_step)
        return nothing
    end

    if mode == REG_MODE_CHAIN
        _step_regularized_chain!(integrator, dt_step, min_distance)
    elseif rb.effective_backend == REG_BACKEND_LIFTED && rb.dims == 2
        _step_regularized_pair_lifted_2d!(integrator, dt_step, min_distance)
    else
        _step_regularized_pair_adaptive!(integrator, dt_step, min_distance)
    end

    return nothing
end

@inline function _clone_diagnostics(
    diagnostics::RegularizationDiagnostics,
    n_steps::Int,
)
    out = RegularizationDiagnostics(
        diagnostics.enabled,
        n_steps,
        diagnostics.requested_backend,
        diagnostics.used_backend,
    )
    out.activation_count = diagnostics.activation_count
    out.deactivation_count = diagnostics.deactivation_count
    out.active_steps = diagnostics.active_steps
    out.pair_steps = diagnostics.pair_steps
    out.adaptive_pair_steps = diagnostics.adaptive_pair_steps
    out.lifted_pair_steps = diagnostics.lifted_pair_steps
    out.chain_steps = diagnostics.chain_steps
    out.unregularized_steps = diagnostics.unregularized_steps
    out.backend_fallback_steps = diagnostics.backend_fallback_steps
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
    rb = buffers.regularization_buffers

    n_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
    t_history = Vector{Float64}(undef, n_steps + 1)
    q_history = [Vector{Float64}(undef, degrees_of_freedom) for _ = 1:(n_steps+1)]
    p_history = [Vector{Float64}(undef, degrees_of_freedom) for _ = 1:(n_steps+1)]

    t_history[1] = prob.tspan[1]
    q_history[1] .= prob.q_initial
    p_history[1] .= prob.p_initial

    requested_backend = prob.regularization.enabled ? prob.regularization.backend : REG_BACKEND_DISABLED
    used_backend = prob.regularization.enabled ? rb.effective_backend : REG_BACKEND_DISABLED
    diagnostics =
        RegularizationDiagnostics(prob.regularization.enabled, n_steps, requested_backend, used_backend)

    if prob.regularization.enabled && rb.backend_fallback && prob.regularization.warn_on_fallback
        @warn "regularization_backend=:lifted_pair is currently supported only for 2D; falling back to :adaptive_cartesian for $(prob.system.dims)D"
    end

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

# Pre-step collision bounce: for each pair closer than bounce_r,
# reflect the relative coordinate through the origin.
# Used for like-charge sub-critical oscillation (ℓ=0, C⁰-continuable).
@inline function _apply_collision_bounces!(
    q::Vector{Float64},
    masses::Vector{Float64},
    dims::Int,
    n_particles::Int,
    bounce_r::Float64,
)
    @inbounds for i = 1:n_particles
        for j = (i+1):n_particles
            r = _current_pair_r(q, dims, i, j)
            if r < bounce_r
                _reflect_pair!(q, masses, dims, i, j)
            end
        end
    end
    return nothing
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

    # Pre-step collision bounce: reflect pairs that are closer than
    # the bounce radius.  This prevents the integrator from entering
    # the region where the implicit midpoint iteration diverges due
    # to the 1/r² force singularity.
    bounce_r = prob.regularization.collision_bounce_radius
    if bounce_r > 0
        _apply_collision_bounces!(
            integrator.q,
            prob.masses,
            prob.system.dims,
            prob.system.n_particles,
            bounce_r,
        )
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
