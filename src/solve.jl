using CommonSolve

# =============================================================================
# CommonSolve Interface Implementation
# =============================================================================

"""
    init(prob, alg=SymmetricProjectionIntegrator()) -> WeberIntegrator

Initialize an integrator for stepped integration. Use `step!()` to advance, `solve!()` to complete.
"""
function CommonSolve.init(prob::WeberProblem{P}, alg::SymmetricProjectionIntegrator=SymmetricProjectionIntegrator()) where P
    degrees_of_freedom = prob.H.degrees_of_freedom

    # Convert params once, cache in buffers (avoids per-step allocation)
    params_vec = _params_to_vector(prob.params)
    buffers = SymmetricProjectionBuffers(degrees_of_freedom, params_vec)

    # Pre-allocate solution storage with all inner vectors
    n_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
    t_history = Vector{Float64}(undef, n_steps + 1)
    q_history = [Vector{Float64}(undef, degrees_of_freedom) for _ in 1:(n_steps + 1)]
    p_history = [Vector{Float64}(undef, degrees_of_freedom) for _ in 1:(n_steps + 1)]

    # Store initial conditions (in-place copy to pre-allocated vectors)
    t_history[1] = prob.tspan[1]
    q_history[1] .= prob.q_initial
    p_history[1] .= prob.p_initial

    WeberIntegrator{P,typeof(alg)}(
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
        p_history
    )
end

"""
    step!(integrator::WeberIntegrator)

Advance the integrator by one time step.

# Returns
- `true` if step succeeded and integration should continue
- `false` if integration is complete (reached t_end)

# Throws
- `NonlinearSolveError` if the nonlinear solver fails to converge
"""
function CommonSolve.step!(integrator::WeberIntegrator{P,SymmetricProjectionIntegrator{S}}) where {P,S}
    # Check if already done
    if integrator.t >= integrator.t_end - eps(integrator.t_end)
        return false
    end

    prob = integrator.prob
    dt = prob.dt
    convergence_tolerance = prob.convergence_tolerance
    maximum_iterations = prob.maximum_iterations

    dq_dt_compiled = prob.H.dq_dt_compiled
    dp_dt_compiled = prob.H.dp_dt_compiled
    params = integrator.buffers.params_vec  # Use cached params (no allocation)

    buffers = integrator.buffers
    degrees_of_freedom = buffers.degrees_of_freedom

    # Extended phase space indices: Z = [Q, X, P, Y] where each component has d elements
    # Q = positions, X = auxiliary positions, P = momenta, Y = auxiliary momenta
    idx_Q_start = 1
    idx_Q_end = degrees_of_freedom
    idx_X_start = degrees_of_freedom + 1
    idx_X_end = 2 * degrees_of_freedom
    idx_P_start = 2 * degrees_of_freedom + 1
    idx_P_end = 3 * degrees_of_freedom
    idx_Y_start = 3 * degrees_of_freedom + 1
    idx_Y_end = 4 * degrees_of_freedom

    constraint_matrix = buffers.constraint_matrix
    extended_state = buffers.extended_state
    extended_state_after_flow = buffers.extended_state_after_flow
    extended_state_result = buffers.extended_state_result
    position_buffer = buffers.position_buffer
    auxiliary_position_buffer = buffers.auxiliary_position_buffer
    momentum_buffer = buffers.momentum_buffer
    auxiliary_momentum_buffer = buffers.auxiliary_momentum_buffer
    constraint_shift = buffers.constraint_shift
    lagrange_multipliers = buffers.lagrange_multipliers
    lagrange_multipliers_previous = buffers.lagrange_multipliers_previous
    residual_buffer = buffers.residual_buffer

    q = integrator.q
    p = integrator.p

    @views begin
        extended_state[idx_Q_start:idx_Q_end] .= q
        extended_state[idx_X_start:idx_X_end] .= q
        extended_state[idx_P_start:idx_P_end] .= p
        extended_state[idx_Y_start:idx_Y_end] .= p
    end

    # Strang splitting flow map ϕ
    function flow_map(Z::Vector{Float64})::Nothing
        @views begin
            Q_component = Z[idx_Q_start:idx_Q_end]
            X_component = Z[idx_X_start:idx_X_end]
            P_component = Z[idx_P_start:idx_P_end]
            Y_component = Z[idx_Y_start:idx_Y_end]

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

    # Projection constraint function
    function projection_residual!(result::Vector{Float64}, multipliers::Vector{Float64})::Nothing
        mul!(constraint_shift, constraint_matrix', multipliers)
        @. extended_state_after_flow = extended_state + constraint_shift
        flow_map(extended_state_after_flow)
        @. extended_state_result = extended_state_after_flow + constraint_shift
        mul!(result, constraint_matrix, extended_state_result)
        return nothing
    end

    # Solve the nonlinear projection problem
    solver = integrator.alg.solver
    success, iterations, final_residual = solve_nonlinear!(
        lagrange_multipliers, projection_residual!, solver, convergence_tolerance, maximum_iterations;
        fu_buffer=residual_buffer,
        u_old_buffer=lagrange_multipliers_previous
    )

    if success
        # Recompute extended_state_after_flow with the final converged lagrange_multipliers to ensure consistency
        # (The solver's last projection_residual! call used previous multipliers, not final)
        mul!(constraint_shift, constraint_matrix', lagrange_multipliers)
        @. extended_state_after_flow = extended_state + constraint_shift
        flow_map(extended_state_after_flow)

        # Apply symmetric projection shift: Z_{n+1} = ϕ(Z_n + A'μ) + A'μ
        @. extended_state_after_flow = extended_state_after_flow + constraint_shift

        @views begin
            integrator.q .= extended_state_after_flow[idx_Q_start:idx_Q_end]
            integrator.p .= extended_state_after_flow[idx_P_start:idx_P_end]
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

    solver_name = string(typeof(solver).name.name)
    throw(NonlinearSolveError(iterations, convergence_tolerance, final_residual, integrator.step_count + 1, integrator.t, solver_name))
end

"""
    solve!(integrator::WeberIntegrator)

Complete integration from current state to t_end.

Returns a `WeberSolution` containing the full trajectory.
"""
function CommonSolve.solve!(integrator::WeberIntegrator)
    retcode = :Success
    try
        while CommonSolve.step!(integrator)
            # Continue stepping
        end
    catch e
        if e isa NonlinearSolveError
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
        retcode
    )
end

"""
    solve(prob, alg=SymmetricProjectionIntegrator()) -> WeberSolution

Solve a Weber electrodynamics problem and return the full trajectory.
"""
function CommonSolve.solve(prob::WeberProblem, alg::WeberAlgorithm=SymmetricProjectionIntegrator())
    integrator = CommonSolve.init(prob, alg)
    CommonSolve.solve!(integrator)
end

# =============================================================================
# Helper Functions
# =============================================================================

"""Convert params to Vector{Float64} for compiled functions."""
function _params_to_vector(params::Vector{Float64})
    params
end

function _params_to_vector(params::AbstractVector)
    Vector{Float64}(params)
end

function _params_to_vector(params::NamedTuple)
    Vector{Float64}(collect(values(params)))
end
