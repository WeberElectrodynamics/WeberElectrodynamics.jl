using CommonSolve

# =============================================================================
# CommonSolve Interface Implementation
# =============================================================================

"""
    init(prob::WeberProblem, alg::WeberAlgorithm=SymmetricProjection())

Initialize an integrator for stepped integration.

Returns a `WeberIntegrator` that can be advanced with `step!()` or
completed with `solve!()`.

# Example
```julia
integrator = init(prob, SymmetricProjection())
for i in 1:100
    step!(integrator)
    println("t = \$(integrator.t), q = \$(integrator.q)")
end
sol = solve!(integrator)
```
"""
function CommonSolve.init(prob::WeberProblem{P}, alg::SymmetricProjection=SymmetricProjection()) where P
    d = prob.H.n_dof
    buffers = IntegratorBuffers(d)

    # Pre-allocate solution storage
    n_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
    t_history = Vector{Float64}(undef, n_steps + 1)
    q_history = Vector{Vector{Float64}}(undef, n_steps + 1)
    p_history = Vector{Vector{Float64}}(undef, n_steps + 1)

    # Store initial conditions
    t_history[1] = prob.tspan[1]
    q_history[1] = copy(prob.q₀)
    p_history[1] = copy(prob.p₀)

    WeberIntegrator{P,typeof(alg)}(
        prob,
        alg,
        prob.tspan[1],
        prob.tspan[2],
        copy(prob.q₀),
        copy(prob.p₀),
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
function CommonSolve.step!(integrator::WeberIntegrator{P,SymmetricProjection{S}}) where {P,S}
    # Check if already done
    if integrator.t >= integrator.t_end - eps(integrator.t_end)
        return false
    end

    prob = integrator.prob
    dt = prob.dt
    tol = prob.tolerance
    max_iter = prob.max_iterations

    qdot_func = prob.H.qdot_func
    pdot_func = prob.H.pdot_func
    params = _params_to_vector(prob.params)

    buffers = integrator.buffers
    d = buffers.d

    # Extended phase space indices
    iQ = 1
    iQ_ = d
    iX = d + 1
    iX_ = 2d
    iP = 2d + 1
    iP_ = 3d
    iY = 3d + 1
    iY_ = 4d

    A = buffers.A
    Z_current = buffers.Z_current
    Z_post_phi = buffers.Z_post_phi
    Z_result = buffers.Z_result
    qs_buf = buffers.qs_buf
    xs_buf = buffers.xs_buf
    ps_buf = buffers.ps_buf
    ys_buf = buffers.ys_buf
    ATμ = buffers.ATμ
    μ = buffers.μ
    μ_old = buffers.μ_old
    f_val = buffers.f_val

    q = integrator.q
    p = integrator.p

    @views begin
        Z_current[iQ:iQ_] .= q
        Z_current[iX:iX_] .= q
        Z_current[iP:iP_] .= p
        Z_current[iY:iY_] .= p
    end

    # Strang splitting flow map ϕ
    function ϕ(Z::Vector{Float64})::Nothing
        @views begin
            qs, xs, ps, ys = Z[iQ:iQ_], Z[iX:iX_], Z[iP:iP_], Z[iY:iY_]

            qdot_func(xs_buf, qs, ys, params)
            pdot_func(ps_buf, qs, ys, params)

            @. xs = xs + xs_buf * (dt / 2)
            @. ps = ps + ps_buf * (dt / 2)

            qdot_func(qs_buf, xs, ps, params)
            pdot_func(ys_buf, xs, ps, params)

            @. qs = qs + qs_buf * dt
            @. ys = ys + ys_buf * dt

            qdot_func(xs_buf, qs, ys, params)
            pdot_func(ps_buf, qs, ys, params)

            @. xs = xs + xs_buf * (dt / 2)
            @. ps = ps + ps_buf * (dt / 2)
        end
        return nothing
    end

    # Projection constraint function
    function f!(result::Vector{Float64}, μ_input::Vector{Float64})::Nothing
        mul!(ATμ, A', μ_input)
        @. Z_post_phi = Z_current + ATμ
        ϕ(Z_post_phi)
        @. Z_result = Z_post_phi + ATμ
        mul!(result, A, Z_result)
        return nothing
    end

    # Solve the nonlinear projection problem
    solver = integrator.alg.solver
    success, iterations, final_residual = solve_nonlinear!(
        μ, f!, solver, tol, max_iter;
        fu_buffer=f_val,
        u_old_buffer=μ_old
    )

    if success
        mul!(ATμ, A', μ)
        @. Z_post_phi = Z_post_phi + ATμ

        @views begin
            integrator.q .= Z_post_phi[iQ:iQ_]
            integrator.p .= Z_post_phi[iP:iP_]
        end

        integrator.step_count += 1
        integrator.t = prob.tspan[1] + integrator.step_count * dt

        # Store in history
        idx = integrator.step_count + 1
        integrator.t_history[idx] = integrator.t
        integrator.q_history[idx] = copy(integrator.q)
        integrator.p_history[idx] = copy(integrator.p)

        return integrator.t < integrator.t_end
    end

    solver_name = string(typeof(solver).name.name)
    throw(NonlinearSolveError(iterations, tol, final_residual, integrator.step_count + 1, integrator.t, solver_name))
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
    solve(prob::WeberProblem, alg::WeberAlgorithm=SymmetricProjection())

Solve a Weber electrodynamics problem.

# Arguments
- `prob`: Problem specification
- `alg`: Integration algorithm (default: `SymmetricProjection()`)

# Returns
`WeberSolution` containing the full trajectory.

# Example
```julia
H = @hamiltonian (q, p, params) -> begin
    m1, m2, k, c = params
    # ... Hamiltonian expression
end

prob = WeberProblem(H, (0.0, 10.0), q₀, p₀; params=[1.0, 0.1, 0.1, 4.0], dt=0.01)
sol = solve(prob)
sol = solve(prob, SymmetricProjection())  # Explicit algorithm
```
"""
function CommonSolve.solve(prob::WeberProblem, alg::WeberAlgorithm=SymmetricProjection())
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
