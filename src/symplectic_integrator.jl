using LinearAlgebra

mutable struct PhaseSpacePoint
    q::Vector{Float64}
    p::Vector{Float64}
end

struct IntegratorSettings
    dt::Float64
    tolerance::Float64
    max_iterations::Int

    function IntegratorSettings(dt::Float64, tolerance::Float64=1e-13, max_iterations::Int=100)
        new(dt, tolerance, max_iterations)
    end
end

mutable struct IntegratorBuffers
    d::Int
    A::Matrix{Float64}
    Zₙ::Vector{Float64}
    Ẑₙ::Vector{Float64}
    Ẑₙ₊₁::Vector{Float64}
    Z_result::Vector{Float64}
    qs_buf::Vector{Float64}
    xs_buf::Vector{Float64}
    ps_buf::Vector{Float64}
    ys_buf::Vector{Float64}
    ATμ::Vector{Float64}
    μ::Vector{Float64}
    μ_old::Vector{Float64}
    f_val::Vector{Float64}

    function IntegratorBuffers(d::Int)
        Id = Matrix{Float64}(I, d, d)
        Zd = zeros(Float64, d, d)
        A = [Id -Id Zd Zd;
             Zd Zd Id -Id]

        new(
            d,
            A,
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, 4d),
            zeros(Float64, 2d),
            Vector{Float64}(undef, 2d),
            Vector{Float64}(undef, 2d)
        )
    end
end

mutable struct IntegratorState
    point::PhaseSpacePoint
    settings::IntegratorSettings
    vector_field::HamiltonianVectorField
    params::Vector{Float64}
    buffers::IntegratorBuffers

    function IntegratorState(point::PhaseSpacePoint, settings::IntegratorSettings,
                             vector_field::HamiltonianVectorField, params::Vector{Float64})
        d = length(point.q)
        buffers = IntegratorBuffers(d)
        new(point, settings, vector_field, params, buffers)
    end
end

struct TimeSpan
    start::Float64
    stop::Float64
end

struct IntegratorSolution
    t::Vector{Float64}
    q::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
end

function step!(state::IntegratorState)::Nothing
    q = state.point.q
    p = state.point.p

    dt = state.settings.dt
    tol = state.settings.tolerance
    max_iter = state.settings.max_iterations

    qdot_func = state.vector_field.qdot_func
    pdot_func = state.vector_field.pdot_func

    params = state.params

    buffers = state.buffers
    d = buffers.d

    Q = 1
    Q_ = d
    X = d + 1
    X_ = 2d
    P = 2d + 1
    P_ = 3d
    Y = 3d + 1
    Y_ = 4d

    A = buffers.A
    Zₙ = buffers.Zₙ
    Ẑₙ = buffers.Ẑₙ
    Ẑₙ₊₁ = buffers.Ẑₙ₊₁
    qs_buf = buffers.qs_buf
    xs_buf = buffers.xs_buf
    ps_buf = buffers.ps_buf
    ys_buf = buffers.ys_buf
    ATμ = buffers.ATμ
    μ = buffers.μ
    μ_old = buffers.μ_old
    f_val = buffers.f_val
    Z_result = buffers.Z_result

    @views begin
        Zₙ[Q:Q_] .= q
        Zₙ[X:X_] .= q
        Zₙ[P:P_] .= p
        Zₙ[Y:Y_] .= p
    end

    function ϕ(Z::Vector{Float64})::Nothing
        @views begin
            qs, xs, ps, ys = Z[Q:Q_], Z[X:X_], Z[P:P_], Z[Y:Y_]

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

    function f!(result::Vector{Float64}, μ_input::Vector{Float64})::Nothing
        mul!(ATμ, A', μ_input)
        @. Ẑₙ = Zₙ + ATμ
        Ẑₙ₊₁ .= Ẑₙ
        ϕ(Ẑₙ₊₁)
        @. Z_result = Ẑₙ₊₁ + ATμ
        mul!(result, A, Z_result)
        return nothing
    end

    for _ in 1:max_iter
        f!(f_val, μ)

        copyto!(μ_old, μ)
        @. μ = μ - 0.25 * f_val

        diff_norm_sq = 0.0
        for i in eachindex(μ, μ_old)
            diff_norm_sq += (μ[i] - μ_old[i])^2
        end

        if sqrt(diff_norm_sq) < tol
            mul!(ATμ, A', μ_old)
            @. Ẑₙ₊₁ = Ẑₙ₊₁ + ATμ

            @views begin
                state.point.q .= Ẑₙ₊₁[Q:Q_]
                state.point.p .= Ẑₙ₊₁[P:P_]
            end

            return nothing
        end
    end

    error("Newton iteration did not converge within $max_iter iterations")
end

function integrate(state::IntegratorState, timespan::TimeSpan; reset_mu::Bool=false)::IntegratorSolution
    dt = state.settings.dt
    t_start = timespan.start
    t_stop = timespan.stop

    n_steps = Int(ceil((t_stop - t_start) / dt))

    if reset_mu
        fill!(state.buffers.μ, 0.0)
    end

    t = Vector{Float64}(undef, n_steps + 1)
    q = Vector{Vector{Float64}}(undef, n_steps + 1)
    p = Vector{Vector{Float64}}(undef, n_steps + 1)

    t[1] = t_start
    q[1] = copy(state.point.q)
    p[1] = copy(state.point.p)

    for step in 1:n_steps
        try
            step!(state)
        catch e
            if e isa ErrorException && contains(e.msg, "Newton iteration")
                error("Integration failed at step $step (t=$(t_start + (step-1)*dt)): $(e.msg)")
            else
                rethrow()
            end
        end

        t[step+1] = t_start + step * dt
        q[step+1] = copy(state.point.q)
        p[step+1] = copy(state.point.p)
    end

    return IntegratorSolution(t, q, p)
end
