using LinearAlgebra

abstract type WeberAlgorithm end

struct SymmetricProjectionIntegrator <: WeberAlgorithm
    relaxation::Float64

    function SymmetricProjectionIntegrator(; relaxation::Real = 0.25)
        @assert 0 < relaxation <= 1 "Relaxation must be in (0, 1], got $relaxation"
        new(Float64(relaxation))
    end
end

struct WeberProblem
    system::WeberSystem
    tspan::Tuple{Float64,Float64}
    q_initial::Vector{Float64}
    p_initial::Vector{Float64}
    masses::Vector{Float64}
    charges::Vector{Float64}
    c::Float64
    params::Vector{Float64}
    dt::Float64
    convergence_tolerance::Float64
    maximum_iterations::Int

    function WeberProblem(
        system::WeberSystem,
        tspan::Tuple{Real,Real},
        q_initial::AbstractVector,
        p_initial::AbstractVector;
        masses::AbstractVector{<:Real},
        charges::AbstractVector{<:Real},
        c::Real,
        dt::Real,
        convergence_tolerance::Real = 1e-13,
        maximum_iterations::Integer = 100,
    )
        n_particles = system.n_particles
        dof = system.degrees_of_freedom
        @assert length(q_initial) == dof "q_initial must have length $dof (got $(length(q_initial)))"
        @assert length(p_initial) == dof "p_initial must have length $dof (got $(length(p_initial)))"
        @assert length(masses) == n_particles "masses must have length $n_particles (got $(length(masses)))"
        @assert length(charges) == n_particles "charges must have length $n_particles (got $(length(charges)))"
        @assert all(m -> m > 0, masses) "All masses must be positive"
        @assert c > 0 "Speed of light must be positive"
        @assert tspan[2] > tspan[1] "End time must be greater than start time"
        @assert dt > 0 "Time step must be positive"
        @assert convergence_tolerance > 0 "Convergence tolerance must be positive"
        @assert maximum_iterations > 0 "Maximum iterations must be positive"

        masses_f64 = Vector{Float64}(masses)
        charges_f64 = Vector{Float64}(charges)
        c_f64 = Float64(c)
        params = vcat(masses_f64, charges_f64, [c_f64])

        new(
            system,
            (Float64(tspan[1]), Float64(tspan[2])),
            Vector{Float64}(q_initial),
            Vector{Float64}(p_initial),
            masses_f64,
            charges_f64,
            c_f64,
            params,
            Float64(dt),
            Float64(convergence_tolerance),
            Int(maximum_iterations),
        )
    end
end

struct WeberSolution
    t::Vector{Float64}
    q::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    prob::WeberProblem
    retcode::Symbol
end

Base.length(sol::WeberSolution) = length(sol.t)
Base.getindex(sol::WeberSolution, i::Int) = (sol.t[i], sol.q[i], sol.p[i])
Base.firstindex(sol::WeberSolution) = 1
Base.lastindex(sol::WeberSolution) = length(sol)

function Base.iterate(sol::WeberSolution, state = 1)
    state > length(sol) && return nothing
    return (sol[state], state + 1)
end

function Base.show(io::IO, sol::WeberSolution)
    print(io, "WeberSolution with $(length(sol)) timesteps (retcode: $(sol.retcode))")
end

function Base.show(io::IO, ::MIME"text/plain", sol::WeberSolution)
    println(io, "WeberSolution")
    println(io, "  retcode: $(sol.retcode)")
    println(io, "  t: $(sol.t[1]) → $(sol.t[end]) ($(length(sol)) points)")
    println(io, "  DOF: $(length(sol.q[1]))")
end

mutable struct SymmetricProjectionBuffers
    d::Int                          # degrees of freedom
    A::Matrix{Float64}              # projection matrix (2d × 4d)
    A_transpose::Transpose{Float64,Matrix{Float64}}  # cached transpose view (4d × 2d)
    Z::Vector{Float64}              # extended state Zₙ (4d)
    Ẑ::Vector{Float64}              # shifted/evolved state (4d)
    Z_result::Vector{Float64}       # final projected state Zₙ₊₁ (4d)
    position_buffer::Vector{Float64}
    auxiliary_position_buffer::Vector{Float64}
    momentum_buffer::Vector{Float64}
    auxiliary_momentum_buffer::Vector{Float64}
    ATμ::Vector{Float64}            # constraint shift A^T μ (4d)
    μ::Vector{Float64}              # Lagrange multipliers (2d)
    μ_prev::Vector{Float64}         # previous iteration μ^(k-1) (2d)
    f_μ::Vector{Float64}            # nonlinear residual f(μ) (2d)
    diff_buffer::Vector{Float64}    # workspace for μ - μ_prev (2d)

    function SymmetricProjectionBuffers(degrees_of_freedom::Int)
        d = degrees_of_freedom
        # Construct projection matrix A directly without intermediate allocations
        A = zeros(Float64, 2d, 4d)
        @inbounds for i = 1:d
            A[i, i] = 1.0           # I_d in top-left
            A[i, d+i] = -1.0      # -I_d in top-middle-left
            A[d+i, 2d+i] = 1.0  # I_d in bottom-middle-right
            A[d+i, 3d+i] = -1.0 # -I_d in bottom-right
        end
        A_transpose = transpose(A)

        new(
            d,
            A,
            A_transpose,
            Vector{Float64}(undef, 4d),  # Z
            Vector{Float64}(undef, 4d),  # Ẑ
            Vector{Float64}(undef, 4d),  # Z_result
            Vector{Float64}(undef, d),   # position_buffer
            Vector{Float64}(undef, d),   # auxiliary_position_buffer
            Vector{Float64}(undef, d),   # momentum_buffer
            Vector{Float64}(undef, d),   # auxiliary_momentum_buffer
            Vector{Float64}(undef, 4d),  # ATμ
            zeros(Float64, 2d),          # μ (initialized to zero)
            Vector{Float64}(undef, 2d),  # μ_prev
            Vector{Float64}(undef, 2d),  # f_μ
            Vector{Float64}(undef, 2d),  # diff_buffer
        )
    end
end

mutable struct WeberIntegrator
    prob::WeberProblem
    alg::SymmetricProjectionIntegrator
    t::Float64
    t_end::Float64
    q::Vector{Float64}
    p::Vector{Float64}
    step_count::Int
    buffers::SymmetricProjectionBuffers
    # Solution accumulator
    t_history::Vector{Float64}
    q_history::Vector{Vector{Float64}}
    p_history::Vector{Vector{Float64}}
end

function Base.show(io::IO, int::WeberIntegrator)
    print(io, "WeberIntegrator at t=$(int.t) (step $(int.step_count))")
end
