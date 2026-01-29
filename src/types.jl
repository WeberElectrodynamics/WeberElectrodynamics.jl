using LinearAlgebra

# =============================================================================
# Algorithm Types
# =============================================================================

"""
    WeberAlgorithm

Abstract type for Weber electrodynamics integration algorithms.
Extensible for future methods (ImplicitMidpoint, StormerVerlet, etc.).
"""
abstract type WeberAlgorithm end

"""
    SymmetricProjectionIntegrator <: WeberAlgorithm

Semi-explicit symplectic integrator using symmetric projection for non-separable Hamiltonians.
Second-order accurate. Uses relaxed fixed-point iteration for the nonlinear projection solve.

# Constructor
    SymmetricProjectionIntegrator(; relaxation=0.25)

# Arguments
- `relaxation::Float64`: Relaxation parameter for fixed-point iteration, must be in (0, 1].
  Default 0.25. Higher values converge faster for well-conditioned problems but may diverge.
"""
struct SymmetricProjectionIntegrator <: WeberAlgorithm
    relaxation::Float64

    function SymmetricProjectionIntegrator(; relaxation::Real=0.25)
        @assert 0 < relaxation <= 1 "Relaxation must be in (0, 1], got $relaxation"
        new(Float64(relaxation))
    end
end

# =============================================================================
# Hamiltonian Type
# =============================================================================

"""
    WeberHamiltonian

Compiled Hamiltonian system with symbolic expressions and fast numeric functions.
Contains both the symbolic representation (for analysis) and compiled functions
(for integration).

# Fields
- `hamiltonian_symbolic`: Symbolic Hamiltonian expression H(q, p)
- `dq_dt_symbolic`: Symbolic dq/dt = ∂H/∂p (Hamilton's first equation)
- `dp_dt_symbolic`: Symbolic dp/dt = -∂H/∂q (Hamilton's second equation)
- `dq_dt_compiled`: Compiled dq/dt function with signature `dq_dt_compiled(out, q, p, params)`
- `dp_dt_compiled`: Compiled dp/dt function with signature `dp_dt_compiled(out, q, p, params)`
- `degrees_of_freedom`: Number of degrees of freedom (dimension of q or p)
- `parameter_names`: Names of parameters (for documentation)
"""
struct WeberHamiltonian{H,QD,PD,QF,PF}
    hamiltonian_symbolic::H
    dq_dt_symbolic::QD
    dp_dt_symbolic::PD
    dq_dt_compiled::QF
    dp_dt_compiled::PF
    degrees_of_freedom::Int
    parameter_names::Vector{Symbol}
end

# =============================================================================
# Problem Type
# =============================================================================

"""
    WeberProblem{P}

Problem specification for Weber electrodynamics simulation.

    WeberProblem(H, tspan, q_initial, p_initial; params, dt, convergence_tolerance=1e-13, maximum_iterations=100)

`H` is a compiled `WeberHamiltonian`, `tspan = (t_start, t_end)`, `params` is Vector or NamedTuple.
"""
struct WeberProblem{P}
    H::WeberHamiltonian
    tspan::Tuple{Float64,Float64}
    q_initial::Vector{Float64}
    p_initial::Vector{Float64}
    params::P
    dt::Float64
    convergence_tolerance::Float64
    maximum_iterations::Int

    function WeberProblem(
        H::WeberHamiltonian,
        tspan::Tuple{Real,Real},
        q_initial::AbstractVector,
        p_initial::AbstractVector;
        params,
        dt::Real,
        convergence_tolerance::Real=1e-13,
        maximum_iterations::Integer=100
    )
        @assert length(q_initial) == length(p_initial) == H.degrees_of_freedom "Initial conditions must match Hamiltonian DOF ($(H.degrees_of_freedom))"
        @assert tspan[2] > tspan[1] "End time must be greater than start time"
        @assert dt > 0 "Time step must be positive"
        @assert convergence_tolerance > 0 "Convergence tolerance must be positive"
        @assert maximum_iterations > 0 "Maximum iterations must be positive"

        P = typeof(params)
        new{P}(
            H,
            (Float64(tspan[1]), Float64(tspan[2])),
            Vector{Float64}(q_initial),
            Vector{Float64}(p_initial),
            params,
            Float64(dt),
            Float64(convergence_tolerance),
            Int(maximum_iterations)
        )
    end
end

# =============================================================================
# Solution Type
# =============================================================================

"""
    WeberSolution{P}

Solution of a Weber electrodynamics problem.

# Fields
- `t`, `q`, `p`: Time points, coordinates, momenta
- `prob`: Reference to the problem
- `retcode`: `:Success`, `:Failure`, or `:MaxIters`

Supports indexing (`sol[i]`) and iteration.
"""
struct WeberSolution{P}
    t::Vector{Float64}
    q::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    prob::WeberProblem{P}
    retcode::Symbol
end

Base.length(sol::WeberSolution) = length(sol.t)
Base.getindex(sol::WeberSolution, i::Int) = (sol.t[i], sol.q[i], sol.p[i])
Base.firstindex(sol::WeberSolution) = 1
Base.lastindex(sol::WeberSolution) = length(sol)

function Base.iterate(sol::WeberSolution, state=1)
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

# =============================================================================
# Integrator Buffers (Internal)
# =============================================================================

"""
    SymmetricProjectionBuffers

Pre-allocated workspace for the symmetric projection integrator.
Internal type - not part of public API.

Field names match notation in docs/theory/SemiExplicitIntegrator.md (Jayawardana-Ohsawa 2023):
- `d`: degrees of freedom
- `A`: projection matrix (2d × 4d)
- `Z`: extended phase space state Zₙ = (q,q,p,p)
- `Ẑ`: shifted/evolved state Ẑₙ₊₁
- `Z_result`: final projected state Zₙ₊₁
- `ATμ`: constraint shift A^T μ
- `μ`: Lagrange multipliers for symmetric projection
- `f_μ`: nonlinear residual f(μ) = A(Φ̂(Zₙ + A^T μ) + A^T μ)
"""
mutable struct SymmetricProjectionBuffers
    d::Int                          # degrees of freedom
    A::Matrix{Float64}              # projection matrix (2d × 4d)
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
    params_vec::Vector{Float64}     # cached parameter vector

    function SymmetricProjectionBuffers(degrees_of_freedom::Int, params_vec::Vector{Float64})
        d = degrees_of_freedom
        # Construct projection matrix A directly without intermediate allocations
        A = zeros(Float64, 2d, 4d)
        @inbounds for i in 1:d
            A[i, i] = 1.0           # I_d in top-left
            A[i, d + i] = -1.0      # -I_d in top-middle-left
            A[d + i, 2d + i] = 1.0  # I_d in bottom-middle-right
            A[d + i, 3d + i] = -1.0 # -I_d in bottom-right
        end

        new(
            d,
            A,
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
            params_vec
        )
    end
end

# =============================================================================
# Integrator Type (Iterator for CommonSolve)
# =============================================================================

"""
    WeberIntegrator{P}

Mutable integrator state for stepped integration via CommonSolve interface.
Access `t`, `q`, `p`, `step_count` during integration. Use `step!()` to advance, `solve!()` to complete.
"""
mutable struct WeberIntegrator{P}
    prob::WeberProblem{P}
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

