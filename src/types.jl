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
    SymmetricProjectionIntegrator{S} <: WeberAlgorithm

Semi-explicit symplectic integrator using symmetric projection for non-separable Hamiltonians.
Second-order accurate. Type parameter `S` is the nonlinear solver (default: `RelaxedFixedPointSolver`).
"""
struct SymmetricProjectionIntegrator{S} <: WeberAlgorithm
    solver::S
end

SymmetricProjectionIntegrator(; solver=RelaxedFixedPointSolver()) = SymmetricProjectionIntegrator(solver)

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
"""
mutable struct SymmetricProjectionBuffers
    degrees_of_freedom::Int
    constraint_matrix::Matrix{Float64}
    extended_state::Vector{Float64}
    extended_state_after_flow::Vector{Float64}
    extended_state_result::Vector{Float64}
    position_buffer::Vector{Float64}
    auxiliary_position_buffer::Vector{Float64}
    momentum_buffer::Vector{Float64}
    auxiliary_momentum_buffer::Vector{Float64}
    constraint_shift::Vector{Float64}
    lagrange_multipliers::Vector{Float64}
    lagrange_multipliers_previous::Vector{Float64}
    residual_buffer::Vector{Float64}
    params_vec::Vector{Float64}  # Cached parameter vector (avoids per-step allocation)

    function SymmetricProjectionBuffers(degrees_of_freedom::Int, params_vec::Vector{Float64})
        d = degrees_of_freedom
        # Construct constraint matrix directly without intermediate allocations
        constraint_matrix = zeros(Float64, 2d, 4d)
        @inbounds for i in 1:d
            constraint_matrix[i, i] = 1.0           # Id in top-left
            constraint_matrix[i, d + i] = -1.0      # -Id in top-middle-left
            constraint_matrix[d + i, 2d + i] = 1.0  # Id in bottom-middle-right
            constraint_matrix[d + i, 3d + i] = -1.0 # -Id in bottom-right
        end

        new(
            d,
            constraint_matrix,
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
            Vector{Float64}(undef, 2d),
            params_vec
        )
    end
end

# =============================================================================
# Integrator Type (Iterator for CommonSolve)
# =============================================================================

"""
    WeberIntegrator{P,A<:WeberAlgorithm}

Mutable integrator state for stepped integration via CommonSolve interface.
Access `t`, `q`, `p`, `step_count` during integration. Use `step!()` to advance, `solve!()` to complete.
"""
mutable struct WeberIntegrator{P,A<:WeberAlgorithm}
    prob::WeberProblem{P}
    alg::A
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

# =============================================================================
# Error Types
# =============================================================================

"""
    NonlinearSolveError <: Exception

Thrown when the nonlinear solver fails to converge during an integration step.

# Fields
- `iterations`: Number of iterations performed
- `convergence_tolerance`: Target convergence tolerance
- `final_residual`: Residual at termination
- `step`: Integration step number where failure occurred
- `time`: Simulation time at failure
- `solver_name`: Name of the solver that failed
"""
struct NonlinearSolveError <: Exception
    iterations::Int
    convergence_tolerance::Float64
    final_residual::Float64
    step::Int
    time::Float64
    solver_name::String
end

function Base.showerror(io::IO, e::NonlinearSolveError)
    print(io, "NonlinearSolveError: $(e.solver_name) failed to converge at step $(e.step) (t=$(e.time)) ")
    print(io, "after $(e.iterations) iterations (residual=$(e.final_residual), convergence_tolerance=$(e.convergence_tolerance))")
end

# Deprecated alias for backward compatibility
const NewtonConvergenceError = NonlinearSolveError
