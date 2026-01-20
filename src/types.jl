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
    SymmetricProjection{S} <: WeberAlgorithm

Semi-explicit symplectic integrator using symmetric projection for non-separable Hamiltonians.
Second-order accurate. Type parameter `S` is the nonlinear solver (default: `RelaxedFixedPoint`).
"""
struct SymmetricProjection{S} <: WeberAlgorithm
    solver::S
end

SymmetricProjection(; solver=RelaxedFixedPoint()) = SymmetricProjection(solver)

# =============================================================================
# Hamiltonian Type
# =============================================================================

"""
    WeberHamiltonian

Compiled Hamiltonian system with symbolic expressions and fast numeric functions.
Contains both the symbolic representation (for analysis) and compiled functions
(for integration).

# Fields
- `H_sym`: Symbolic Hamiltonian expression
- `qdot_sym`: Symbolic q̇ = ∂H/∂p
- `pdot_sym`: Symbolic ṗ = -∂H/∂q
- `qdot_func`: Compiled q̇ function with signature `qdot_func(out, q, p, params)`
- `pdot_func`: Compiled ṗ function with signature `pdot_func(out, q, p, params)`
- `n_dof`: Number of degrees of freedom
- `param_names`: Names of parameters (for documentation)
"""
struct WeberHamiltonian{H,QD,PD,QF,PF}
    H_sym::H
    qdot_sym::QD
    pdot_sym::PD
    qdot_func::QF
    pdot_func::PF
    n_dof::Int
    param_names::Vector{Symbol}
end

# =============================================================================
# Problem Type
# =============================================================================

"""
    WeberProblem{P}

Problem specification for Weber electrodynamics simulation.

    WeberProblem(H, tspan, q₀, p₀; params, dt, tolerance=1e-13, max_iterations=100)

`H` is a compiled `WeberHamiltonian`, `tspan = (t_start, t_end)`, `params` is Vector or NamedTuple.
"""
struct WeberProblem{P}
    H::WeberHamiltonian
    tspan::Tuple{Float64,Float64}
    q₀::Vector{Float64}
    p₀::Vector{Float64}
    params::P
    dt::Float64
    tolerance::Float64
    max_iterations::Int

    function WeberProblem(
        H::WeberHamiltonian,
        tspan::Tuple{Real,Real},
        q₀::AbstractVector,
        p₀::AbstractVector;
        params,
        dt::Real,
        tolerance::Real=1e-13,
        max_iterations::Integer=100
    )
        @assert length(q₀) == length(p₀) == H.n_dof "Initial conditions must match Hamiltonian DOF ($(H.n_dof))"
        @assert tspan[2] > tspan[1] "End time must be greater than start time"
        @assert dt > 0 "Time step must be positive"
        @assert tolerance > 0 "Tolerance must be positive"
        @assert max_iterations > 0 "Max iterations must be positive"

        P = typeof(params)
        new{P}(
            H,
            (Float64(tspan[1]), Float64(tspan[2])),
            Vector{Float64}(q₀),
            Vector{Float64}(p₀),
            params,
            Float64(dt),
            Float64(tolerance),
            Int(max_iterations)
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
    IntegratorBuffers

Pre-allocated workspace for the symmetric projection integrator.
Internal type - not part of public API.
"""
mutable struct IntegratorBuffers
    d::Int
    A::Matrix{Float64}
    Z_current::Vector{Float64}
    Z_post_phi::Vector{Float64}
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
    buffers::IntegratorBuffers
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
- `tolerance`: Target tolerance
- `final_residual`: Residual at termination
- `step`: Integration step number where failure occurred
- `time`: Simulation time at failure
- `solver_name`: Name of the solver that failed
"""
struct NonlinearSolveError <: Exception
    iterations::Int
    tolerance::Float64
    final_residual::Float64
    step::Int
    time::Float64
    solver_name::String
end

function Base.showerror(io::IO, e::NonlinearSolveError)
    print(io, "NonlinearSolveError: $(e.solver_name) failed to converge at step $(e.step) (t=$(e.time)) ")
    print(io, "after $(e.iterations) iterations (residual=$(e.final_residual), tolerance=$(e.tolerance))")
end

# Deprecated alias for backward compatibility
const NewtonConvergenceError = NonlinearSolveError
