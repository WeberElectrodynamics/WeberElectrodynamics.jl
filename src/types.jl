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
    SymmetricProjection <: WeberAlgorithm

Semi-explicit symplectic integrator using symmetric projection for
non-separable Hamiltonians. Second-order accurate, preserves symplectic
structure via extended phase space projection.

This is the default algorithm.
"""
struct SymmetricProjection <: WeberAlgorithm end

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

# Constructor
    WeberProblem(H::WeberHamiltonian, tspan, q₀, p₀; params, dt, tolerance=1e-13, max_iterations=100)

# Arguments
- `H`: Compiled Hamiltonian (from `@hamiltonian` macro)
- `tspan`: Time span as `(t_start, t_end)`
- `q₀`: Initial generalized coordinates
- `p₀`: Initial conjugate momenta
- `params`: Physical parameters (Vector or NamedTuple)
- `dt`: Time step size
- `tolerance`: Newton iteration convergence tolerance (default: 1e-13)
- `max_iterations`: Maximum Newton iterations per step (default: 100)

# Example
```julia
H = @hamiltonian (q, p, params) -> kinetic_energy(p, params) + potential_energy(q, params)
prob = WeberProblem(H, (0.0, 10.0), q₀, p₀; params=[m1, m2, k, c], dt=0.01)
sol = solve(prob)
```
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
- `t`: Time points
- `q`: Generalized coordinates at each time point
- `p`: Conjugate momenta at each time point
- `prob`: Reference to the problem
- `retcode`: Solution status (`:Success`, `:Failure`, `:MaxIters`)

# Indexing
- `sol[i]` returns `(t, q, p)` at step `i`
- `sol.t`, `sol.q`, `sol.p` for direct array access

# Iteration
```julia
for (t, q, p) in sol
    # process each time point
end
```
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

# Fields (user-accessible)
- `t`: Current time
- `q`: Current generalized coordinates
- `p`: Current conjugate momenta
- `step_count`: Number of steps taken

# Usage
```julia
integrator = init(prob, SymmetricProjection())
while integrator.t < t_end
    step!(integrator)
    # Access integrator.q, integrator.p, integrator.t
end
sol = solve!(integrator)  # Or continue stepping
```
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
    NewtonConvergenceError <: Exception

Thrown when Newton iteration fails to converge during integration step.
"""
struct NewtonConvergenceError <: Exception
    iterations::Int
    tolerance::Float64
    final_residual::Float64
    step::Int
    time::Float64
end

function Base.showerror(io::IO, e::NewtonConvergenceError)
    print(io, "NewtonConvergenceError: failed to converge at step $(e.step) (t=$(e.time)) ")
    print(io, "after $(e.iterations) iterations (residual=$(e.final_residual), tolerance=$(e.tolerance))")
end
