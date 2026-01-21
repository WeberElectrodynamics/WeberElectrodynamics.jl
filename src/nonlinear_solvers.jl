using SciMLBase: AbstractNonlinearAlgorithm, NonlinearProblem, ReturnCode
using SimpleNonlinearSolve: AbstractSimpleNonlinearSolveAlgorithm, solve as nl_solve
using LinearAlgebra: norm

# =============================================================================
# RelaxedFixedPointSolver
# =============================================================================

"""
    RelaxedFixedPointSolver <: AbstractNonlinearAlgorithm

Relaxed fixed-point iteration: `x_{n+1} = x_n - relaxation * f(x_n)`.
Default solver for `SymmetricProjectionIntegrator`. Constructor: `RelaxedFixedPointSolver(; relaxation=0.25)`.
"""
struct RelaxedFixedPointSolver <: AbstractNonlinearAlgorithm
    relaxation::Float64

    function RelaxedFixedPointSolver(; relaxation::Real=0.25)
        @assert 0 < relaxation <= 1 "Relaxation must be in (0, 1], got $relaxation"
        new(Float64(relaxation))
    end
end

function Base.show(io::IO, alg::RelaxedFixedPointSolver)
    print(io, "RelaxedFixedPointSolver(relaxation=$(alg.relaxation))")
end

# =============================================================================
# In-place Nonlinear Solve Interface
# =============================================================================

"""
    solve_nonlinear!(u, f!, alg, convergence_tolerance, maximum_iterations; fu_buffer, u_old_buffer)

In-place nonlinear solve for Weber integrator.

# Arguments
- `u`: Initial guess (modified in-place to contain solution)
- `f!`: In-place residual function `f!(result, u)`
- `alg`: Nonlinear solver algorithm
- `convergence_tolerance`: Absolute tolerance for convergence
- `maximum_iterations`: Maximum iterations

# Keyword Arguments
- `fu_buffer`: Pre-allocated buffer for residual (same size as `u`)
- `u_old_buffer`: Pre-allocated buffer for previous iterate (same size as `u`)

# Returns
`(success::Bool, iterations::Int, residual::Float64)`
"""
function solve_nonlinear! end

# Dispatch for RelaxedFixedPointSolver (direct implementation, no allocations)
function solve_nonlinear!(
    u::AbstractVector{Float64},
    f!::F,
    alg::RelaxedFixedPointSolver,
    convergence_tolerance::Float64,
    maximum_iterations::Int;
    fu_buffer::AbstractVector{Float64},
    u_old_buffer::AbstractVector{Float64}
) where {F<:Function}
    relaxation = alg.relaxation

    for iter in 1:maximum_iterations
        f!(fu_buffer, u)

        copyto!(u_old_buffer, u)
        @. u = u - relaxation * fu_buffer

        # Convergence based on step size
        diff_norm_sq = zero(Float64)
        @inbounds for i in eachindex(u, u_old_buffer)
            diff_norm_sq += (u[i] - u_old_buffer[i])^2
        end
        diff_norm = sqrt(diff_norm_sq)

        if diff_norm < convergence_tolerance
            # Verify constraint satisfaction (not just iteration stationarity)
            f!(fu_buffer, u)
            residual = norm(fu_buffer)
            if residual < convergence_tolerance
                return (true, iter, residual)
            end
            # Continue iterating if residual still large
        end
    end

    # Failed to converge - compute final residual
    f!(fu_buffer, u)
    final_residual = norm(fu_buffer)
    return (false, maximum_iterations, final_residual)
end

# Dispatch for SimpleNonlinearSolve algorithms
# Note: Only derivative-free solvers (SimpleBroyden, SimpleDFSane, SimpleKlement, etc.)
# are supported. Newton-type solvers requiring automatic differentiation (SimpleNewtonRaphson,
# SimpleTrustRegion) will fail because the internal projection function uses pre-allocated
# Float64 buffers incompatible with ForwardDiff's dual numbers.
function solve_nonlinear!(
    u::AbstractVector{Float64},
    f!::Function,
    alg::AbstractSimpleNonlinearSolveAlgorithm,
    convergence_tolerance::Float64,
    maximum_iterations::Int;
    fu_buffer::AbstractVector{Float64},
    _u_old_buffer::Union{AbstractVector{Float64}, Nothing}=nothing  # unused, avoids allocation
)
    # Wrap f! to match SciML signature f!(result, u, p)
    wrapped_f! = (result, u_inner, _p) -> f!(result, u_inner)

    # Create NonlinearProblem with wrapped in-place function
    prob = NonlinearProblem{true}(wrapped_f!, copy(u), nothing)

    sol = nl_solve(prob, alg; abstol=convergence_tolerance, maxiters=maximum_iterations)

    # Copy solution back to u
    copyto!(u, sol.u)

    # Compute residual
    f!(fu_buffer, u)
    residual = norm(fu_buffer)

    success = sol.retcode == ReturnCode.Success
    iterations = maximum_iterations  # SimpleNonlinearSolve doesn't always expose iteration count

    return (success, iterations, residual)
end
