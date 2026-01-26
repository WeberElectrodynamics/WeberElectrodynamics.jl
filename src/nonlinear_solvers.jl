using LinearAlgebra: norm

"""
    solve_relaxed_fixed_point!(u, f!, relaxation, convergence_tolerance, maximum_iterations; fu_buffer, u_old_buffer)

In-place relaxed fixed-point solve: `x_{n+1} = x_n - relaxation * f(x_n)`.

# Arguments
- `u`: Initial guess (modified in-place to contain solution)
- `f!`: In-place residual function `f!(result, u)`
- `relaxation`: Relaxation parameter in (0, 1]
- `convergence_tolerance`: Absolute tolerance for convergence
- `maximum_iterations`: Maximum iterations

# Returns
`(success::Bool, iterations::Int, residual::Float64)`
"""
function solve_relaxed_fixed_point!(
    u::AbstractVector{Float64},
    f!::F,
    relaxation::Float64,
    convergence_tolerance::Float64,
    maximum_iterations::Int;
    fu_buffer::AbstractVector{Float64},
    u_old_buffer::AbstractVector{Float64}
) where {F<:Function}
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
