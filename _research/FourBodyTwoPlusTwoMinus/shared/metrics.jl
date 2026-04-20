"""
Diagnostic utilities for 4-body 2+/2− Weber runs.
All functions take flat length-(N*d) vectors or a full `HamiltonianSolution`.
"""

using LinearAlgebra

function unflatten(v::AbstractVector, n::Int, d::Int)
    return reshape(copy(v), d, n)'  # returns n×d
end

"""
    pair_distances(q, n, d) -> Vector{Float64}
Returns the distances in the order (1,2),(1,3),...(n-1,n).
"""
function pair_distances(q::AbstractVector, n::Int, d::Int)
    X = unflatten(q, n, d)
    out = Float64[]
    for i = 1:n, j = (i+1):n
        push!(out, norm(view(X, i, :) .- view(X, j, :)))
    end
    return out
end

"""
    pair_rdots(q, p, masses, n, d) -> Vector{Float64}
Radial velocities ṙ_ij for each pair (order matches pair_distances).
"""
function pair_rdots(q, p, masses, n::Int, d::Int)
    X = unflatten(q, n, d)
    P = unflatten(p, n, d)
    V = P ./ masses
    out = Float64[]
    for i = 1:n, j = (i+1):n
        dx = view(X, i, :) .- view(X, j, :)
        dv = view(V, i, :) .- view(V, j, :)
        r = norm(dx)
        push!(out, (dx ⋅ dv) / r)
    end
    return out
end

"""
    critical_radii(charges, masses, c) -> Matrix{Float64}
Returns N×N matrix ρ_ij = q_i q_j / (μ_ij c²) for each pair; entries are
negative for unlike-charge pairs (no critical radius for attraction).
"""
function critical_radii(charges, masses, c)
    n = length(charges)
    rho = zeros(n, n)
    for i = 1:n, j = (i+1):n
        mu = (masses[i] * masses[j]) / (masses[i] + masses[j])
        rho[i, j] = charges[i] * charges[j] / (mu * c^2)
        rho[j, i] = rho[i, j]
    end
    return rho
end

"""
    total_linear_momentum(p, n, d) -> Vector
"""
total_linear_momentum(p, n, d) = vec(sum(unflatten(p, n, d), dims = 1))

"""
    total_angular_momentum(q, p, n, d) -> scalar (2D) or Vec3 (3D)
"""
function total_angular_momentum(q, p, n::Int, d::Int)
    X = unflatten(q, n, d)
    P = unflatten(p, n, d)
    if d == 2
        L = 0.0
        for i = 1:n
            L += X[i, 1] * P[i, 2] - X[i, 2] * P[i, 1]
        end
        return L
    elseif d == 3
        L = zeros(3)
        for i = 1:n
            L .+= cross(view(X, i, :), view(P, i, :))
        end
        return L
    else
        return 0.0
    end
end

"""
    bound_indicator(sol; escape_radius=20.0) -> Bool
Naive: sample sol every stride; return false if any particle's distance from COM
ever exceeds escape_radius, or retcode != :Success.
"""
function bound_indicator(sol; escape_radius::Float64 = 20.0, stride::Int = 10)
    sol.retcode == :Success || return false
    n = sol.prob.system.n_particles
    d = sol.prob.system.dims
    masses = sol.prob.masses
    for k = 1:stride:length(sol.t)
        X = unflatten(sol.q[k], n, d)
        com = (masses' * X) ./ sum(masses)
        for i = 1:n
            if norm(view(X, i, :) .- vec(com)) > escape_radius
                return false
            end
        end
    end
    return true
end

"""
    jacobi_coordinates_2p2m(q, masses) -> (R_com, r_plus, r_minus, r_inter)
Adapted to 2+/2− splitting: particles 1,2 are positive, 3,4 negative.
- r_plus  = x₂ − x₁  (inside the + dimer)
- r_minus = x₄ − x₃  (inside the − dimer)
- r_inter = center_minus − center_plus
"""
function jacobi_coordinates_2p2m(q::AbstractVector, masses::AbstractVector)
    X = unflatten(q, 4, length(q) ÷ 4)
    com_plus = (masses[1] * X[1, :] + masses[2] * X[2, :]) / (masses[1] + masses[2])
    com_minus = (masses[3] * X[3, :] + masses[4] * X[4, :]) / (masses[3] + masses[4])
    R_com = (sum(masses) > 0 ?
             (masses[1]*X[1,:] + masses[2]*X[2,:] + masses[3]*X[3,:] + masses[4]*X[4,:]) /
             sum(masses) : zeros(size(X, 2)))
    r_plus = X[2, :] .- X[1, :]
    r_minus = X[4, :] .- X[3, :]
    r_inter = com_minus .- com_plus
    return (R_com, r_plus, r_minus, r_inter)
end
