# Phase space data computation

using LinearAlgebra

"""
Phase space data for (r, ṙ) portraits.

# Fields
- `t`, `r`, `rdot`: Time, separation distance, radial velocity
- `theta`: Angle in 2D (optional)
- `L`: Angular momentum (optional)
"""
struct PhaseSpaceData
    t::Vector{Float64}
    r::Vector{Float64}
    rdot::Vector{Float64}
    theta::Union{Vector{Float64}, Nothing}
    L::Union{Vector{Float64}, Nothing}
end

"""
    compute_phase_space_data(sol, n_particles, dims, masses; particle_pair=(1,2), stride=1, ...) -> PhaseSpaceData

Extract phase space coordinates (r, ṙ) for a particle pair from a `WeberSolution`.
"""
function compute_phase_space_data(sol::WeberSolution, n_particles::Int, dims::Int,
                                   masses::Vector{Float64};
                                   particle_pair::Tuple{Int,Int}=(1, 2),
                                   stride::Int=1,
                                   compute_angle::Bool=true,
                                   compute_angular_momentum::Bool=true)::PhaseSpaceData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end

    i, j = particle_pair
    if i < 1 || i > n_particles || j < 1 || j > n_particles
        throw(ArgumentError("particle indices must be in 1:$n_particles, got ($i, $j)"))
    end

    indices = 1:stride:length(sol.t)
    n_points = length(indices)

    t = sol.t[indices]
    r_vals = zeros(n_points)
    rdot_vals = zeros(n_points)
    theta_vals = compute_angle && dims >= 2 ? zeros(n_points) : nothing
    L_vals = compute_angular_momentum && dims >= 2 ? zeros(n_points) : nothing

    m_i, m_j = masses[i], masses[j]

    # Pre-allocate buffers outside loop (avoids per-iteration allocations)
    pos_i = zeros(dims)
    pos_j = zeros(dims)
    vel_i = zeros(dims)
    vel_j = zeros(dims)
    rel_pos = zeros(dims)
    rel_vel = zeros(dims)

    for (k, idx) in enumerate(indices)
        q = sol.q[idx]
        p = sol.p[idx]

        # Extract positions and velocities (in-place)
        @inbounds for d in 1:dims
            pos_i[d] = q[(i-1)*dims + d]
            pos_j[d] = q[(j-1)*dims + d]
            vel_i[d] = p[(i-1)*dims + d] / m_i
            vel_j[d] = p[(j-1)*dims + d] / m_j
        end

        # Relative position and velocity (in-place)
        @. rel_pos = pos_i - pos_j
        @. rel_vel = vel_i - vel_j

        # Separation distance
        r = norm(rel_pos)
        r_vals[k] = r

        # Radial velocity: rdot = (r . v) / |r|
        if r > eps(Float64)
            rdot_vals[k] = dot(rel_pos, rel_vel) / r
        else
            rdot_vals[k] = 0.0
        end

        # Angle (2D only)
        if !isnothing(theta_vals)
            theta_vals[k] = atan(rel_pos[2], rel_pos[1])
        end

        # Angular momentum (2D: L = x*vy - y*vx, 3D: compute norm directly)
        if !isnothing(L_vals)
            if dims == 2
                L_vals[k] = rel_pos[1] * rel_vel[2] - rel_pos[2] * rel_vel[1]
            elseif dims == 3
                # Compute cross product norm directly (avoids allocation)
                cx = rel_pos[2] * rel_vel[3] - rel_pos[3] * rel_vel[2]
                cy = rel_pos[3] * rel_vel[1] - rel_pos[1] * rel_vel[3]
                cz = rel_pos[1] * rel_vel[2] - rel_pos[2] * rel_vel[1]
                L_vals[k] = sqrt(cx^2 + cy^2 + cz^2)
            end
        end
    end

    return PhaseSpaceData(t, r_vals, rdot_vals, theta_vals, L_vals)
end
