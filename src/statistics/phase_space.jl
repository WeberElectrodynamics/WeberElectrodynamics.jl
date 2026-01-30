using LinearAlgebra

struct PhaseSpaceData
    t::Vector{Float64}
    separation_distance::Vector{Float64}
    radial_velocity::Vector{Float64}
    theta::Union{Vector{Float64},Nothing}
    angular_momentum::Union{Vector{Float64},Nothing}
end

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
    separation_distance_vals = zeros(n_points)
    radial_velocity_vals = zeros(n_points)
    theta_vals = compute_angle && dims >= 2 ? zeros(n_points) : nothing
    angular_momentum_vals = compute_angular_momentum && dims >= 2 ? zeros(n_points) : nothing

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
            pos_i[d] = q[(i-1)*dims+d]
            pos_j[d] = q[(j-1)*dims+d]
            vel_i[d] = p[(i-1)*dims+d] / m_i
            vel_j[d] = p[(j-1)*dims+d] / m_j
        end

        # Relative position and velocity (in-place)
        @. rel_pos = pos_i - pos_j
        @. rel_vel = vel_i - vel_j

        # Separation distance
        r = norm(rel_pos)
        separation_distance_vals[k] = r

        # Radial velocity: rdot = (r . v) / |r|
        if r > eps(Float64)
            radial_velocity_vals[k] = dot(rel_pos, rel_vel) / r
        else
            radial_velocity_vals[k] = 0.0
        end

        # Angle (2D only)
        if !isnothing(theta_vals)
            theta_vals[k] = atan(rel_pos[2], rel_pos[1])
        end

        # Angular momentum (2D: L = x*vy - y*vx, 3D: compute norm directly)
        if !isnothing(angular_momentum_vals)
            if dims == 2
                angular_momentum_vals[k] = rel_pos[1] * rel_vel[2] - rel_pos[2] * rel_vel[1]
            elseif dims == 3
                # Compute cross product norm directly (avoids allocation)
                cx = rel_pos[2] * rel_vel[3] - rel_pos[3] * rel_vel[2]
                cy = rel_pos[3] * rel_vel[1] - rel_pos[1] * rel_vel[3]
                cz = rel_pos[1] * rel_vel[2] - rel_pos[2] * rel_vel[1]
                angular_momentum_vals[k] = sqrt(cx^2 + cy^2 + cz^2)
            end
        end
    end

    return PhaseSpaceData(t, separation_distance_vals, radial_velocity_vals, theta_vals, angular_momentum_vals)
end
