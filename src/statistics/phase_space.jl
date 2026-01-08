# Phase space data computation

using LinearAlgebra

"""
Phase space data for (r, rdot) portraits.

# Fields
- `t::Vector{Float64}`: Time points
- `r::Vector{Float64}`: Separation distance between particle pair
- `rdot::Vector{Float64}`: Radial velocity dr/dt
- `theta::Union{Vector{Float64}, Nothing}`: Angle in 2D (optional)
- `L::Union{Vector{Float64}, Nothing}`: Angular momentum (optional)
"""
struct PhaseSpaceData
    t::Vector{Float64}
    r::Vector{Float64}
    rdot::Vector{Float64}
    theta::Union{Vector{Float64}, Nothing}
    L::Union{Vector{Float64}, Nothing}
end

"""
    compute_phase_space_data(sol::IntegratorSolution, n_particles::Int, dims::Int,
                              masses::Vector{Float64};
                              particle_pair::Tuple{Int,Int}=(1,2),
                              stride::Int=1,
                              compute_angle::Bool=true,
                              compute_angular_momentum::Bool=true) -> PhaseSpaceData

Extract phase space coordinates (r, rdot) for a particle pair.

# Arguments
- `sol::WeberSolution`: Integration solution
- `n_particles::Int`: Number of particles
- `dims::Int`: Spatial dimensions (1, 2, or 3)
- `masses::Vector{Float64}`: Particle masses
- `particle_pair::Tuple{Int,Int}`: Which particle pair to analyze (default: (1,2))
- `stride::Int`: Downsampling factor
- `compute_angle::Bool`: Compute angle θ (2D only)
- `compute_angular_momentum::Bool`: Compute angular momentum L

# Returns
- `PhaseSpaceData`: Phase space data for the particle pair
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

    for (k, idx) in enumerate(indices)
        q = sol.q[idx]
        p = sol.p[idx]

        # Extract positions
        pos_i = [q[(i-1)*dims + d] for d in 1:dims]
        pos_j = [q[(j-1)*dims + d] for d in 1:dims]

        # Extract momenta and compute velocities
        vel_i = [p[(i-1)*dims + d] / m_i for d in 1:dims]
        vel_j = [p[(j-1)*dims + d] / m_j for d in 1:dims]

        # Relative position and velocity
        rel_pos = pos_i - pos_j
        rel_vel = vel_i - vel_j

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

        # Angular momentum (2D: L = x*vy - y*vx, 3D: |r × v|)
        if !isnothing(L_vals)
            if dims == 2
                L_vals[k] = rel_pos[1] * rel_vel[2] - rel_pos[2] * rel_vel[1]
            elseif dims == 3
                L_vec = cross(rel_pos, rel_vel)
                L_vals[k] = norm(L_vec)
            end
        end
    end

    return PhaseSpaceData(t, r_vals, rdot_vals, theta_vals, L_vals)
end
