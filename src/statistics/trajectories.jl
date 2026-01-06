struct TrajectoryData
    trajectories::Vector{Matrix{Float64}}
    initial_positions::Vector{Vector{Float64}}
    final_positions::Vector{Vector{Float64}}
    n_particles::Int
    dims::Int
end

function create_trajectory_data(sol::IntegratorSolution, n_particles::Int, dims::Int; stride::Int=1)::TrajectoryData
    indices = 1:stride:length(sol.t)
    n_points = length(indices)

    trajectories = Vector{Matrix{Float64}}(undef, n_particles)
    initial_positions = Vector{Vector{Float64}}(undef, n_particles)
    final_positions = Vector{Vector{Float64}}(undef, n_particles)

    for particle in 1:n_particles
        traj = zeros(n_points, dims)
        @inbounds for (i, idx) in enumerate(indices)
            for d in 1:dims
                coord_idx = (particle - 1) * dims + d
                traj[i, d] = sol.q[idx][coord_idx]
            end
        end
        trajectories[particle] = traj

        initial_positions[particle] = [sol.q[1][(particle - 1) * dims + d] for d in 1:dims]
        final_positions[particle] = [sol.q[end][(particle - 1) * dims + d] for d in 1:dims]
    end

    return TrajectoryData(trajectories, initial_positions, final_positions, n_particles, dims)
end
