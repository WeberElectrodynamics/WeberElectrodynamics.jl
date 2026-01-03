struct TrajectoryData
    trajectories::Vector{Matrix{Float64}}
    initial_positions::Vector{Vector{Float64}}
    final_positions::Vector{Vector{Float64}}
    n_particles::Int
end

function create_trajectory_data(sol::IntegratorSolution, n_particles::Int; stride::Int=1)::TrajectoryData
    indices = 1:stride:length(sol.t)
    n_points = length(indices)

    trajectories = Vector{Matrix{Float64}}(undef, n_particles)
    initial_positions = Vector{Vector{Float64}}(undef, n_particles)
    final_positions = Vector{Vector{Float64}}(undef, n_particles)

    for particle in 1:n_particles
        # Assumes 2D system: each particle has (x, y) coordinates
        x_idx = (particle - 1) * 2 + 1
        y_idx = (particle - 1) * 2 + 2

        traj = zeros(n_points, 2)
        @inbounds for (i, idx) in enumerate(indices)
            traj[i, 1] = sol.q[idx][x_idx]
            traj[i, 2] = sol.q[idx][y_idx]
        end
        trajectories[particle] = traj

        initial_positions[particle] = [sol.q[1][x_idx], sol.q[1][y_idx]]
        final_positions[particle] = [sol.q[end][x_idx], sol.q[end][y_idx]]
    end

    return TrajectoryData(trajectories, initial_positions, final_positions, n_particles)
end
