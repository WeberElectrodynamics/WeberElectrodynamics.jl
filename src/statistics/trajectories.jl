"""
    TrajectoryData

Per-particle position trajectories extracted from a `WeberSolution`.

# Fields
- `trajectories::Vector{Matrix{Float64}}`: One `(n_points, dims)` matrix per
  particle, where rows are time steps and columns are spatial dimensions.
- `initial_positions::Vector{Vector{Float64}}`: Position of each particle at
  the first selected time step.
- `final_positions::Vector{Vector{Float64}}`: Position of each particle at
  the last selected time step.
- `n_particles::Int`: Number of particles.
- `dims::Int`: Spatial dimension.
"""
struct TrajectoryData
    trajectories::Vector{Matrix{Float64}}
    initial_positions::Vector{Vector{Float64}}
    final_positions::Vector{Vector{Float64}}
    n_particles::Int
    dims::Int
end

"""
    compute_trajectory_data(sol, n_particles, dims; stride=1) -> TrajectoryData

Extract per-particle position trajectories from a `WeberSolution`.

# Arguments
- `sol::WeberSolution`: Completed simulation result.
- `n_particles::Int`: Number of particles encoded in `sol`.
- `dims::Int`: Spatial dimension (must match `sol.prob.system.dims`).

# Keywords
- `stride=1`: Downsample factor; every `stride`-th timestep is included.

# Returns
- `TrajectoryData` with one `(n_points, dims)` trajectory matrix per particle.
"""
function compute_trajectory_data(
    sol::WeberSolution,
    n_particles::Int,
    dims::Int;
    stride::Int = 1,
)::TrajectoryData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(sol.t) < 1
        throw(ArgumentError("solution must have at least 1 time point"))
    end
    expected_dim = n_particles * dims
    actual_dim = length(sol.q[1])
    if actual_dim != expected_dim
        throw(
            ArgumentError(
                "dimension mismatch: n_particles=$n_particles × dims=$dims = $expected_dim, but solution has dimension $actual_dim",
            ),
        )
    end

    indices = 1:stride:length(sol.t)
    n_points = length(indices)

    # Pre-allocate all trajectory matrices upfront (avoids per-particle allocations)
    trajectories = [Matrix{Float64}(undef, n_points, dims) for _ = 1:n_particles]
    initial_positions = [Vector{Float64}(undef, dims) for _ = 1:n_particles]
    final_positions = [Vector{Float64}(undef, dims) for _ = 1:n_particles]

    initial_idx = indices[1]
    final_idx = indices[end]

    for particle = 1:n_particles
        traj = trajectories[particle]
        @inbounds for (i, idx) in enumerate(indices)
            @inbounds for d = 1:dims
                coord_idx = (particle - 1) * dims + d
                traj[i, d] = sol.q[idx][coord_idx]
            end
        end

        # In-place assignment to pre-allocated vectors
        @inbounds for d = 1:dims
            coord_idx = (particle - 1) * dims + d
            initial_positions[particle][d] = sol.q[initial_idx][coord_idx]
            final_positions[particle][d] = sol.q[final_idx][coord_idx]
        end
    end

    return TrajectoryData(
        trajectories,
        initial_positions,
        final_positions,
        n_particles,
        dims,
    )
end
