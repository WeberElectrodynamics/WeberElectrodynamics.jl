using LinearAlgebra

"""
Weber force timeseries between particle pairs.

# Fields
- `forces`: Dict mapping `(i,j)` to force vectors over time
- `t`: Time points
- `n_particles`, `dims`: System dimensions
"""
struct ForceData
    forces::Dict{Tuple{Int,Int}, Vector{Vector{Float64}}}
    t::Vector{Float64}
    n_particles::Int
    dims::Int
end

"""Pre-allocated buffers for force computation (internal)."""
struct ForceComputationBuffers
    r::Vector{Float64}
    v::Vector{Float64}
    a::Vector{Float64}
    r_hat::Vector{Float64}
end

function ForceComputationBuffers(dims::Int)
    ForceComputationBuffers(zeros(dims), zeros(dims), zeros(dims), zeros(dims))
end

"""
    compute_force_timeseries(sol, n_particles, dims, masses, charges, c; stride=1) -> ForceData

Compute Weber forces between all particle pairs from a solution.
"""
function compute_force_timeseries(sol::WeberSolution,
                                  n_particles::Int,
                                  dims::Int,
                                  masses::Vector{Float64},
                                  charges::Vector{Float64},
                                  c::Float64;
                                  stride::Int=1)::ForceData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(sol.t) < 2
        throw(ArgumentError("solution must have at least 2 time points for force computation"))
    end
    if length(masses) != n_particles
        throw(ArgumentError("masses length ($(length(masses))) must equal n_particles ($n_particles)"))
    end
    if length(charges) != n_particles
        throw(ArgumentError("charges length ($(length(charges))) must equal n_particles ($n_particles)"))
    end
    expected_dim = n_particles * dims
    actual_dim = length(sol.q[1])
    if actual_dim != expected_dim
        throw(ArgumentError("dimension mismatch: n_particles=$n_particles × dims=$dims = $expected_dim, but solution has dimension $actual_dim"))
    end
    if c <= 0
        throw(ArgumentError("speed of light c must be positive, got $c"))
    end

    indices = 1:stride:length(sol.t)
    n_steps = length(indices)

    if n_steps < 2
        throw(ArgumentError("stride=$stride results in only $n_steps points; need at least 2 for acceleration computation"))
    end

    dt = sol.t[indices[2]] - sol.t[indices[1]]

    t_forces = sol.t[indices[1:end-1]]
    n_force_steps = length(t_forces)

    velocities = Vector{Vector{Vector{Float64}}}(undef, n_particles)
    for particle in 1:n_particles
        m = masses[particle]

        vels = Vector{Vector{Float64}}(undef, n_steps)
        @inbounds for (i, idx) in enumerate(indices)
            vels[i] = Vector{Float64}(undef, dims)
            for d in 1:dims
                p_idx = (particle - 1) * dims + d
                vels[i][d] = sol.p[idx][p_idx] / m
            end
        end
        velocities[particle] = vels
    end

    accelerations = Vector{Vector{Vector{Float64}}}(undef, n_particles)
    for particle in 1:n_particles
        accels = Vector{Vector{Float64}}(undef, n_force_steps)
        @inbounds for t in 1:n_force_steps
            accels[t] = Vector{Float64}(undef, dims)
            @. accels[t] = (velocities[particle][t+1] - velocities[particle][t]) / dt
        end
        accelerations[particle] = accels
    end

    positions = Vector{Vector{Vector{Float64}}}(undef, n_particles)
    for particle in 1:n_particles
        pos = Vector{Vector{Float64}}(undef, n_steps)
        @inbounds for (i, idx) in enumerate(indices)
            pos[i] = Vector{Float64}(undef, dims)
            for d in 1:dims
                q_idx = (particle - 1) * dims + d
                pos[i][d] = sol.q[idx][q_idx]
            end
        end
        positions[particle] = pos
    end

    forces = Dict{Tuple{Int,Int}, Vector{Vector{Float64}}}()
    buf = ForceComputationBuffers(dims)

    for i in 1:n_particles
        for j in (i+1):n_particles
            qi = charges[i]
            qj = charges[j]

            force_series = Vector{Vector{Float64}}(undef, n_force_steps)

            @inbounds for t in 1:n_force_steps
                r_i = positions[i][t]
                r_j = positions[j][t]
                v_i = velocities[i][t]
                v_j = velocities[j][t]
                a_i = accelerations[i][t]
                a_j = accelerations[j][t]

                @. buf.r = r_i - r_j
                @. buf.v = v_i - v_j
                @. buf.a = a_i - a_j

                r_norm = norm(buf.r)
                @. buf.r_hat = buf.r / r_norm

                v_dot_v = dot(buf.v, buf.v)
                r_dot_a = dot(buf.r, buf.a)
                r_hat_dot_v = dot(buf.r_hat, buf.v)

                factor = (qi * qj / (r_norm^2)) * (1.0 + (1.0 / c^2) * (v_dot_v + r_dot_a - 1.5 * r_hat_dot_v^2))

                force_series[t] = Vector{Float64}(undef, dims)
                @. force_series[t] = factor * buf.r_hat
            end

            forces[(i, j)] = force_series

            neg_series = Vector{Vector{Float64}}(undef, n_force_steps)
            @inbounds for t in 1:n_force_steps
                neg_series[t] = Vector{Float64}(undef, dims)
                @. neg_series[t] = -force_series[t]
            end
            forces[(j, i)] = neg_series
        end
    end

    return ForceData(forces, t_forces, n_particles, dims)
end

"""
Newton's third law violation analysis for Weber forces.

# Fields
- `pair_violations`: `|F_ij + F_ji|` timeseries per pair
- `t`: Time points
- `max_violations`, `mean_violations`, `rms_violations`: Statistics per pair
- `global_max_violation`: Maximum across all pairs
- `n_pairs`: Number of particle pairs
"""
struct NewtonsThirdLawData
    pair_violations::Dict{Tuple{Int,Int}, Vector{Float64}}
    t::Vector{Float64}
    max_violations::Dict{Tuple{Int,Int}, Float64}
    mean_violations::Dict{Tuple{Int,Int}, Float64}
    rms_violations::Dict{Tuple{Int,Int}, Float64}
    global_max_violation::Float64
    n_pairs::Int
end

"""
    check_newtons_third_law(force_data) -> NewtonsThirdLawData

Check Newton's third law violations in computed Weber forces.
"""
function check_newtons_third_law(force_data::ForceData)::NewtonsThirdLawData
    n = force_data.n_particles
    n_times = length(force_data.t)

    pair_violations = Dict{Tuple{Int,Int}, Vector{Float64}}()
    max_violations = Dict{Tuple{Int,Int}, Float64}()
    mean_violations = Dict{Tuple{Int,Int}, Float64}()
    rms_violations = Dict{Tuple{Int,Int}, Float64}()

    global_max = 0.0
    n_pairs = 0

    for i in 1:n
        for j in (i+1):n
            n_pairs += 1

            F_ij = force_data.forces[(i, j)]
            F_ji = force_data.forces[(j, i)]

            violations = Vector{Float64}(undef, n_times)
            @inbounds for t in 1:n_times
                violations[t] = norm(F_ij[t] + F_ji[t])
            end

            pair_violations[(i, j)] = violations
            max_violations[(i, j)] = maximum(violations)
            mean_violations[(i, j)] = sum(violations) / length(violations)
            rms_violations[(i, j)] = sqrt(sum(x -> x^2, violations) / length(violations))

            if max_violations[(i, j)] > global_max
                global_max = max_violations[(i, j)]
            end
        end
    end

    return NewtonsThirdLawData(
        pair_violations,
        force_data.t,
        max_violations,
        mean_violations,
        rms_violations,
        global_max,
        n_pairs
    )
end
