using LinearAlgebra

struct ForceData
    forces::Dict{Tuple{Int,Int},Vector{Vector{Float64}}}
    t::Vector{Float64}
    n_particles::Int
    dims::Int
end

struct ForceComputationBuffers
    r::Vector{Float64}
    v::Vector{Float64}
    a::Vector{Float64}
    r_hat::Vector{Float64}
end

function ForceComputationBuffers(dims::Int)
    ForceComputationBuffers(zeros(dims), zeros(dims), zeros(dims), zeros(dims))
end

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

    # Use 3D arrays instead of nested vectors (avoids per-timestep allocations)
    velocities = Array{Float64}(undef, dims, n_steps, n_particles)
    @inbounds for particle in 1:n_particles
        m = masses[particle]
        for (i, idx) in enumerate(indices)
            for d in 1:dims
                p_idx = (particle - 1) * dims + d
                velocities[d, i, particle] = sol.p[idx][p_idx] / m
            end
        end
    end

    accelerations = Array{Float64}(undef, dims, n_force_steps, n_particles)
    @inbounds for particle in 1:n_particles
        for t in 1:n_force_steps
            for d in 1:dims
                accelerations[d, t, particle] = (velocities[d, t+1, particle] - velocities[d, t, particle]) / dt
            end
        end
    end

    positions = Array{Float64}(undef, dims, n_steps, n_particles)
    @inbounds for particle in 1:n_particles
        for (i, idx) in enumerate(indices)
            for d in 1:dims
                q_idx = (particle - 1) * dims + d
                positions[d, i, particle] = sol.q[idx][q_idx]
            end
        end
    end

    forces = Dict{Tuple{Int,Int},Vector{Vector{Float64}}}()
    buf = ForceComputationBuffers(dims)

    for i in 1:n_particles
        for j in (i+1):n_particles
            qi = charges[i]
            qj = charges[j]

            # Pre-allocate all force vectors before time loop
            force_series = [Vector{Float64}(undef, dims) for _ in 1:n_force_steps]
            neg_series = [Vector{Float64}(undef, dims) for _ in 1:n_force_steps]

            @inbounds for t in 1:n_force_steps
                # Use views into 3D arrays
                r_i = @view positions[:, t, i]
                r_j = @view positions[:, t, j]
                v_i = @view velocities[:, t, i]
                v_j = @view velocities[:, t, j]
                a_i = @view accelerations[:, t, i]
                a_j = @view accelerations[:, t, j]

                @. buf.r = r_i - r_j
                @. buf.v = v_i - v_j
                @. buf.a = a_i - a_j

                r_norm = norm(buf.r)
                @. buf.r_hat = buf.r / r_norm

                v_dot_v = dot(buf.v, buf.v)
                r_dot_a = dot(buf.r, buf.a)
                r_hat_dot_v = dot(buf.r_hat, buf.v)

                factor = (qi * qj / (r_norm^2)) * (1.0 + (1.0 / c^2) * (v_dot_v + r_dot_a - 1.5 * r_hat_dot_v^2))

                # In-place assignment to pre-allocated vectors
                @. force_series[t] = factor * buf.r_hat
                @. neg_series[t] = -force_series[t]
            end

            forces[(i, j)] = force_series
            forces[(j, i)] = neg_series
        end
    end

    return ForceData(forces, t_forces, n_particles, dims)
end

struct NewtonsThirdLawData
    pair_violations::Dict{Tuple{Int,Int},Vector{Float64}}
    t::Vector{Float64}
    max_violations::Dict{Tuple{Int,Int},Float64}
    mean_violations::Dict{Tuple{Int,Int},Float64}
    rms_violations::Dict{Tuple{Int,Int},Float64}
    global_max_violation::Float64
    n_pairs::Int
end

function check_newtons_third_law(force_data::ForceData)::NewtonsThirdLawData
    n = force_data.n_particles
    n_times = length(force_data.t)
    dims = force_data.dims

    pair_violations = Dict{Tuple{Int,Int},Vector{Float64}}()
    max_violations = Dict{Tuple{Int,Int},Float64}()
    mean_violations = Dict{Tuple{Int,Int},Float64}()
    rms_violations = Dict{Tuple{Int,Int},Float64}()

    global_max = 0.0
    n_pairs = 0

    # Pre-allocate buffer for sum computation (avoids per-iteration allocation)
    sum_buf = zeros(dims)

    for i in 1:n
        for j in (i+1):n
            n_pairs += 1

            F_ij = force_data.forces[(i, j)]
            F_ji = force_data.forces[(j, i)]

            violations = Vector{Float64}(undef, n_times)
            sum_val = 0.0
            sum_sq = 0.0
            max_val = 0.0

            @inbounds for t in 1:n_times
                # Compute norm without temporary allocation
                @. sum_buf = F_ij[t] + F_ji[t]
                v = norm(sum_buf)
                violations[t] = v
                sum_val += v
                sum_sq += v * v
                if v > max_val
                    max_val = v
                end
            end

            pair_violations[(i, j)] = violations
            max_violations[(i, j)] = max_val
            mean_violations[(i, j)] = sum_val / n_times
            rms_violations[(i, j)] = sqrt(sum_sq / n_times)

            if max_val > global_max
                global_max = max_val
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
