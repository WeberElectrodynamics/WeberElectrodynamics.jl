using LinearAlgebra

struct ForceData
    forces::Dict{Tuple{Int,Int}, Vector{Vector{Float64}}}
    t::Vector{Float64}
    n_particles::Int
end

struct ForceComputationBuffers
    r::Vector{Float64}
    v::Vector{Float64}
    a::Vector{Float64}
    r_hat::Vector{Float64}
end

function ForceComputationBuffers()
    ForceComputationBuffers(zeros(2), zeros(2), zeros(2), zeros(2))
end

function compute_force_timeseries(sol::IntegratorSolution,
                                  n_particles::Int,
                                  masses::Vector{Float64},
                                  charges::Vector{Float64},
                                  c::Float64;
                                  stride::Int=1)::ForceData
    @assert length(masses) == n_particles
    @assert length(charges) == n_particles

    indices = 1:stride:length(sol.t)
    n_steps = length(indices)

    dt = sol.t[indices[2]] - sol.t[indices[1]]

    t_forces = sol.t[indices[1:end-1]]
    n_force_steps = length(t_forces)

    velocities = Vector{Vector{Vector{Float64}}}(undef, n_particles)
    for particle in 1:n_particles
        x_idx = (particle - 1) * 2 + 1
        y_idx = (particle - 1) * 2 + 2
        m = masses[particle]

        vels = Vector{Vector{Float64}}(undef, n_steps)
        @inbounds for (i, idx) in enumerate(indices)
            px = sol.p[idx][x_idx]
            py = sol.p[idx][y_idx]
            vels[i] = Vector{Float64}(undef, 2)
            vels[i][1] = px / m
            vels[i][2] = py / m
        end
        velocities[particle] = vels
    end

    accelerations = Vector{Vector{Vector{Float64}}}(undef, n_particles)
    for particle in 1:n_particles
        accels = Vector{Vector{Float64}}(undef, n_force_steps)
        @inbounds for t in 1:n_force_steps
            accels[t] = Vector{Float64}(undef, 2)
            @. accels[t] = (velocities[particle][t+1] - velocities[particle][t]) / dt
        end
        accelerations[particle] = accels
    end

    positions = Vector{Vector{Vector{Float64}}}(undef, n_particles)
    for particle in 1:n_particles
        x_idx = (particle - 1) * 2 + 1
        y_idx = (particle - 1) * 2 + 2

        pos = Vector{Vector{Float64}}(undef, n_steps)
        @inbounds for (i, idx) in enumerate(indices)
            pos[i] = Vector{Float64}(undef, 2)
            pos[i][1] = sol.q[idx][x_idx]
            pos[i][2] = sol.q[idx][y_idx]
        end
        positions[particle] = pos
    end

    forces = Dict{Tuple{Int,Int}, Vector{Vector{Float64}}}()
    buf = ForceComputationBuffers()

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

                force_series[t] = Vector{Float64}(undef, 2)
                @. force_series[t] = factor * buf.r_hat
            end

            forces[(i, j)] = force_series

            neg_series = Vector{Vector{Float64}}(undef, n_force_steps)
            @inbounds for t in 1:n_force_steps
                neg_series[t] = Vector{Float64}(undef, 2)
                @. neg_series[t] = -force_series[t]
            end
            forces[(j, i)] = neg_series
        end
    end

    return ForceData(forces, t_forces, n_particles)
end

struct NewtonsThirdLawData
    pair_violations::Dict{Tuple{Int,Int}, Vector{Float64}}
    t::Vector{Float64}
    max_violations::Dict{Tuple{Int,Int}, Float64}
    mean_violations::Dict{Tuple{Int,Int}, Float64}
    rms_violations::Dict{Tuple{Int,Int}, Float64}
    global_max_violation::Float64
    n_pairs::Int
end

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
