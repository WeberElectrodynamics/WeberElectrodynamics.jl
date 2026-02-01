using LinearAlgebra

struct ForceStatistics
    min::Float64
    max::Float64
    mean::Float64
    range::Float64
end

struct PairForceData
    # Metadata
    t::Vector{Float64}
    dims::Int
    pair::Tuple{Int,Int}

    # Total force (vector per timestep)
    force::Vector{Vector{Float64}}
    magnitude::Vector{Float64}
    stats::ForceStatistics

    # Shared Coulomb term (base for both decompositions)
    coulomb::Vector{Vector{Float64}}

    # Vector form decomposition (3 additional terms, Coulomb shared)
    # F = coulomb * (1 + (v·v + r·a - 1.5*(r̂·v)²) / c²)
    vector_term_vv::Vector{Vector{Float64}}
    vector_term_ra::Vector{Vector{Float64}}
    vector_term_rv2::Vector{Vector{Float64}}

    # Radial form decomposition (2 additional terms, Coulomb shared)
    # F = coulomb * (1 - ṙ²/(2c²) + r·r̈/c²)
    radial_term_rdot2::Vector{Vector{Float64}}
    radial_term_rddot::Vector{Vector{Float64}}
end

function compute_pair_force_timeseries(
    sol::WeberSolution,
    pair::Tuple{Int,Int},
    n_particles::Int,
    dims::Int,
    masses::Vector{Float64},
    charges::Vector{Float64},
    c::Float64;
    stride::Int = 1,
)::PairForceData
    i, j = pair

    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(sol.t) < 2
        throw(
            ArgumentError(
                "solution must have at least 2 time points for force computation",
            ),
        )
    end
    if i < 1 || i > n_particles || j < 1 || j > n_particles
        throw(ArgumentError("pair indices must be in range 1:$n_particles, got ($i, $j)"))
    end
    if i == j
        throw(ArgumentError("pair indices must be different, got ($i, $j)"))
    end
    if length(masses) != n_particles
        throw(
            ArgumentError(
                "masses length ($(length(masses))) must equal n_particles ($n_particles)",
            ),
        )
    end
    if length(charges) != n_particles
        throw(
            ArgumentError(
                "charges length ($(length(charges))) must equal n_particles ($n_particles)",
            ),
        )
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
    if c <= 0
        throw(ArgumentError("speed of light c must be positive, got $c"))
    end

    indices = 1:stride:length(sol.t)
    n_steps = length(indices)

    if n_steps < 2
        throw(
            ArgumentError(
                "stride=$stride results in only $n_steps points; need at least 2 for acceleration computation",
            ),
        )
    end

    dt = sol.t[indices[2]] - sol.t[indices[1]]
    t_forces = sol.t[indices[1:(end-1)]]
    n_force_steps = length(t_forces)

    # Extract positions and velocities for the pair
    mi = masses[i]
    mj = masses[j]
    qi_start = (i - 1) * dims
    qj_start = (j - 1) * dims

    # Build position and velocity arrays for both particles
    positions_i = Array{Float64}(undef, dims, n_steps)
    positions_j = Array{Float64}(undef, dims, n_steps)
    velocities_i = Array{Float64}(undef, dims, n_steps)
    velocities_j = Array{Float64}(undef, dims, n_steps)

    @inbounds for (step_idx, sol_idx) in enumerate(indices)
        for d = 1:dims
            positions_i[d, step_idx] = sol.q[sol_idx][qi_start+d]
            positions_j[d, step_idx] = sol.q[sol_idx][qj_start+d]
            velocities_i[d, step_idx] = sol.p[sol_idx][qi_start+d] / mi
            velocities_j[d, step_idx] = sol.p[sol_idx][qj_start+d] / mj
        end
    end

    # Compute accelerations via finite difference
    accelerations_i = Array{Float64}(undef, dims, n_force_steps)
    accelerations_j = Array{Float64}(undef, dims, n_force_steps)

    @inbounds for t = 1:n_force_steps
        for d = 1:dims
            accelerations_i[d, t] = (velocities_i[d, t+1] - velocities_i[d, t]) / dt
            accelerations_j[d, t] = (velocities_j[d, t+1] - velocities_j[d, t]) / dt
        end
    end

    # Pre-allocate output arrays
    force = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]
    magnitude = Vector{Float64}(undef, n_force_steps)
    coulomb = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]
    vector_term_vv = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]
    vector_term_ra = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]
    vector_term_rv2 = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]
    radial_term_rdot2 = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]
    radial_term_rddot = [Vector{Float64}(undef, dims) for _ = 1:n_force_steps]

    # Working buffers
    r_vec = zeros(dims)
    v_vec = zeros(dims)
    a_vec = zeros(dims)
    r_hat = zeros(dims)

    c2 = c * c
    k = charges[i] * charges[j]

    sum_mag = 0.0
    min_mag = Inf
    max_mag = -Inf

    @inbounds for t = 1:n_force_steps
        # Relative quantities
        for d = 1:dims
            r_vec[d] = positions_i[d, t] - positions_j[d, t]
            v_vec[d] = velocities_i[d, t] - velocities_j[d, t]
            a_vec[d] = accelerations_i[d, t] - accelerations_j[d, t]
        end

        r = norm(r_vec)
        @. r_hat = r_vec / r

        # Dot products
        v_dot_v = dot(v_vec, v_vec)
        r_dot_a = dot(r_vec, a_vec)
        rhat_dot_v = dot(r_hat, v_vec)

        # Radial quantities
        r_dot = dot(r_vec, v_vec) / r
        r_ddot = (r_dot_a + v_dot_v - rhat_dot_v^2) / r

        # Coulomb base vector: (k / r²) * r̂
        coulomb_coeff = k / (r * r)
        @. coulomb[t] = coulomb_coeff * r_hat

        # Vector form terms
        # term_vv = coulomb * (v·v / c²)
        vv_coeff = v_dot_v / c2
        @. vector_term_vv[t] = coulomb[t] * vv_coeff

        # term_ra = coulomb * (r·a / c²)
        ra_coeff = r_dot_a / c2
        @. vector_term_ra[t] = coulomb[t] * ra_coeff

        # term_rv2 = coulomb * (-1.5 * (r̂·v)² / c²)
        rv2_coeff = -1.5 * rhat_dot_v^2 / c2
        @. vector_term_rv2[t] = coulomb[t] * rv2_coeff

        # Radial form terms
        # term_rdot2 = coulomb * (-ṙ² / (2c²))
        rdot2_coeff = -r_dot^2 / (2 * c2)
        @. radial_term_rdot2[t] = coulomb[t] * rdot2_coeff

        # term_rddot = coulomb * (r * r̈ / c²)
        rddot_coeff = r * r_ddot / c2
        @. radial_term_rddot[t] = coulomb[t] * rddot_coeff

        # Total force (using vector form)
        @. force[t] = coulomb[t] + vector_term_vv[t] + vector_term_ra[t] + vector_term_rv2[t]

        # Magnitude
        mag = norm(force[t])
        magnitude[t] = mag
        sum_mag += mag
        if mag < min_mag
            min_mag = mag
        end
        if mag > max_mag
            max_mag = mag
        end
    end

    mean_mag = sum_mag / n_force_steps
    stats = ForceStatistics(min_mag, max_mag, mean_mag, max_mag - min_mag)

    return PairForceData(
        t_forces,
        dims,
        pair,
        force,
        magnitude,
        stats,
        coulomb,
        vector_term_vv,
        vector_term_ra,
        vector_term_rv2,
        radial_term_rdot2,
        radial_term_rddot,
    )
end
