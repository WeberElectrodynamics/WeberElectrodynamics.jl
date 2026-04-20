using LinearAlgebra

"""
    ForceStatistics

Summary statistics over a force magnitude timeseries.

# Fields
- `min::Float64`: Minimum force magnitude.
- `max::Float64`: Maximum force magnitude.
- `mean::Float64`: Mean force magnitude.
- `range::Float64`: `max − min`.
"""
struct ForceStatistics
    min::Float64
    max::Float64
    mean::Float64
    range::Float64
end

"""
    PhaseSpaceData

Pair (r, ṙ) phase-space timeseries extracted alongside force data.

# Fields
- `t::Vector{Float64}`: Time points.
- `separation_distance::Vector{Float64}`: Interparticle separation r(t).
- `radial_velocity::Vector{Float64}`: Radial velocity ṙ(t) = dr/dt.
- `theta::Union{Vector{Float64},Nothing}`: Polar angle θ for 2D/3D; `nothing` in 1D.
- `angular_momentum::Union{Vector{Float64},Nothing}`: Angular momentum magnitude;
  `nothing` in 1D.
"""
struct PhaseSpaceData
    t::Vector{Float64}
    separation_distance::Vector{Float64}
    radial_velocity::Vector{Float64}
    theta::Union{Vector{Float64},Nothing}
    angular_momentum::Union{Vector{Float64},Nothing}
end

"""
    PairForceData

Comprehensive Weber force timeseries for a single particle pair (i, j).

Forces are decomposed into a vector form (3 velocity/acceleration correction terms)
and a radial form (2 terms), both sharing the same κ-scaled Coulomb base.

# Fields
- `t::Vector{Float64}`: Time points.
- `dims::Int`: Spatial dimension.
- `pair::Tuple{Int,Int}`: Particle indices (i, j).
- `kappa::Float64`: Zöllner coupling factor κ_ij (1.0 = standard Weber).
- `charge_product::Float64`: k = qᵢqⱼ; sign determines repulsion/attraction.
- `force::Vector{Vector{Float64}}`: Total force vector at each time point.
- `magnitude::Vector{Float64}`: |F| at each time point.
- `stats::ForceStatistics`: Min/max/mean/range of |F|.
- `coulomb::Vector{Vector{Float64}}`: κ-scaled Coulomb base force vector.
- `vector_term_vv`, `vector_term_ra`, `vector_term_rv2`: Three correction terms
  in the vector form decomposition.
- `radial_term_rdot2`, `radial_term_rddot`: Two correction terms in the radial
  form decomposition.
- `zollner_extra_magnitude::Vector{Float64}`: |(κ−1)·F_Coulomb| (zero for
  standard Weber without Zöllner).
- `phase_space::PhaseSpaceData`: Pair phase-space portrait (r, ṙ, θ, L).
"""
struct PairForceData
    # Metadata
    t::Vector{Float64}
    dims::Int
    pair::Tuple{Int,Int}
    # Zöllner extension field (1.0 = standard Weber, no Zöllner)
    kappa::Float64           # coupling factor κ_ij
    charge_product::Float64  # k = q_i * q_j, determines repulsion/attraction sign

    # Total force (vector per timestep)
    force::Vector{Vector{Float64}}
    magnitude::Vector{Float64}
    stats::ForceStatistics

    # Shared Coulomb term (base for both decompositions); scaled by κ
    coulomb::Vector{Vector{Float64}}

    # Vector form decomposition (3 additional terms, Coulomb shared)
    # F = κ * coulomb_base * (1 + (v·v + r·a - 1.5*(r̂·v)²) / c²)
    vector_term_vv::Vector{Vector{Float64}}
    vector_term_ra::Vector{Vector{Float64}}
    vector_term_rv2::Vector{Vector{Float64}}

    # Radial form decomposition (2 additional terms, Coulomb shared)
    # F = κ * coulomb_base * (1 - ṙ²/(2c²) + r·r̈/c²)
    radial_term_rdot2::Vector{Vector{Float64}}
    radial_term_rddot::Vector{Vector{Float64}}

    # Zöllner extra force magnitude: |(κ-1) * F_coulomb_base| per timestep
    zollner_extra_magnitude::Vector{Float64}

    # Phase space data (computed at same timesteps, reuses r and ṙ)
    phase_space::PhaseSpaceData
end

"""
    compute_pair_force_timeseries(sol, pair, n_particles, dims, masses, charges, c; stride=1) -> PairForceData

Compute the Weber force decomposition for one particle pair over a simulation.

Both the vector form and the radial form of the Weber force are evaluated, together
with phase-space data (r, ṙ, θ, L).

# Arguments
- `sol::HamiltonianSolution`: Completed simulation.
- `pair::Tuple{Int,Int}`: Particle pair indices (order does not matter; i ≠ j).
- `n_particles::Int`, `dims::Int`: System geometry.
- `masses`, `charges::Vector{Float64}`: Physical parameters.
- `c::Float64`: Speed of light.

# Keywords
- `stride=1`: Downsample factor; every `stride`-th timestep is included.

# Returns
- `PairForceData` with force vectors, decompositions, statistics, and phase-space data.
"""
function compute_pair_force_timeseries(
    sol::HamiltonianSolution,
    pair::Tuple{Int,Int},
    n_particles::Int,
    dims::Int,
    masses::AbstractVector{Float64},
    charges::AbstractVector{Float64},
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

    # Phase space arrays (reuse computed r and ṙ)
    separation_distance = Vector{Float64}(undef, n_force_steps)
    radial_velocity_arr = Vector{Float64}(undef, n_force_steps)
    theta_arr = dims >= 2 ? Vector{Float64}(undef, n_force_steps) : nothing
    angular_momentum_arr = dims >= 2 ? Vector{Float64}(undef, n_force_steps) : nothing

    # Working buffers
    r_vec = zeros(dims)
    v_vec = zeros(dims)
    a_vec = zeros(dims)
    r_hat = zeros(dims)

    c2 = c * c
    k = charges[i] * charges[j]

    # Zöllner coupling for this pair (1.0 = standard Weber)
    kappa_ij = kappas(sol.prob)[_pair_index(i, j, n_particles)]

    zollner_extra_magnitude = Vector{Float64}(undef, n_force_steps)

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

        # Phase space data (reuse r and r_dot)
        separation_distance[t] = r
        radial_velocity_arr[t] = r_dot
        if !isnothing(theta_arr)
            theta_arr[t] = atan(r_vec[2], r_vec[1])
        end
        if !isnothing(angular_momentum_arr)
            if dims == 2
                angular_momentum_arr[t] = r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1]
            elseif dims == 3
                cx = r_vec[2] * v_vec[3] - r_vec[3] * v_vec[2]
                cy = r_vec[3] * v_vec[1] - r_vec[1] * v_vec[3]
                cz = r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1]
                angular_momentum_arr[t] = sqrt(cx^2 + cy^2 + cz^2)
            end
        end

        # Coulomb base vector, scaled by κ: (κ * k / r²) * r̂
        coulomb_coeff = kappa_ij * k / (r * r)
        @. coulomb[t] = coulomb_coeff * r_hat

        # Zöllner extra force magnitude: |(κ-1) * k / r²|
        zollner_extra_magnitude[t] = abs((kappa_ij - 1.0) * k / (r * r))

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

    phase_space = PhaseSpaceData(
        t_forces,
        separation_distance,
        radial_velocity_arr,
        theta_arr,
        angular_momentum_arr,
    )

    return PairForceData(
        t_forces,
        dims,
        pair,
        kappa_ij,
        k,
        force,
        magnitude,
        stats,
        coulomb,
        vector_term_vv,
        vector_term_ra,
        vector_term_rv2,
        radial_term_rdot2,
        radial_term_rddot,
        zollner_extra_magnitude,
        phase_space,
    )
end
