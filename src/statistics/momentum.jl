"""
    MomentumData

Total linear and angular momentum timeseries from a simulation.

Momentum conservation (constant |P| and |L|) is the primary quality indicator
for simulations without external forces.

# Fields
- `t::Vector{Float64}`: Time points.
- `linear_momentum::Vector{Vector{Float64}}`: Total momentum vector P(t).
- `linear_momentum_components::Matrix{Float64}`: `(n_points, dims)` matrix of
  individual momentum components.
- `linear_momentum_magnitude::Vector{Float64}`: |P(t)|.
- `angular_momentum`: `nothing` (1D), scalar `Vector{Float64}` (2D Lz values),
  or `Vector{Vector{Float64}}` (3D, full L vector per timestep).
- `angular_momentum_magnitude::Union{Nothing,Vector{Float64}}`: |L(t)|;
  `nothing` in 1D.
- `n_particles::Int`, `dims::Int`: System geometry.
"""
struct MomentumData
    t::Vector{Float64}
    linear_momentum::Vector{Vector{Float64}}
    linear_momentum_components::Matrix{Float64}
    linear_momentum_magnitude::Vector{Float64}
    angular_momentum::Union{Nothing,Vector{Float64},Vector{Vector{Float64}}}
    angular_momentum_magnitude::Union{Nothing,Vector{Float64}}
    n_particles::Int
    dims::Int
end

"""
    compute_momentum_timeseries(solution; stride=1) -> MomentumData

Compute total linear and angular momentum timeseries from a `HamiltonianSolution`.

# Keywords
- `stride=1`: Downsample factor; every `stride`-th timestep is included.

# Returns
- `MomentumData` with linear and angular momentum at each selected timestep.
"""
function compute_momentum_timeseries(
    solution::HamiltonianSolution;
    stride::Int = 1,
)::MomentumData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(solution.t) < 1
        throw(ArgumentError("solution must have at least 1 time point"))
    end

    # Extract problem data
    prob = solution.prob
    n = n_particles(prob)
    d = dims(prob)

    # Compute indices and allocate
    indices = 1:stride:length(solution.t)
    n_points = length(indices)
    t = solution.t[indices]

    # Pre-allocate arrays
    linear_momentum = [Vector{Float64}(undef, d) for _ = 1:n_points]
    linear_momentum_components = Matrix{Float64}(undef, n_points, d)
    linear_momentum_magnitude = Vector{Float64}(undef, n_points)

    # Angular momentum: dimension-dependent allocation
    if d == 1
        angular_momentum = nothing
        angular_momentum_magnitude = nothing
    elseif d == 2
        angular_momentum = Vector{Float64}(undef, n_points)
        angular_momentum_magnitude = Vector{Float64}(undef, n_points)
    else  # d == 3
        angular_momentum = [Vector{Float64}(undef, 3) for _ = 1:n_points]
        angular_momentum_magnitude = Vector{Float64}(undef, n_points)
    end
    P = Vector{Float64}(undef, d)

    # Main computation loop
    @inbounds for (pt_idx, sol_idx) in enumerate(indices)
        q = solution.q[sol_idx]
        p = solution.p[sol_idx]

        # Linear momentum: P_total = sum_i p_i
        fill!(P, 0.0)
        for i = 1:n
            for k = 1:d
                coord_idx = (i - 1) * d + k
                P[k] += p[coord_idx]
            end
        end
        P_mag_sq = 0.0
        for k = 1:d
            P_k = P[k]
            linear_momentum[pt_idx][k] = P_k
            linear_momentum_components[pt_idx, k] = P_k
            P_mag_sq += P_k * P_k
        end
        linear_momentum_magnitude[pt_idx] = sqrt(P_mag_sq)

        # Angular momentum: L_total = sum_i r_i x p_i
        if d == 2
            # 2D: L_z = sum_i (x_i * p_y_i - y_i * p_x_i)
            L_z = 0.0
            for i = 1:n
                x_idx = (i - 1) * d + 1
                y_idx = (i - 1) * d + 2
                L_z += q[x_idx] * p[y_idx] - q[y_idx] * p[x_idx]
            end
            angular_momentum[pt_idx] = L_z
            angular_momentum_magnitude[pt_idx] = abs(L_z)
        elseif d == 3
            # 3D: L = sum_i r_i x p_i
            Lx, Ly, Lz = 0.0, 0.0, 0.0
            for i = 1:n
                x_idx = (i - 1) * d + 1
                y_idx = (i - 1) * d + 2
                z_idx = (i - 1) * d + 3
                # Cross product components: r x p
                Lx += q[y_idx] * p[z_idx] - q[z_idx] * p[y_idx]
                Ly += q[z_idx] * p[x_idx] - q[x_idx] * p[z_idx]
                Lz += q[x_idx] * p[y_idx] - q[y_idx] * p[x_idx]
            end
            angular_momentum[pt_idx][1] = Lx
            angular_momentum[pt_idx][2] = Ly
            angular_momentum[pt_idx][3] = Lz
            angular_momentum_magnitude[pt_idx] = sqrt(Lx^2 + Ly^2 + Lz^2)
        end
    end

    return MomentumData(
        t,
        linear_momentum,
        linear_momentum_components,
        linear_momentum_magnitude,
        angular_momentum,
        angular_momentum_magnitude,
        n,
        d,
    )
end

function Base.show(io::IO, data::MomentumData)
    print(
        io,
        "MomentumData($(data.n_particles) particles, $(data.dims)D, $(length(data.t)) timesteps)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", data::MomentumData)
    println(io, "MomentumData")
    println(io, "  Particles: $(data.n_particles)")
    println(io, "  Dimensions: $(data.dims)")
    println(io, "  Timesteps: $(length(data.t))")
    println(io, "  t: $(data.t[1]) -> $(data.t[end])")
    println(
        io,
        "  |P| range: $(minimum(data.linear_momentum_magnitude)) -> $(maximum(data.linear_momentum_magnitude))",
    )
    if !isnothing(data.angular_momentum_magnitude)
        println(
            io,
            "  |L| range: $(minimum(data.angular_momentum_magnitude)) -> $(maximum(data.angular_momentum_magnitude))",
        )
    end
end
