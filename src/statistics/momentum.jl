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

function compute_momentum_timeseries(solution::WeberSolution; stride::Int=1)::MomentumData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(solution.t) < 1
        throw(ArgumentError("solution must have at least 1 time point"))
    end

    # Extract problem data
    prob = solution.prob
    system = prob.system
    n_particles = system.n_particles
    dims = system.dims

    # Compute indices and allocate
    indices = 1:stride:length(solution.t)
    n_points = length(indices)
    t = solution.t[indices]

    # Pre-allocate arrays
    linear_momentum = [Vector{Float64}(undef, dims) for _ = 1:n_points]
    linear_momentum_components = Matrix{Float64}(undef, n_points, dims)
    linear_momentum_magnitude = Vector{Float64}(undef, n_points)

    # Angular momentum: dimension-dependent allocation
    if dims == 1
        angular_momentum = nothing
        angular_momentum_magnitude = nothing
    elseif dims == 2
        angular_momentum = Vector{Float64}(undef, n_points)
        angular_momentum_magnitude = angular_momentum  # Same array for 2D (Lz is scalar)
    else  # dims == 3
        angular_momentum = [Vector{Float64}(undef, 3) for _ = 1:n_points]
        angular_momentum_magnitude = Vector{Float64}(undef, n_points)
    end

    # Main computation loop
    @inbounds for (pt_idx, sol_idx) in enumerate(indices)
        q = solution.q[sol_idx]
        p = solution.p[sol_idx]

        # Linear momentum: P_total = sum_i p_i
        P = zeros(dims)
        for i = 1:n_particles
            for d = 1:dims
                coord_idx = (i - 1) * dims + d
                P[d] += p[coord_idx]
            end
        end
        linear_momentum[pt_idx] .= P
        for d = 1:dims
            linear_momentum_components[pt_idx, d] = P[d]
        end
        linear_momentum_magnitude[pt_idx] = sqrt(sum(P .^ 2))

        # Angular momentum: L_total = sum_i r_i x p_i
        if dims == 2
            # 2D: L_z = sum_i (x_i * p_y_i - y_i * p_x_i)
            L_z = 0.0
            for i = 1:n_particles
                x_idx = (i - 1) * dims + 1
                y_idx = (i - 1) * dims + 2
                L_z += q[x_idx] * p[y_idx] - q[y_idx] * p[x_idx]
            end
            angular_momentum[pt_idx] = L_z
        elseif dims == 3
            # 3D: L = sum_i r_i x p_i
            Lx, Ly, Lz = 0.0, 0.0, 0.0
            for i = 1:n_particles
                x_idx = (i - 1) * dims + 1
                y_idx = (i - 1) * dims + 2
                z_idx = (i - 1) * dims + 3
                # Cross product components: r x p
                Lx += q[y_idx] * p[z_idx] - q[z_idx] * p[y_idx]
                Ly += q[z_idx] * p[x_idx] - q[x_idx] * p[z_idx]
                Lz += q[x_idx] * p[y_idx] - q[y_idx] * p[x_idx]
            end
            angular_momentum[pt_idx] .= [Lx, Ly, Lz]
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
        n_particles,
        dims,
    )
end

function Base.show(io::IO, data::MomentumData)
    print(io, "MomentumData($(data.n_particles) particles, $(data.dims)D, $(length(data.t)) timesteps)")
end

function Base.show(io::IO, ::MIME"text/plain", data::MomentumData)
    println(io, "MomentumData")
    println(io, "  Particles: $(data.n_particles)")
    println(io, "  Dimensions: $(data.dims)")
    println(io, "  Timesteps: $(length(data.t))")
    println(io, "  t: $(data.t[1]) -> $(data.t[end])")
    println(io, "  |P| range: $(minimum(data.linear_momentum_magnitude)) -> $(maximum(data.linear_momentum_magnitude))")
    if !isnothing(data.angular_momentum_magnitude)
        println(io, "  |L| range: $(minimum(data.angular_momentum_magnitude)) -> $(maximum(data.angular_momentum_magnitude))")
    end
end
