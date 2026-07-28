"""
    PairEnergyData

Energy decomposition timeseries for a single particle pair (i, j).

# Fields
- `pair::Tuple{Int,Int}`: Particle indices (i < j).
- `coulomb_term::Vector{Float64}`: qᵢqⱼ/r — Coulomb part of the potential.
- `velocity_term::Vector{Float64}`: −qᵢqⱼ/r·ṙ²/(2c²) — velocity correction.
- `total_pair_potential::Vector{Float64}`: Sum of all pair potential contributions.
- `radial_velocity::Vector{Float64}`: **Physical** radial velocity ṙ = dr/dt,
  obtained from the canonical momenta through the coupled pair solve.
"""
struct PairEnergyData
    pair::Tuple{Int,Int}
    coulomb_term::Vector{Float64}
    velocity_term::Vector{Float64}
    total_pair_potential::Vector{Float64}
    radial_velocity::Vector{Float64}
end

"""
    EnergyStatistics

Summary statistics for energy conservation over a simulation.

Local error is the step-to-step absolute change |E_t − E_{t−1}|.
Global error ratio is E_t / E_initial (ideal value: 1.0).
Global error percentage is |E_t − E_initial| / |E_initial| × 100.

# Fields
- `local_error_max`, `local_error_min`, `local_error_avg::Float64`:
  Maximum, minimum, and mean step-to-step energy change.
- `global_error_ratio_max`, `global_error_ratio_min`, `global_error_ratio_avg::Float64`:
  Maximum, minimum, and mean of E_t / E_initial.
- `global_error_percent_max`, `global_error_percent_min`, `global_error_percent_avg::Float64`:
  Maximum, minimum, and mean global error as a percentage.
"""
struct EnergyStatistics
    # Local error: |E_t - E_{t-1}| (absolute values)
    local_error_max::Float64
    local_error_min::Float64
    local_error_avg::Float64

    # Global error ratio: E_t / E_initial (~1.0 means good conservation)
    global_error_ratio_max::Float64
    global_error_ratio_min::Float64
    global_error_ratio_avg::Float64

    # Global error percentage: |(E_t - E_initial) / E_initial| * 100
    global_error_percent_max::Float64
    global_error_percent_min::Float64
    global_error_percent_avg::Float64
end

"""
    EnergyData

Complete energy timeseries and conservation statistics for an n-body simulation.

# Fields
- `t::Vector{Float64}`: Time points.
- `total_energy::Vector{Float64}`: Total Hamiltonian H = KE + V at each time point.
- `kinetic_energy::Vector{Float64}`: Total **physical** kinetic energy
  `Σ ½ mᵢ|vᵢ|²`, computed from the physical velocities recovered from the
  canonical momenta. It is not `Σ |pᵢ|²/(2mᵢ)`, which is only the Coulomb-limit
  value.
- `total_potential_energy::Vector{Float64}`: Sum of all pair Weber potentials.
- `pair_energies::Dict{Tuple{Int,Int},PairEnergyData}`: Per-pair energy decomposition.
- `statistics::EnergyStatistics`: Energy conservation quality metrics.
- `hamiltonian_validation_error::Vector{Float64}`: |H_manual − H_compiled| at each
  time point, cross-checking the manual decomposition against the compiled Hamiltonian.
- `n_particles::Int`, `n_pairs::Int`: System size.
"""
struct EnergyData
    t::Vector{Float64}
    total_energy::Vector{Float64}
    kinetic_energy::Vector{Float64}
    total_potential_energy::Vector{Float64}
    pair_energies::Dict{Tuple{Int,Int},PairEnergyData}
    statistics::EnergyStatistics
    hamiltonian_validation_error::Vector{Float64}
    n_particles::Int
    n_pairs::Int
end

"""
    compute_total_kinetic_energy(v, masses, dims) -> Float64

Physical kinetic energy `Σ ½ mᵢ|vᵢ|²` from **physical velocities** `v`.

`v` must come from [`physical_velocities`](@ref); passing canonical momenta
divided by masses is wrong for any velocity-dependent Hamiltonian.
"""
function compute_total_kinetic_energy(
    v::AbstractVector{Float64},
    masses::AbstractVector{Float64},
    dims::Int,
)::Float64
    KE = 0.0
    @inbounds for i in eachindex(masses)
        v_start = (i - 1) * dims + 1
        v_squared = 0.0
        @inbounds for d = 0:(dims-1)
            v_squared += v[v_start+d]^2
        end
        KE += 0.5 * masses[i] * v_squared
    end
    return KE
end

"""
    compute_pair_weber_components(q, v, i, j, charges, c, dims)
        -> (coulomb_term, velocity_term, rdot)

Weber pair potential split for pair `(i, j)` from positions `q` and **physical
velocities** `v`:

`U_ij = qᵢqⱼ/r · (1 − ṙ²/(2c²))`, returned as `(qᵢqⱼ/r, −qᵢqⱼṙ²/(2c²r), ṙ)`.
"""
function compute_pair_weber_components(
    q::AbstractVector{Float64},
    v::AbstractVector{Float64},
    i::Int,
    j::Int,
    charges::AbstractVector{Float64},
    c::Float64,
    dims::Int,
)::Tuple{Float64,Float64,Float64}
    qi_start = (i - 1) * dims + 1
    qj_start = (j - 1) * dims + 1

    # Compute separation distance
    r_squared = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        r_squared += dq^2
    end
    r = sqrt(r_squared)

    r_dot_v = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        r_dot_v += dq * (v[qi_start+d] - v[qj_start+d])
    end
    rdot = r_dot_v / r

    # Coulomb term: qᵢqⱼ/r
    @inbounds coulomb_term = charges[i] * charges[j] / r

    # Velocity-dependent term: -qᵢqⱼ/r * ṙ²/(2c²)
    velocity_term = -coulomb_term * rdot^2 / (2 * c^2)

    return (coulomb_term, velocity_term, rdot)
end

function compute_energy_statistics(total_energy::Vector{Float64})::EnergyStatistics
    n = length(total_energy)
    E_initial = total_energy[1]

    if n < 2
        return EnergyStatistics(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
    end

    # Local error: |E_t - E_{t-1}|
    local_error_max = 0.0
    local_error_min = Inf
    local_error_sum = 0.0
    @inbounds for i = 2:n
        err = abs(total_energy[i] - total_energy[i-1])
        local_error_max = max(local_error_max, err)
        local_error_min = min(local_error_min, err)
        local_error_sum += err
    end
    local_error_avg = local_error_sum / (n - 1)

    # Global error ratio: E_t / E_initial
    # Handle near-zero initial energy
    if abs(E_initial) < 100 * eps(Float64)
        return EnergyStatistics(
            local_error_max,
            local_error_min,
            local_error_avg,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
        )
    end

    global_ratio_max = -Inf
    global_ratio_min = Inf
    global_ratio_sum = 0.0
    global_percent_max = 0.0
    global_percent_min = Inf
    global_percent_sum = 0.0

    @inbounds for i = 2:n
        ratio = total_energy[i] / E_initial
        global_ratio_max = max(global_ratio_max, ratio)
        global_ratio_min = min(global_ratio_min, ratio)
        global_ratio_sum += ratio

        percent = abs((total_energy[i] - E_initial) / E_initial) * 100.0
        global_percent_max = max(global_percent_max, percent)
        global_percent_min = min(global_percent_min, percent)
        global_percent_sum += percent
    end

    global_ratio_avg = global_ratio_sum / (n - 1)
    global_percent_avg = global_percent_sum / (n - 1)

    return EnergyStatistics(
        local_error_max,
        local_error_min,
        local_error_avg,
        global_ratio_max,
        global_ratio_min,
        global_ratio_avg,
        global_percent_max,
        global_percent_min,
        global_percent_avg,
    )
end

"""
    compute_energy_timeseries(sol; stride=1) -> EnergyData

Compute the full energy decomposition timeseries from a `HamiltonianSolution`.

Evaluates physical kinetic energy, per-pair Weber potentials, and the compiled
Hamiltonian at each selected timestep.

The decomposition is the velocity-space energy
`E = Σ ½ mᵢ|vᵢ|² + Σ qᵢqⱼ/rᵢⱼ (1 − ṙᵢⱼ²/(2c²))`, built from the physical
velocities recovered from the canonical momenta. Because that equals the exact
canonical Hamiltonian by the Legendre transform, `hamiltonian_validation_error`
is a genuine independent cross-check of the compiled `H(q, p)` rather than a
restatement of it.

# Keywords
- `stride=1`: Downsample factor; every `stride`-th timestep is included.

# Returns
- `EnergyData` with all energy components and conservation statistics.
"""
function compute_energy_timeseries(
    solution::HamiltonianSolution;
    stride::Int = 1,
)::EnergyData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(solution.t) < 1
        throw(ArgumentError("solution must have at least 1 time point"))
    end

    # Extract problem data
    prob = solution.prob
    system = prob.system
    n = n_particles(prob)
    d = dims(prob)
    ms = masses(prob)
    params_vec = params(prob)
    hamiltonian_compiled = system.hamiltonian_compiled

    # Compute indices and allocate
    indices = 1:stride:length(solution.t)
    n_points = length(indices)
    t = solution.t[indices]

    total_energy = Vector{Float64}(undef, n_points)
    kinetic_energy = Vector{Float64}(undef, n_points)
    total_potential_energy = Vector{Float64}(undef, n_points)
    hamiltonian_validation = Vector{Float64}(undef, n_points)

    # Compute n_pairs without allocating pairs vector
    n_pairs = n * (n - 1) ÷ 2

    weber_term_or_nothing = has_term(system, :weber) ? get_term(system, :weber) : nothing
    weber_decomp =
        weber_term_or_nothing === nothing ? nothing :
        weber_term_or_nothing.pair_decomposition
    weber_kinetic =
        weber_term_or_nothing === nothing ? nothing : weber_term_or_nothing.kinetic_energy
    has_weber_decomposition = weber_decomp !== nothing
    pairs = pair_indices(system)

    # Initialize pair energy storage only when the system supplies the built-in
    # Weber pair decomposition. Generic Hamiltonians still get a valid
    # compiled-Hamiltonian energy timeseries, but arbitrary NamedTerm shapes
    # cannot be losslessly projected into PairEnergyData.
    pair_energies = sizehint!(Dict{Tuple{Int,Int},PairEnergyData}(), n_pairs)
    if has_weber_decomposition
        for i = 1:n
            for j = (i+1):n
                pair_energies[(i, j)] = PairEnergyData(
                    (i, j),
                    Vector{Float64}(undef, n_points),
                    Vector{Float64}(undef, n_points),
                    Vector{Float64}(undef, n_points),
                    Vector{Float64}(undef, n_points),
                )
            end
        end
    end

    # Main computation loop
    @inbounds for (pt_idx, sol_idx) in enumerate(indices)
        q = solution.q[sol_idx]
        p = solution.p[sol_idx]
        t_pt = solution.t[sol_idx]
        H_compiled = hamiltonian_compiled(q, p, t_pt, params_vec)

        if has_weber_decomposition
            # Physical kinetic energy from the recovered physical velocities.
            KE = weber_kinetic(q, p, params_vec)
            kinetic_energy[pt_idx] = KE

            # Pair-wise potential energies, one coupled solve per timestep.
            wc = weber_decomp(q, p, params_vec)
            PE_total = 0.0
            for (a, (i, j)) in enumerate(pairs)
                coulomb = wc.coulomb[a]
                velocity = wc.velocity[a]
                pair_data = pair_energies[(i, j)]
                pair_data.coulomb_term[pt_idx] = coulomb
                pair_data.velocity_term[pt_idx] = velocity
                pair_data.total_pair_potential[pt_idx] = coulomb + velocity
                pair_data.radial_velocity[pt_idx] = wc.rdot[a]
                PE_total += coulomb + velocity
            end
            total_potential_energy[pt_idx] = PE_total
            total_energy[pt_idx] = KE + PE_total
        else
            # Generic custom Hamiltonians have no physical-velocity map, so the
            # canonical kinetic split is the only meaningful one available.
            KE = 0.0
            for i = 1:n
                p_start = (i - 1) * d + 1
                acc = 0.0
                for dd = 0:(d-1)
                    acc += p[p_start+dd]^2
                end
                KE += acc / (2 * ms[i])
            end
            kinetic_energy[pt_idx] = KE
            total_energy[pt_idx] = H_compiled
            total_potential_energy[pt_idx] = H_compiled - KE
        end

        # Validate against compiled Hamiltonian
        hamiltonian_validation[pt_idx] = abs(total_energy[pt_idx] - H_compiled)
    end

    # Compute statistics
    statistics = compute_energy_statistics(total_energy)

    return EnergyData(
        t,
        total_energy,
        kinetic_energy,
        total_potential_energy,
        pair_energies,
        statistics,
        hamiltonian_validation,
        n,
        n_pairs,
    )
end

function Base.show(io::IO, data::EnergyData)
    print(
        io,
        "EnergyData($(data.n_particles) particles, $(data.n_pairs) pairs, $(length(data.t)) timesteps)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", data::EnergyData)
    println(io, "EnergyData")
    println(io, "  Particles: $(data.n_particles)")
    println(io, "  Pairs: $(data.n_pairs)")
    println(io, "  Timesteps: $(length(data.t))")
    println(io, "  t: $(data.t[1]) → $(data.t[end])")
    println(io, "  Statistics:")
    s = data.statistics
    println(
        io,
        "    Local error  (max/min/avg): $(s.local_error_max) / $(s.local_error_min) / $(s.local_error_avg)",
    )
    println(
        io,
        "    Global ratio (max/min/avg): $(s.global_error_ratio_max) / $(s.global_error_ratio_min) / $(s.global_error_ratio_avg)",
    )
    println(
        io,
        "    Global %     (max/min/avg): $(s.global_error_percent_max) / $(s.global_error_percent_min) / $(s.global_error_percent_avg)",
    )
    println(
        io,
        "  Hamiltonian validation: max error = $(maximum(data.hamiltonian_validation_error))",
    )
end

function _summary_min_pair_distance(sol::HamiltonianSolution, stride::Int)
    n = n_particles(sol.prob)
    d = dims(sol.prob)
    n < 2 && return Inf

    min_r = Inf
    @inbounds for sol_idx = 1:stride:length(sol.t)
        q = sol.q[sol_idx]
        for i = 1:n
            i0 = (i - 1) * d
            for j = (i+1):n
                j0 = (j - 1) * d
                r2 = 0.0
                for k = 1:d
                    dx = q[i0+k] - q[j0+k]
                    r2 += dx * dx
                end
                min_r = min(min_r, sqrt(r2))
            end
        end
    end
    return min_r
end

function _summary_angular_drift(momentum)
    if isnothing(momentum.angular_momentum_magnitude)
        return (max = nothing, relative_max = nothing)
    end

    initial = momentum.angular_momentum_magnitude[1]
    max_drift = 0.0
    @inbounds for value in momentum.angular_momentum_magnitude
        max_drift = max(max_drift, abs(value - initial))
    end
    relative = abs(initial) <= 100 * eps(Float64) ? NaN : max_drift / abs(initial)
    return (max = max_drift, relative_max = relative)
end

"""
    conservation_summary(sol; stride=1) -> NamedTuple

Return a compact conservation summary for notebooks, smoke tests, and
regression checks.

The result contains energy error statistics, maximum linear and angular
momentum drift, the minimum sampled pair distance, and the solution's
regularization diagnostics.
"""
function conservation_summary(sol::HamiltonianSolution; stride::Int = 1)
    stride > 0 || throw(ArgumentError("stride must be positive, got $stride"))

    energy = compute_energy_timeseries(sol; stride = stride)
    momentum = compute_momentum_timeseries(sol; stride = stride)

    linear_initial = momentum.linear_momentum_magnitude[1]
    linear_drift_max = 0.0
    @inbounds for value in momentum.linear_momentum_magnitude
        linear_drift_max = max(linear_drift_max, abs(value - linear_initial))
    end
    linear_relative =
        abs(linear_initial) <= 100 * eps(Float64) ? NaN :
        linear_drift_max / abs(linear_initial)

    angular = _summary_angular_drift(momentum)
    validation_max = maximum(energy.hamiltonian_validation_error)

    return (
        n_points = length(energy.t),
        energy = (
            local_error_max = energy.statistics.local_error_max,
            global_error_percent_max = energy.statistics.global_error_percent_max,
            hamiltonian_validation_max = validation_max,
        ),
        momentum = (
            linear_drift_max = linear_drift_max,
            linear_relative_drift_max = linear_relative,
            angular_drift_max = angular.max,
            angular_relative_drift_max = angular.relative_max,
        ),
        min_pair_distance = _summary_min_pair_distance(sol, stride),
        regularization = sol.regularization,
    )
end
