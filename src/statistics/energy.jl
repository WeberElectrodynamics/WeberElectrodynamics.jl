"""
Energy analysis for Weber electrodynamics simulations.

Provides pair-wise energy decomposition and comprehensive error statistics.
"""

struct PairEnergyData
    pair::Tuple{Int,Int}
    coulomb_term::Vector{Float64}
    velocity_term::Vector{Float64}
    total_pair_potential::Vector{Float64}
    radial_velocity::Vector{Float64}
end

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
    compute_total_kinetic_energy(p, masses, dims) -> Float64

Compute total kinetic energy: Σᵢ pᵢ²/(2mᵢ)
"""
function compute_total_kinetic_energy(
    p::AbstractVector{Float64},
    masses::AbstractVector{Float64},
    dims::Int,
)::Float64
    KE = 0.0
    @inbounds for i in eachindex(masses)
        p_start = (i - 1) * dims + 1
        p_squared = 0.0
        @inbounds for d = 0:(dims-1)
            p_squared += p[p_start+d]^2
        end
        KE += p_squared / (2 * masses[i])
    end
    return KE
end

"""
    compute_pair_weber_components(q, p, i, j, masses, charges, c, dims) -> (coulomb, velocity, rdot)

Compute Weber potential components for particle pair (i, j).

Returns:
- coulomb_term: qᵢqⱼ/r
- velocity_term: -qᵢqⱼ/r * ṙ²/(2c²)
- rdot: radial velocity ṙ = (r⃗·v⃗)/r
"""
function compute_pair_weber_components(
    q::AbstractVector{Float64},
    p::AbstractVector{Float64},
    i::Int,
    j::Int,
    masses::AbstractVector{Float64},
    charges::AbstractVector{Float64},
    c::Float64,
    dims::Int,
)::Tuple{Float64,Float64,Float64}
    qi_start = (i - 1) * dims + 1
    qj_start = (j - 1) * dims + 1

    # Cache mass lookups
    @inbounds mi = masses[i]
    @inbounds mj = masses[j]

    # Compute separation distance
    r_squared = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        r_squared += dq^2
    end
    r = sqrt(r_squared)

    # Compute r⃗·v⃗ where v⃗ = vᵢ - vⱼ
    r_dot_v = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        vi_d = p[qi_start+d] / mi
        vj_d = p[qj_start+d] / mj
        dv = vi_d - vj_d
        r_dot_v += dq * dv
    end
    rdot = r_dot_v / r

    # Coulomb term: qᵢqⱼ/r
    @inbounds k = charges[i] * charges[j]
    coulomb_term = k / r

    # Velocity-dependent term: -qᵢqⱼ/r * ṙ²/(2c²)
    c_squared = c^2
    velocity_term = -coulomb_term * rdot^2 / (2 * c_squared)

    return (coulomb_term, velocity_term, rdot)
end

"""
    compute_energy_statistics(total_energy) -> EnergyStatistics

Compute local and global error statistics from energy timeseries.
"""
function compute_energy_statistics(total_energy::Vector{Float64})::EnergyStatistics
    n = length(total_energy)
    E_initial = total_energy[1]

    if n < 2
        return EnergyStatistics(
            0.0, 0.0, 0.0,
            1.0, 1.0, 1.0,
            0.0, 0.0, 0.0,
        )
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
            local_error_max, local_error_min, local_error_avg,
            NaN, NaN, NaN,
            NaN, NaN, NaN,
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
        local_error_max, local_error_min, local_error_avg,
        global_ratio_max, global_ratio_min, global_ratio_avg,
        global_percent_max, global_percent_min, global_percent_avg,
    )
end

"""
    compute_energy_timeseries(solution; stride=1) -> EnergyData

Compute comprehensive energy analysis for a Weber simulation.

Automatically extracts particle count and parameters from the solution.
Computes pair-wise energy decomposition and validates against the compiled Hamiltonian.

# Arguments
- `solution::WeberSolution`: Solution from `solve(prob)`
- `stride::Int=1`: Sample every `stride`-th timestep

# Returns
`EnergyData` containing:
- Time and energy timeseries
- Pair-wise energy decomposition (Coulomb, velocity terms, radial velocity)
- Error statistics (local and global)
- Validation error vs compiled Hamiltonian
"""
function compute_energy_timeseries(solution::WeberSolution; stride::Int = 1)::EnergyData
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
    masses = prob.masses
    charges = prob.charges
    c = prob.c
    params = prob.params
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
    n_pairs = n_particles * (n_particles - 1) ÷ 2

    # Initialize pair energy storage with pre-sized Dict
    pair_energies = sizehint!(Dict{Tuple{Int,Int},PairEnergyData}(), n_pairs)
    for i = 1:n_particles
        for j = (i+1):n_particles
            pair_energies[(i, j)] = PairEnergyData(
                (i, j),
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
            )
        end
    end

    # Main computation loop
    @inbounds for (pt_idx, sol_idx) in enumerate(indices)
        q = solution.q[sol_idx]
        p = solution.p[sol_idx]

        # Kinetic energy
        KE = compute_total_kinetic_energy(p, masses, dims)
        kinetic_energy[pt_idx] = KE

        # Pair-wise potential energies
        PE_total = 0.0
        for i = 1:n_particles
            for j = (i+1):n_particles
                coulomb, velocity, rdot =
                    compute_pair_weber_components(q, p, i, j, masses, charges, c, dims)
                pair_data = pair_energies[(i, j)]
                pair_data.coulomb_term[pt_idx] = coulomb
                pair_data.velocity_term[pt_idx] = velocity
                pair_data.total_pair_potential[pt_idx] = coulomb + velocity
                pair_data.radial_velocity[pt_idx] = rdot
                PE_total += coulomb + velocity
            end
        end
        total_potential_energy[pt_idx] = PE_total
        total_energy[pt_idx] = KE + PE_total

        # Validate against compiled Hamiltonian
        H_compiled = hamiltonian_compiled(q, p, params)
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
        n_particles,
        n_pairs,
    )
end

function Base.show(io::IO, data::EnergyData)
    print(io, "EnergyData($(data.n_particles) particles, $(data.n_pairs) pairs, $(length(data.t)) timesteps)")
end

function Base.show(io::IO, ::MIME"text/plain", data::EnergyData)
    println(io, "EnergyData")
    println(io, "  Particles: $(data.n_particles)")
    println(io, "  Pairs: $(data.n_pairs)")
    println(io, "  Timesteps: $(length(data.t))")
    println(io, "  t: $(data.t[1]) → $(data.t[end])")
    println(io, "  Statistics:")
    s = data.statistics
    println(io, "    Local error  (max/min/avg): $(s.local_error_max) / $(s.local_error_min) / $(s.local_error_avg)")
    println(io, "    Global ratio (max/min/avg): $(s.global_error_ratio_max) / $(s.global_error_ratio_min) / $(s.global_error_ratio_avg)")
    println(io, "    Global %     (max/min/avg): $(s.global_error_percent_max) / $(s.global_error_percent_min) / $(s.global_error_percent_avg)")
    println(io, "  Hamiltonian validation: max error = $(maximum(data.hamiltonian_validation_error))")
end
