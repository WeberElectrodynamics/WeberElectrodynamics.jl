"""
    PairEnergyData

Energy decomposition timeseries for a single particle pair (i, j).

# Fields
- `pair::Tuple{Int,Int}`: Particle indices (i < j).
- `kappa::Float64`: Zöllner coupling factor κ_ij (1.0 for standard Weber).
- `coulomb_term::Vector{Float64}`: κ·qᵢqⱼ/r — Coulomb part of the potential.
- `velocity_term::Vector{Float64}`: −κ·qᵢqⱼ/r·ṙ²/(2c²) — velocity correction.
- `zollner_extra_potential::Vector{Float64}`: (κ−1)·qᵢqⱼ/r·(...) — extra
  potential from the Zöllner mismatch (zero when Zöllner is disabled).
- `total_pair_potential::Vector{Float64}`: Sum of all pair potential contributions.
- `radial_velocity::Vector{Float64}`: Radial velocity ṙ = dr/dt.
"""
struct PairEnergyData
    pair::Tuple{Int,Int}
    kappa::Float64
    coulomb_term::Vector{Float64}
    velocity_term::Vector{Float64}
    # Zöllner extension field (zero when Zöllner is disabled)
    zollner_extra_potential::Vector{Float64}
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
- `kinetic_energy::Vector{Float64}`: Total kinetic energy.
- `total_potential_energy::Vector{Float64}`: Sum of all pair Weber potentials.
- `total_zollner_residual::Vector{Float64}`: Summed Zöllner extra potential (zero
  when Zöllner extension is disabled).
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
    # Zöllner extension field (zero when Zöllner is disabled)
    total_zollner_residual::Vector{Float64}
    pair_energies::Dict{Tuple{Int,Int},PairEnergyData}
    statistics::EnergyStatistics
    hamiltonian_validation_error::Vector{Float64}
    n_particles::Int
    n_pairs::Int
end

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

function compute_pair_weber_components(
    q::AbstractVector{Float64},
    p::AbstractVector{Float64},
    i::Int,
    j::Int,
    masses::AbstractVector{Float64},
    charges::AbstractVector{Float64},
    c::Float64,
    dims::Int,
    kappa::Float64 = 1.0,
)::Tuple{Float64,Float64,Float64,Float64}
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

    r_dot_v = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        vi_d = p[qi_start+d] / mi
        vj_d = p[qj_start+d] / mj
        dv = vi_d - vj_d
        r_dot_v += dq * dv
    end
    rdot = r_dot_v / r

    # Base charge product
    @inbounds base_k = charges[i] * charges[j]

    # κ-scaled Coulomb term: κ·qᵢqⱼ/r
    coulomb_term = kappa * base_k / r

    # Velocity-dependent term: -κ·qᵢqⱼ/r * ṙ²/(2c²)
    c_squared = c^2
    velocity_term = -coulomb_term * rdot^2 / (2 * c_squared)

    # Zöllner extra potential: (κ-1)·qᵢqⱼ/r · (1 - ṙ²/(2c²))
    # This is the emergent 'gravitational' residual from the mismatch.
    zollner_extra = (kappa - 1.0) * base_k / r * (1.0 - rdot^2 / (2 * c_squared))

    return (coulomb_term, velocity_term, rdot, zollner_extra)
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

Evaluates kinetic energy, per-pair Weber potentials, and the compiled Hamiltonian
at each selected timestep. The `hamiltonian_validation_error` field cross-checks
the manual decomposition against the Symbolics-compiled function.

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
    κs = kappas(prob)
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

    # Initialize pair energy storage with pre-sized Dict
    pair_energies = sizehint!(Dict{Tuple{Int,Int},PairEnergyData}(), n_pairs)
    for i = 1:n
        for j = (i+1):n
            kappa_ij = κs[_pair_index(i, j, n)]
            pair_energies[(i, j)] = PairEnergyData(
                (i, j),
                kappa_ij,
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
                Vector{Float64}(undef, n_points),
            )
        end
    end

    total_zollner_residual = zeros(Float64, n_points)

    weber_decomp = get_term(system, :weber).pair_decomposition
    zollner_decomp = get_term(system, :zollner).pair_decomposition

    # Main computation loop
    @inbounds for (pt_idx, sol_idx) in enumerate(indices)
        q = solution.q[sol_idx]
        p = solution.p[sol_idx]
        t_pt = solution.t[sol_idx]

        # Kinetic energy
        KE = compute_total_kinetic_energy(p, ms, d)
        kinetic_energy[pt_idx] = KE

        # Pair-wise potential energies
        PE_total = 0.0
        zollner_sum = 0.0
        for i = 1:n
            for j = (i+1):n
                wc = weber_decomp(i, j, q, p, params_vec)
                zc = zollner_decomp(i, j, q, p, params_vec)
                coulomb = wc.coulomb
                velocity = wc.velocity
                rdot = wc.rdot
                zollner_extra = zc.zollner_extra
                pair_data = pair_energies[(i, j)]
                pair_data.coulomb_term[pt_idx] = coulomb
                pair_data.velocity_term[pt_idx] = velocity
                pair_data.zollner_extra_potential[pt_idx] = zollner_extra
                pair_data.total_pair_potential[pt_idx] = coulomb + velocity
                pair_data.radial_velocity[pt_idx] = rdot
                PE_total += coulomb + velocity
                zollner_sum += zollner_extra
            end
        end
        total_potential_energy[pt_idx] = PE_total
        total_energy[pt_idx] = KE + PE_total
        total_zollner_residual[pt_idx] = zollner_sum

        # Validate against compiled Hamiltonian
        H_compiled = hamiltonian_compiled(q, p, t_pt, params_vec)
        hamiltonian_validation[pt_idx] = abs(total_energy[pt_idx] - H_compiled)
    end

    # Compute statistics
    statistics = compute_energy_statistics(total_energy)

    return EnergyData(
        t,
        total_energy,
        kinetic_energy,
        total_potential_energy,
        total_zollner_residual,
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
