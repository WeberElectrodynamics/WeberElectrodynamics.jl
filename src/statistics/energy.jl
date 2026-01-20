"""
Energy timeseries data from a solution.

# Fields
- `t`, `total`: Time points and total energy
- `kinetic`, `potential`: Optional component energies
- `max_local_error`: Max deviation from initial energy
- `relative_energy_range`: (E_max - E_min) / E_initial
"""
struct EnergyData
    t::Vector{Float64}
    total::Vector{Float64}
    kinetic::Union{Vector{Float64}, Nothing}
    potential::Union{Vector{Float64}, Nothing}
    max_local_error::Float64
    relative_energy_range::Float64
end

"""
    compute_energy_timeseries(sol, total_energy_func, [KE_func], [PE_func], params; stride=1) -> EnergyData

Compute energy timeseries from a solution. Energy functions have signature `(q, p, params, t)`.
"""
function compute_energy_timeseries(solution::WeberSolution,
                                   total_energy_func::Function,
                                   KE_func::Union{Function, Nothing}=nothing,
                                   PE_func::Union{Function, Nothing}=nothing,
                                   params=nothing;
                                   stride::Int=1)::EnergyData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(solution.t) < 1
        throw(ArgumentError("solution must have at least 1 time point"))
    end

    indices = 1:stride:length(solution.t)
    n_points = length(indices)
    t = solution.t[indices]

    total = zeros(n_points)
    for (i, idx) in enumerate(indices)
        total[i] = total_energy_func(solution.q[idx], solution.p[idx], params, solution.t[idx])
    end

    kinetic = if !isnothing(KE_func)
        ke = zeros(n_points)
        for (i, idx) in enumerate(indices)
            ke[i] = KE_func(solution.q[idx], solution.p[idx], params, solution.t[idx])
        end
        ke
    else
        nothing
    end

    potential = if !isnothing(PE_func)
        pe = zeros(n_points)
        for (i, idx) in enumerate(indices)
            pe[i] = PE_func(solution.q[idx], solution.p[idx], params, solution.t[idx])
        end
        pe
    else
        nothing
    end

    E_initial = total[1]
    max_local_error = maximum(abs.(total .- E_initial))

    E_max = maximum(total)
    E_min = minimum(total)

    relative_energy_range = if abs(E_initial) < 100 * eps(Float64)
        NaN
    else
        (E_max - E_min) / E_initial
    end

    return EnergyData(t, total, kinetic, potential, max_local_error, relative_energy_range)
end
