struct EnergyData
    t::Vector{Float64}
    total_energy::Vector{Float64}
    kinetic_energy::Union{Vector{Float64},Nothing}
    potential_energy::Union{Vector{Float64},Nothing}
    max_local_error::Float64
    relative_energy_range::Float64
end

function compute_energy_timeseries(
    solution::WeberSolution,
    total_energy_func::Function,
    kinetic_energy_func::Union{Function,Nothing} = nothing,
    potential_energy_func::Union{Function,Nothing} = nothing,
    params = nothing;
    stride::Int = 1,
)::EnergyData
    if stride <= 0
        throw(ArgumentError("stride must be positive, got $stride"))
    end
    if length(solution.t) < 1
        throw(ArgumentError("solution must have at least 1 time point"))
    end

    indices = 1:stride:length(solution.t)
    n_points = length(indices)
    t = solution.t[indices]

    total_energy = zeros(n_points)
    for (i, idx) in enumerate(indices)
        total_energy[i] =
            total_energy_func(solution.q[idx], solution.p[idx], params, solution.t[idx])
    end

    kinetic_energy = if !isnothing(kinetic_energy_func)
        ke = zeros(n_points)
        for (i, idx) in enumerate(indices)
            ke[i] = kinetic_energy_func(
                solution.q[idx],
                solution.p[idx],
                params,
                solution.t[idx],
            )
        end
        ke
    else
        nothing
    end

    potential_energy = if !isnothing(potential_energy_func)
        pe = zeros(n_points)
        for (i, idx) in enumerate(indices)
            pe[i] = potential_energy_func(
                solution.q[idx],
                solution.p[idx],
                params,
                solution.t[idx],
            )
        end
        pe
    else
        nothing
    end

    E_initial = total_energy[1]
    # Single-pass reduction without temporary allocation
    max_local_error = mapreduce(x -> abs(x - E_initial), max, total_energy)

    # Single pass for both min and max
    E_min, E_max = extrema(total_energy)

    relative_energy_range = if abs(E_initial) < 100 * eps(Float64)
        NaN
    else
        (E_max - E_min) / E_initial
    end

    return EnergyData(
        t,
        total_energy,
        kinetic_energy,
        potential_energy,
        max_local_error,
        relative_energy_range,
    )
end
