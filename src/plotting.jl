using Plots

function plot_trajectories(data::TrajectoryData)::Plots.Plot
    if data.dims == 1
        return _plot_trajectories_1d(data)
    elseif data.dims == 2
        return _plot_trajectories_2d(data)
    elseif data.dims == 3
        return _plot_trajectories_3d(data)
    else
        error("Unsupported dimension: $(data.dims). Only 1D, 2D, 3D supported.")
    end
end

function _plot_trajectories_1d(data::TrajectoryData)::Plots.Plot
    plt = plot(xlabel="Time Index", ylabel="x", legend=:topright)

    for particle in 1:data.n_particles
        plot!(plt, data.trajectories[particle][:, 1],
            label="Particle $particle", color=particle, linewidth=1.5)

        scatter!(plt, [1], [data.initial_positions[particle][1]],
            marker=:circle, markersize=8, color=particle,
            label=(particle == 1 ? "Initial" : ""))
        scatter!(plt, [size(data.trajectories[particle], 1)], [data.final_positions[particle][1]],
            marker=:square, markersize=8, color=particle,
            label=(particle == 1 ? "Final" : ""))
    end

    return plt
end

function _plot_trajectories_2d(data::TrajectoryData)::Plots.Plot
    plt = plot(xlabel="x", ylabel="y", aspect_ratio=:equal, legend=:topright)

    for particle in 1:data.n_particles
        plot!(plt, data.trajectories[particle][:, 1], data.trajectories[particle][:, 2],
            label="", color=particle, linewidth=1.5)

        scatter!(plt, [data.initial_positions[particle][1]], [data.initial_positions[particle][2]],
            marker=:circle, markersize=8, color=particle,
            label=(particle == 1 ? "Initial" : ""))
        scatter!(plt, [data.final_positions[particle][1]], [data.final_positions[particle][2]],
            marker=:square, markersize=8, color=particle,
            label=(particle == 1 ? "Final" : ""))
    end

    return plt
end

function _plot_trajectories_3d(data::TrajectoryData)::Plots.Plot
    plt = plot(xlabel="x", ylabel="y", zlabel="z", legend=:topright)

    for particle in 1:data.n_particles
        plot!(plt, data.trajectories[particle][:, 1],
              data.trajectories[particle][:, 2],
              data.trajectories[particle][:, 3],
            label="", color=particle, linewidth=1.5)

        scatter!(plt, [data.initial_positions[particle][1]],
                [data.initial_positions[particle][2]],
                [data.initial_positions[particle][3]],
            marker=:circle, markersize=8, color=particle,
            label=(particle == 1 ? "Initial" : ""))
        scatter!(plt, [data.final_positions[particle][1]],
                [data.final_positions[particle][2]],
                [data.final_positions[particle][3]],
            marker=:square, markersize=8, color=particle,
            label=(particle == 1 ? "Final" : ""))
    end

    return plt
end

function plot_energy(data::EnergyData)::Plots.Plot
    plt = plot(xlabel="Time", ylabel="Energy", legend=:topright)

    plot!(plt, data.t, data.total, label="Total", linewidth=2, color=:black)

    if !isnothing(data.kinetic)
        plot!(plt, data.t, data.kinetic, label="Kinetic", linewidth=1.5, color=:blue)
    end

    if !isnothing(data.potential)
        plot!(plt, data.t, data.potential, label="Potential", linewidth=1.5, color=:red)
    end

    return plt
end

function plot_newtons_third_law(data::NewtonsThirdLawData)::Plots.Plot
    plt = plot(xlabel="Time", ylabel="||F_ij + F_ji||", legend=:topright,
               title="Newton's 3rd Law (max: $(round(data.global_max_violation, sigdigits=3)))")

    for (pair, violations) in data.pair_violations
        i, j = pair
        plot!(plt, data.t, violations, label="($i,$j)", linewidth=1.5)
    end

    return plt
end
