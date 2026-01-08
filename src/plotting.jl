using Plots
using LinearAlgebra

const PLOT_DEFAULTS = (
    framestyle = :box,
    grid = true,
    gridalpha = 0.3,
    tickfontsize = 10,
    guidefontsize = 12,
    legendfontsize = 9,
    dpi = 150,
    margin = 5Plots.mm
)

# Public API

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

function plot_energy(data::EnergyData)::Plots.Plot
    E0 = data.total[1]
    relative_error = abs.((data.total .- E0) ./ E0)
    relative_error = max.(relative_error, eps(Float64))

    p1 = plot(; xlabel="", ylabel="Energy", legend=:topright,
              size=(600, 300), PLOT_DEFAULTS...)
    plot!(p1, data.t, data.total, label="Total", linewidth=2, color=:black)
    if !isnothing(data.kinetic)
        plot!(p1, data.t, data.kinetic, label="Kinetic", linewidth=1.5, color=:steelblue)
    end
    if !isnothing(data.potential)
        plot!(p1, data.t, data.potential, label="Potential", linewidth=1.5, color=:firebrick)
    end

    max_err = maximum(relative_error)
    p2 = plot(; xlabel="Time", ylabel="|ΔE/E₀|", legend=:topright,
              yscale=:log10, size=(600, 200), PLOT_DEFAULTS...)
    plot!(p2, data.t, relative_error, label="", linewidth=1.5, color=:black)
    hline!(p2, [max_err], linestyle=:dash, color=:gray, label="max = $(round(max_err, sigdigits=2))")

    plt = plot(p1, p2, layout=grid(2, 1, heights=[0.65, 0.35]), size=(600, 450))
    return plt
end

function plot_forces(data::ForceData)::Plots.Plot
    n_times = length(data.t)
    n = data.n_particles

    plt = plot(; xlabel="Time", ylabel="|F|",
               legend=:topright,
               size=(600, 400), palette=:tab10, PLOT_DEFAULTS...)

    for i in 1:n
        for j in (i+1):n
            F_ij = data.forces[(i, j)]

            magnitudes = Vector{Float64}(undef, n_times)
            for t in 1:n_times
                magnitudes[t] = norm(F_ij[t])
            end

            plot!(plt, data.t, magnitudes, label="F($i,$j)", linewidth=1.5)
        end
    end

    return plt
end

function plot_phase_space(data::PhaseSpaceData)::Plots.Plot
    plt = plot(; xlabel="r", ylabel="dr/dt",
               size=(500, 500), PLOT_DEFAULTS...)

    plot!(plt, data.r, data.rdot, label="", linewidth=1.5, color=:black)

    scatter!(plt, [data.r[1]], [data.rdot[1]],
        marker=:circle, markersize=6, color=:steelblue, label="Initial")
    scatter!(plt, [data.r[end]], [data.rdot[end]],
        marker=:square, markersize=6, color=:firebrick, label="Final")

    return plt
end

# Internal helpers

function _plot_trajectories_1d(data::TrajectoryData)::Plots.Plot
    plt = plot(; xlabel="Time Index", ylabel="x", legend=:topright,
               size=(600, 400), palette=:tab10, PLOT_DEFAULTS...)

    for particle in 1:data.n_particles
        plot!(plt, data.trajectories[particle][:, 1],
            label="Particle $particle", linewidth=1.5)

        scatter!(plt, [1], [data.initial_positions[particle][1]],
            marker=:circle, markersize=6,
            label=(particle == 1 ? "Initial" : ""))
        scatter!(plt, [size(data.trajectories[particle], 1)], [data.final_positions[particle][1]],
            marker=:square, markersize=6,
            label=(particle == 1 ? "Final" : ""))
    end

    return plt
end

function _plot_trajectories_2d(data::TrajectoryData)::Plots.Plot
    plt = plot(; xlabel="x", ylabel="y", aspect_ratio=:equal, legend=:topright,
               size=(500, 500), palette=:tab10, PLOT_DEFAULTS...)

    for particle in 1:data.n_particles
        plot!(plt, data.trajectories[particle][:, 1], data.trajectories[particle][:, 2],
            label="Particle $particle", linewidth=1.5)

        scatter!(plt, [data.initial_positions[particle][1]], [data.initial_positions[particle][2]],
            marker=:circle, markersize=6,
            label=(particle == 1 ? "Initial" : ""))
        scatter!(plt, [data.final_positions[particle][1]], [data.final_positions[particle][2]],
            marker=:square, markersize=6,
            label=(particle == 1 ? "Final" : ""))
    end

    return plt
end

function _plot_trajectories_3d(data::TrajectoryData)::Plots.Plot
    plt = plot(; xlabel="x", ylabel="y", zlabel="z", legend=:topright,
               size=(600, 500), palette=:tab10, PLOT_DEFAULTS...)

    for particle in 1:data.n_particles
        plot!(plt, data.trajectories[particle][:, 1],
              data.trajectories[particle][:, 2],
              data.trajectories[particle][:, 3],
            label="Particle $particle", linewidth=1.5)

        scatter!(plt, [data.initial_positions[particle][1]],
                [data.initial_positions[particle][2]],
                [data.initial_positions[particle][3]],
            marker=:circle, markersize=6,
            label=(particle == 1 ? "Initial" : ""))
        scatter!(plt, [data.final_positions[particle][1]],
                [data.final_positions[particle][2]],
                [data.final_positions[particle][3]],
            marker=:square, markersize=6,
            label=(particle == 1 ? "Final" : ""))
    end

    return plt
end
