module WeberElectrodynamicsPlotsExt

using WeberElectrodynamics
using Plots
using LinearAlgebra
using Printf

# Publication figure sizes (pixels at 300 DPI)
const SINGLE_COLUMN_WIDTH = 1050   # 3.5 inches at 300 DPI
const DOUBLE_COLUMN_WIDTH = 2100   # 7.0 inches at 300 DPI

const PLOT_DEFAULTS = (
    framestyle = :box,
    grid = true,
    gridalpha = 0.2,
    tickfontsize = 10,
    guidefontsize = 11,
    legendfontsize = 9,
    titlefontsize = 12,
    dpi = 300,
    linewidth = 1.5,
    markersize = 4,
    margin = 8Plots.mm,
    left_margin = 15Plots.mm,
    bottom_margin = 12Plots.mm,
    top_margin = 8Plots.mm,
    right_margin = 15Plots.mm,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
)

# =============================================================================
# Helper Functions
# =============================================================================

"""Format a number in scientific notation for plot labels."""
function _format_scientific(x::Float64; sigdigits::Int = 2)::String
    if x == 0.0
        return "0"
    elseif isnan(x) || isinf(x)
        return string(x)
    end
    return @sprintf("%.2e", x)
end

"""Compute local errors |E_t - E_{t-1}| from energy timeseries."""
function _compute_local_errors(energy::Vector{Float64})::Vector{Float64}
    n = length(energy)
    errors = Vector{Float64}(undef, n - 1)
    @inbounds for i = 2:n
        errors[i-1] = abs(energy[i] - energy[i-1])
    end
    return errors
end

"""Compute global error percentage |(E_t - E_0)/E_0| * 100."""
function _compute_global_error_percent(energy::Vector{Float64})::Vector{Float64}
    E0 = energy[1]
    if abs(E0) < 100 * eps(Float64)
        return fill(NaN, length(energy))
    end
    return abs.((energy .- E0) ./ E0) .* 100.0
end

"""Standard figure size for single-panel plots (4:3 aspect ratio)."""
_single_panel_size() = (SINGLE_COLUMN_WIDTH, round(Int, SINGLE_COLUMN_WIDTH * 0.75))

"""Standard figure size for multi-panel vertical layouts."""
_multi_panel_size(n::Int) = (SINGLE_COLUMN_WIDTH, round(Int, SINGLE_COLUMN_WIDTH * 0.5 * n))

"""Standard figure size for square plots (phase space, 2D trajectories)."""
_square_size() = (SINGLE_COLUMN_WIDTH, SINGLE_COLUMN_WIDTH)

# =============================================================================
# Public API
# =============================================================================

"""Plot particle trajectories (1D, 2D, or 3D based on `data.dims`)."""
function WeberElectrodynamics.plot_trajectories(data::TrajectoryData)::Plots.Plot
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

"""
    plot_energy(data::EnergyData) -> Plot

Plot total n-body energy timeseries with comprehensive error visualization.

Two-panel layout:
1. Energy components (Total, Kinetic, Potential)
2. Relative error on log scale with max/avg reference lines
"""
function WeberElectrodynamics.plot_energy(data::EnergyData)::Plots.Plot
    E0 = data.total_energy[1]
    relative_error = abs.((data.total_energy .- E0) ./ E0)
    relative_error = max.(relative_error, eps(Float64))

    # Panel 1: Energy components
    p1 = plot(;
        title = "Energy Conservation",
        xlabel = "",
        ylabel = "Energy",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p1, data.t, data.total_energy, label = "Total E", linewidth = 2, color = :black)
    plot!(
        p1,
        data.t,
        data.kinetic_energy,
        label = "Kinetic T",
        linewidth = 1.5,
        color = :steelblue,
    )
    plot!(
        p1,
        data.t,
        data.total_potential_energy,
        label = "Potential U",
        linewidth = 1.5,
        color = :firebrick,
    )

    # Panel 2: Relative error with enhanced annotations
    max_err = maximum(relative_error)
    avg_err = length(relative_error) > 1 ? sum(relative_error[2:end]) / (length(relative_error) - 1) : max_err
    max_idx = argmax(relative_error)

    p2 = plot(;
        title = "Relative Energy Error",
        xlabel = "Time t",
        ylabel = "|ΔE/E₀|",
        legend = :outertopright,
        yscale = :log10,
        PLOT_DEFAULTS...,
    )
    plot!(p2, data.t, relative_error, label = "", linewidth = 1.5, color = :black)
    hline!(
        p2,
        [max_err],
        linestyle = :dash,
        linewidth = 1,
        color = :firebrick,
        label = "max = $(_format_scientific(max_err))",
    )
    hline!(
        p2,
        [avg_err],
        linestyle = :dot,
        linewidth = 1,
        color = :gray,
        label = "avg = $(_format_scientific(avg_err))",
    )
    # Vertical marker at max error time
    vline!(p2, [data.t[max_idx]], linestyle = :dash, color = :gray, alpha = 0.5, label = "")

    plt = plot(p1, p2, layout = grid(2, 1, heights = [0.6, 0.4]), size = _multi_panel_size(2))
    return plt
end

"""
    plot_pair_energy(data::EnergyData, pair::Tuple{Int,Int}) -> Plot

Plot energy components and radial dynamics for a specific particle pair.

Two-panel layout:
1. Coulomb term, velocity term, and total pair potential
2. Radial velocity with approach/separation shading

# Arguments
- `data::EnergyData`: Energy data from `compute_energy_timeseries`
- `pair::Tuple{Int,Int}`: Particle pair indices, e.g., `(1, 2)`
"""
function WeberElectrodynamics.plot_pair_energy(
    data::EnergyData,
    pair::Tuple{Int,Int},
)::Plots.Plot
    # Validation
    if !haskey(data.pair_energies, pair)
        available = sort(collect(keys(data.pair_energies)))
        throw(ArgumentError("Pair $pair not found. Available pairs: $available"))
    end

    pair_data = data.pair_energies[pair]
    i, j = pair

    # Panel 1: Energy components
    p1 = plot(;
        title = "Pair ($i,$j) Weber Potential",
        xlabel = "",
        ylabel = "Energy",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(
        p1,
        data.t,
        pair_data.coulomb_term,
        label = "Coulomb qᵢqⱼ/r",
        linewidth = 1.5,
        color = :firebrick,
    )
    plot!(
        p1,
        data.t,
        pair_data.velocity_term,
        label = "Velocity term",
        linewidth = 1.5,
        color = :steelblue,
    )
    plot!(
        p1,
        data.t,
        pair_data.total_pair_potential,
        label = "Total",
        linewidth = 2,
        color = :black,
    )

    # Panel 2: Radial velocity with sign-based shading
    rdot = pair_data.radial_velocity

    p2 = plot(;
        title = "Radial Dynamics",
        xlabel = "Time t",
        ylabel = "Radial Velocity ṙ",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )

    # Fill between curve and zero based on sign
    # Positive rdot (separating) - firebrick
    # Negative rdot (approaching) - steelblue
    rdot_pos = [v >= 0 ? v : 0.0 for v in rdot]
    rdot_neg = [v < 0 ? v : 0.0 for v in rdot]

    plot!(
        p2,
        data.t,
        rdot_pos,
        fillrange = 0,
        fillalpha = 0.3,
        fillcolor = :firebrick,
        linewidth = 0,
        label = "",
    )
    plot!(
        p2,
        data.t,
        rdot_neg,
        fillrange = 0,
        fillalpha = 0.3,
        fillcolor = :steelblue,
        linewidth = 0,
        label = "",
    )

    # Overlay the actual line
    plot!(p2, data.t, rdot, label = "ṙ", linewidth = 1.5, color = :black)

    # Zero reference line
    hline!(p2, [0.0], linestyle = :dash, color = :gray, linewidth = 1, label = "")

    # Legend entries for shading (invisible data points)
    plot!(
        p2,
        Float64[],
        Float64[],
        fillrange = 0,
        fillalpha = 0.3,
        fillcolor = :steelblue,
        label = "Approaching (ṙ<0)",
        linewidth = 0,
    )
    plot!(
        p2,
        Float64[],
        Float64[],
        fillrange = 0,
        fillalpha = 0.3,
        fillcolor = :firebrick,
        label = "Separating (ṙ>0)",
        linewidth = 0,
    )

    plt = plot(p1, p2, layout = grid(2, 1, heights = [0.6, 0.4]), size = _multi_panel_size(2))
    return plt
end

"""
    plot_energy_errors(data::EnergyData) -> Plot

Comprehensive error analysis visualization.

Three-panel layout:
1. Local error timeseries |E_t - E_{t-1}| (log scale)
2. Global error percentage over time (log scale)
3. Hamiltonian validation error (log scale)
"""
function WeberElectrodynamics.plot_energy_errors(data::EnergyData)::Plots.Plot
    # Compute derived quantities
    local_errors = _compute_local_errors(data.total_energy)
    local_errors_plot = max.(local_errors, eps(Float64))
    global_percent = _compute_global_error_percent(data.total_energy)
    global_percent_plot = max.(global_percent, eps(Float64))
    h_error = max.(data.hamiltonian_validation_error, eps(Float64))

    stats = data.statistics

    # Panel 1: Local error
    p1 = plot(;
        title = "Local Energy Error",
        xlabel = "",
        ylabel = "|Eₜ - Eₜ₋₁|",
        yscale = :log10,
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p1, data.t[2:end], local_errors_plot, label = "", linewidth = 1, color = :steelblue)
    hline!(
        p1,
        [stats.local_error_max],
        linestyle = :dash,
        color = :firebrick,
        linewidth = 1,
        label = "max = $(_format_scientific(stats.local_error_max))",
    )
    hline!(
        p1,
        [stats.local_error_avg],
        linestyle = :dot,
        color = :gray,
        linewidth = 1,
        label = "avg = $(_format_scientific(stats.local_error_avg))",
    )

    # Panel 2: Global error percentage
    p2 = plot(;
        title = "Global Energy Drift",
        xlabel = "",
        ylabel = "Global Error (%)",
        yscale = :log10,
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p2, data.t, global_percent_plot, label = "", linewidth = 1.5, color = :firebrick)
    hline!(
        p2,
        [stats.global_error_percent_max],
        linestyle = :dash,
        color = :black,
        linewidth = 1,
        label = "max = $(_format_scientific(stats.global_error_percent_max))%",
    )

    # Panel 3: Hamiltonian validation
    max_h_err = maximum(data.hamiltonian_validation_error)

    p3 = plot(;
        title = "Hamiltonian Validation",
        xlabel = "Time t",
        ylabel = "|H_computed - H_compiled|",
        yscale = :log10,
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p3, data.t, h_error, label = "", linewidth = 1, color = :black)
    hline!(
        p3,
        [max.(max_h_err, eps(Float64))],
        linestyle = :dash,
        color = :gray,
        linewidth = 1,
        label = "max = $(_format_scientific(max_h_err))",
    )

    plt = plot(p1, p2, p3, layout = grid(3, 1, heights = [0.35, 0.35, 0.30]), size = _multi_panel_size(3))
    return plt
end

"""Plot force magnitudes between particle pairs over time."""
function WeberElectrodynamics.plot_forces(data::ForceData)::Plots.Plot
    n_times = length(data.t)
    n = data.n_particles

    plt = plot(;
        title = "Weber Force Magnitudes",
        xlabel = "Time t",
        ylabel = "|F|",
        legend = :outertopright,
        size = _single_panel_size(),
        palette = :tab10,
        PLOT_DEFAULTS...,
    )

    for i = 1:n
        for j = (i+1):n
            F_ij = data.forces[(i, j)]

            magnitudes = Vector{Float64}(undef, n_times)
            for t = 1:n_times
                magnitudes[t] = norm(F_ij[t])
            end

            plot!(plt, data.t, magnitudes, label = "F($i,$j)", linewidth = 1.5)
        end
    end

    return plt
end

"""Plot (r, ṙ) phase space portrait."""
function WeberElectrodynamics.plot_phase_space(data::PhaseSpaceData)::Plots.Plot
    plt = plot(;
        title = "Phase Space Portrait",
        xlabel = "Separation r",
        ylabel = "Radial Velocity ṙ",
        size = _square_size(),
        legend = :outertopright,
        legendfontsize = 10,
        PLOT_DEFAULTS...,
    )

    plot!(
        plt,
        data.separation_distance,
        data.radial_velocity,
        label = "",
        linewidth = 1.5,
        color = :black,
    )

    scatter!(
        plt,
        [data.separation_distance[1]],
        [data.radial_velocity[1]],
        marker = :circle,
        markersize = 6,
        color = :steelblue,
        label = "Initial",
    )
    scatter!(
        plt,
        [data.separation_distance[end]],
        [data.radial_velocity[end]],
        marker = :square,
        markersize = 6,
        color = :firebrick,
        label = "Final",
    )

    return plt
end

# =============================================================================
# Internal helpers
# =============================================================================

function _plot_trajectories_1d(data::TrajectoryData)::Plots.Plot
    plt = plot(;
        title = "1D Particle Trajectories",
        xlabel = "Time Index",
        ylabel = "x",
        legend = :outertopright,
        size = _single_panel_size(),
        palette = :tab10,
        PLOT_DEFAULTS...,
    )

    for particle = 1:data.n_particles
        plot!(
            plt,
            data.trajectories[particle][:, 1],
            label = "Particle $particle",
            linewidth = 1.5,
        )

        scatter!(
            plt,
            [1],
            [data.initial_positions[particle][1]],
            marker = :circle,
            markersize = 6,
            label = (particle == 1 ? "Initial" : ""),
        )
        scatter!(
            plt,
            [size(data.trajectories[particle], 1)],
            [data.final_positions[particle][1]],
            marker = :square,
            markersize = 6,
            label = (particle == 1 ? "Final" : ""),
        )
    end

    return plt
end

function _plot_trajectories_2d(data::TrajectoryData)::Plots.Plot
    plt = plot(;
        title = "2D Particle Trajectories",
        xlabel = "x",
        ylabel = "y",
        aspect_ratio = :equal,
        legend = :outertopright,
        size = _square_size(),
        palette = :tab10,
        PLOT_DEFAULTS...,
    )

    for particle = 1:data.n_particles
        plot!(
            plt,
            data.trajectories[particle][:, 1],
            data.trajectories[particle][:, 2],
            label = "Particle $particle",
            linewidth = 1.5,
        )

        scatter!(
            plt,
            [data.initial_positions[particle][1]],
            [data.initial_positions[particle][2]],
            marker = :circle,
            markersize = 6,
            label = (particle == 1 ? "Initial" : ""),
        )
        scatter!(
            plt,
            [data.final_positions[particle][1]],
            [data.final_positions[particle][2]],
            marker = :square,
            markersize = 6,
            label = (particle == 1 ? "Final" : ""),
        )
    end

    return plt
end

function _plot_trajectories_3d(data::TrajectoryData)::Plots.Plot
    plt = plot(;
        title = "3D Particle Trajectories",
        xlabel = "x",
        ylabel = "y",
        zlabel = "z",
        legend = :outertopright,
        size = (SINGLE_COLUMN_WIDTH, round(Int, SINGLE_COLUMN_WIDTH * 0.85)),
        palette = :tab10,
        PLOT_DEFAULTS...,
    )

    for particle = 1:data.n_particles
        plot!(
            plt,
            data.trajectories[particle][:, 1],
            data.trajectories[particle][:, 2],
            data.trajectories[particle][:, 3],
            label = "Particle $particle",
            linewidth = 1.5,
        )

        scatter!(
            plt,
            [data.initial_positions[particle][1]],
            [data.initial_positions[particle][2]],
            [data.initial_positions[particle][3]],
            marker = :circle,
            markersize = 6,
            label = (particle == 1 ? "Initial" : ""),
        )
        scatter!(
            plt,
            [data.final_positions[particle][1]],
            [data.final_positions[particle][2]],
            [data.final_positions[particle][3]],
            marker = :square,
            markersize = 6,
            label = (particle == 1 ? "Final" : ""),
        )
    end

    return plt
end

end # module
