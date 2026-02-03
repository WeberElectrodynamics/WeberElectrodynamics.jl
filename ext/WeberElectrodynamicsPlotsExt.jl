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

"""
    plot_pair_forces(data::PairForceData) -> Plot

Plot comprehensive Weber force analysis for a particle pair.

Four-panel vertical layout (1×4):
1. Force magnitude |F| with min/max/mean/range statistics
2. Force components (Fx, Fy, [Fz]) as stacked colored lines
3. Vector form decomposition (4 terms: Coulomb, v·v, r·a, -(r̂·v)²)
4. Radial form decomposition (3 terms: Coulomb, -ṙ²/2c², r·r̈/c²)

Decomposition terms are shown as signed scalars (positive = repulsion, negative = attraction).
"""
function WeberElectrodynamics.plot_pair_forces(data::PairForceData)::Plots.Plot
    n_times = length(data.t)
    dims = data.dims
    i, j = data.pair

    # Color palette for force components
    component_colors = [:steelblue, :firebrick, :forestgreen]
    component_labels = ["Fx", "Fy", "Fz"]

    # Color palette for decomposition terms
    decomp_colors = [:black, :steelblue, :firebrick, :forestgreen]

    # =========================================================================
    # Panel 1: Force Magnitude with Statistics
    # =========================================================================
    p1 = plot(;
        title = "Force Magnitude |F| — Pair ($i,$j)",
        xlabel = "",
        ylabel = "|F|",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p1, data.t, data.magnitude, label = "", linewidth = 1.5, color = :black)

    # Statistics reference lines
    stats = data.stats
    hline!(
        p1,
        [stats.max],
        linestyle = :dash,
        linewidth = 1,
        color = :firebrick,
        label = "max = $(_format_scientific(stats.max))",
    )
    hline!(
        p1,
        [stats.min],
        linestyle = :dash,
        linewidth = 1,
        color = :steelblue,
        label = "min = $(_format_scientific(stats.min))",
    )
    hline!(
        p1,
        [stats.mean],
        linestyle = :dot,
        linewidth = 1,
        color = :gray,
        label = "mean = $(_format_scientific(stats.mean))",
    )

    # =========================================================================
    # Panel 2: Force Components (Fx, Fy, Fz)
    # =========================================================================
    p2 = plot(;
        title = "Force Components",
        xlabel = "",
        ylabel = "Force",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )

    for d = 1:dims
        component = [data.force[t][d] for t = 1:n_times]
        plot!(
            p2,
            data.t,
            component,
            label = component_labels[d],
            linewidth = 1.5,
            color = component_colors[d],
        )
    end
    hline!(p2, [0.0], linestyle = :dash, color = :gray, linewidth = 0.5, label = "")

    # =========================================================================
    # Panel 3: Vector Form Decomposition (4 terms as signed scalars)
    # =========================================================================
    p3 = plot(;
        title = "Vector Form: F = Coulomb·(1 + (v·v + r·a - 1.5(r̂·v)²)/c²)",
        xlabel = "",
        ylabel = "Signed Force",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )

    # Signed scalars: positive = repulsion, negative = attraction
    # sign(k) gives physics; sign(dot(term, coulomb)) gives coefficient sign
    k_sign = sign(data.charge_product)
    coulomb_signed = [k_sign * norm(data.coulomb[t]) for t = 1:n_times]
    vv_signed = [sign(dot(data.vector_term_vv[t], data.coulomb[t])) * k_sign * norm(data.vector_term_vv[t]) for t = 1:n_times]
    ra_signed = [sign(dot(data.vector_term_ra[t], data.coulomb[t])) * k_sign * norm(data.vector_term_ra[t]) for t = 1:n_times]
    rv2_signed = [sign(dot(data.vector_term_rv2[t], data.coulomb[t])) * k_sign * norm(data.vector_term_rv2[t]) for t = 1:n_times]

    plot!(p3, data.t, coulomb_signed, label = "Coulomb", linewidth = 1.5, color = decomp_colors[1])
    plot!(p3, data.t, vv_signed, label = "v·v/c²", linewidth = 1.5, color = decomp_colors[2])
    plot!(p3, data.t, ra_signed, label = "r·a/c²", linewidth = 1.5, color = decomp_colors[3])
    plot!(p3, data.t, rv2_signed, label = "-1.5(r̂·v)²/c²", linewidth = 1.5, color = decomp_colors[4])
    hline!(p3, [0.0], linestyle = :dash, color = :gray, linewidth = 0.5, label = "")

    # =========================================================================
    # Panel 4: Radial Form Decomposition (3 terms as signed scalars)
    # =========================================================================
    p4 = plot(;
        title = "Radial Form: F = Coulomb·(1 - ṙ²/(2c²) + r·r̈/c²)",
        xlabel = "Time t",
        ylabel = "Signed Force",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )

    # Reuse coulomb_signed and k_sign from above
    rdot2_signed = [sign(dot(data.radial_term_rdot2[t], data.coulomb[t])) * k_sign * norm(data.radial_term_rdot2[t]) for t = 1:n_times]
    rddot_signed = [sign(dot(data.radial_term_rddot[t], data.coulomb[t])) * k_sign * norm(data.radial_term_rddot[t]) for t = 1:n_times]

    plot!(p4, data.t, coulomb_signed, label = "Coulomb", linewidth = 1.5, color = decomp_colors[1])
    plot!(p4, data.t, rdot2_signed, label = "-ṙ²/(2c²)", linewidth = 1.5, color = decomp_colors[2])
    plot!(p4, data.t, rddot_signed, label = "r·r̈/c²", linewidth = 1.5, color = decomp_colors[3])
    hline!(p4, [0.0], linestyle = :dash, color = :gray, linewidth = 0.5, label = "")

    # =========================================================================
    # Combine into 1×4 vertical layout
    # =========================================================================
    plt = plot(
        p1, p2, p3, p4,
        layout = grid(4, 1, heights = [0.25, 0.25, 0.25, 0.25]),
        size = (SINGLE_COLUMN_WIDTH, round(Int, SINGLE_COLUMN_WIDTH * 2)),
    )

    return plt
end

"""
    plot_phase_space(data::PairForceData) -> Plot

Plot (r, ṙ) phase space portrait for a particle pair.

Shows the trajectory in separation distance vs radial velocity space,
with markers for initial and final states.
"""
function WeberElectrodynamics.plot_phase_space(data::PairForceData)::Plots.Plot
    ps = data.phase_space
    i, j = data.pair

    plt = plot(;
        title = "Phase Space Portrait — Pair ($i,$j)",
        xlabel = "Separation r",
        ylabel = "Radial Velocity ṙ",
        size = _square_size(),
        legend = :outertopright,
        legendfontsize = 10,
        PLOT_DEFAULTS...,
    )

    plot!(
        plt,
        ps.separation_distance,
        ps.radial_velocity,
        label = "",
        linewidth = 1.5,
        color = :black,
    )

    scatter!(
        plt,
        [ps.separation_distance[1]],
        [ps.radial_velocity[1]],
        marker = :circle,
        markersize = 6,
        color = :steelblue,
        label = "Initial",
    )
    scatter!(
        plt,
        [ps.separation_distance[end]],
        [ps.radial_velocity[end]],
        marker = :square,
        markersize = 6,
        color = :firebrick,
        label = "Final",
    )

    return plt
end

"""
    plot_momentum(data::MomentumData) -> Plot

Plot total momentum timeseries with comprehensive visualization.

Two-panel layout:
1. Linear momentum components (Px, Py, [Pz]) with magnitude |P|
2. Angular momentum (2D: scalar Lz, 3D: components and magnitude)
   For 1D systems, shows only the linear momentum panel.

Conservation of momentum is indicated by horizontal lines staying constant.
"""
function WeberElectrodynamics.plot_momentum(data::MomentumData)::Plots.Plot
    dims = data.dims

    # Color palette for momentum components
    component_colors = [:steelblue, :firebrick, :forestgreen]
    component_labels = ["Px", "Py", "Pz"]

    # =========================================================================
    # Panel 1: Linear Momentum Components + Magnitude
    # =========================================================================
    p1 = plot(;
        title = "Total Linear Momentum",
        xlabel = dims == 1 || isnothing(data.angular_momentum) ? "Time t" : "",
        ylabel = "Momentum P",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )

    # Plot each component
    for d = 1:dims
        plot!(
            p1,
            data.t,
            data.linear_momentum_components[:, d],
            label = component_labels[d],
            linewidth = 1.5,
            color = component_colors[d],
        )
    end

    # Plot magnitude
    plot!(
        p1,
        data.t,
        data.linear_momentum_magnitude,
        label = "|P|",
        linewidth = 2,
        color = :black,
        linestyle = :dash,
    )

    # Zero reference line
    hline!(p1, [0.0], linestyle = :dot, color = :gray, linewidth = 0.5, label = "")

    # For 1D systems, return single panel
    if dims == 1 || isnothing(data.angular_momentum)
        return plot(p1, size = _single_panel_size())
    end

    # =========================================================================
    # Panel 2: Angular Momentum
    # =========================================================================
    if dims == 2
        # 2D: scalar angular momentum Lz
        p2 = plot(;
            title = "Total Angular Momentum (z-component)",
            xlabel = "Time t",
            ylabel = "Angular Momentum Lz",
            legend = :outertopright,
            PLOT_DEFAULTS...,
        )

        plot!(
            p2,
            data.t,
            data.angular_momentum,
            label = "Lz",
            linewidth = 1.5,
            color = :black,
        )

        # Zero reference line
        hline!(p2, [0.0], linestyle = :dot, color = :gray, linewidth = 0.5, label = "")

    else  # dims == 3
        # 3D: angular momentum vector components + magnitude
        L_labels = ["Lx", "Ly", "Lz"]

        p2 = plot(;
            title = "Total Angular Momentum",
            xlabel = "Time t",
            ylabel = "Angular Momentum L",
            legend = :outertopright,
            PLOT_DEFAULTS...,
        )

        for d = 1:3
            L_component = [data.angular_momentum[t][d] for t = 1:length(data.t)]
            plot!(
                p2,
                data.t,
                L_component,
                label = L_labels[d],
                linewidth = 1.5,
                color = component_colors[d],
            )
        end

        # Plot magnitude
        plot!(
            p2,
            data.t,
            data.angular_momentum_magnitude,
            label = "|L|",
            linewidth = 2,
            color = :black,
            linestyle = :dash,
        )

        # Zero reference line
        hline!(p2, [0.0], linestyle = :dot, color = :gray, linewidth = 0.5, label = "")
    end

    # =========================================================================
    # Combine into 2-panel layout
    # =========================================================================
    plt = plot(
        p1, p2,
        layout = grid(2, 1, heights = [0.5, 0.5]),
        size = _multi_panel_size(2),
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
