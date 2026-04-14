module WeberElectrodynamicsPlotsExt

using WeberElectrodynamics
using WeberElectrodynamics: @sprintf, norm, dot
using Plots
using LaTeXStrings

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

"""
    plot_trajectories(data::TrajectoryData) -> Plot

Plot particle trajectories in 1D, 2D, or 3D based on `data.dims`.

- **1D**: Position vs. time for each particle.
- **2D**: x–y scatter with equal aspect ratio and start/end markers.
- **3D**: x–y–z line plot with start/end markers.

# Arguments
- `data::TrajectoryData`: Trajectory data from `compute_trajectory_data`.

# Returns
- A `Plots.Plot` object.
"""
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

Plot total n-body energy timeseries showing Total, Kinetic, and Potential components.
"""
function WeberElectrodynamics.plot_energy(data::EnergyData)::Plots.Plot
    p1 = plot(;
        title = "Energy Conservation",
        xlabel = L"t",
        ylabel = L"E\ (\mathrm{energy})",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p1, data.t, data.total_energy, label = L"H = T + U", linewidth = 2, color = :black)
    plot!(
        p1,
        data.t,
        data.kinetic_energy,
        label = L"T\ (\mathrm{kinetic})",
        linewidth = 1.5,
        color = :steelblue,
    )
    plot!(
        p1,
        data.t,
        data.total_potential_energy,
        label = L"U\ (\mathrm{potential})",
        linewidth = 1.5,
        color = :firebrick,
    )

    return plot(p1; size = _multi_panel_size(1))
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
        ylabel = L"U_{ij}\ (\mathrm{pair\ potential})",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(
        p1,
        data.t,
        pair_data.coulomb_term,
        label = L"q_i q_j / r\ (\mathrm{Coulomb})",
        linewidth = 1.5,
        color = :firebrick,
    )
    plot!(
        p1,
        data.t,
        pair_data.velocity_term,
        label = L"-q_i q_j\, \dot r^{\,2}/(2 c^2 r)\ (\mathrm{velocity})",
        linewidth = 1.5,
        color = :steelblue,
    )
    plot!(
        p1,
        data.t,
        pair_data.total_pair_potential,
        label = L"U_{ij}\ (\mathrm{total})",
        linewidth = 2,
        color = :black,
    )

    # Panel 2: Radial velocity with sign-based shading
    rdot = pair_data.radial_velocity

    p2 = plot(;
        title = "Radial Dynamics",
        xlabel = L"t",
        ylabel = L"\dot r\ (\mathrm{radial\ velocity})",
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
    plot!(p2, data.t, rdot, label = L"\dot r", linewidth = 1.5, color = :black)

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
        label = L"\dot r < 0\ (\mathrm{approaching})",
        linewidth = 0,
    )
    plot!(
        p2,
        Float64[],
        Float64[],
        fillrange = 0,
        fillalpha = 0.3,
        fillcolor = :firebrick,
        label = L"\dot r > 0\ (\mathrm{separating})",
        linewidth = 0,
    )

    plt = plot(p1, p2, layout = grid(2, 1, heights = [0.6, 0.4]), size = _multi_panel_size(2))
    return plt
end

"""
    plot_energy_errors(data::EnergyData) -> Plot

Energy error diagnostics.

Two-panel layout:
1. Local error timeseries |E_t - E_{t-1}| (log scale)
2. Relative energy error |ΔE/E₀| (log scale)
"""
function WeberElectrodynamics.plot_energy_errors(data::EnergyData)::Plots.Plot
    local_errors = _compute_local_errors(data.total_energy)
    local_errors_plot = max.(local_errors, eps(Float64))

    E0 = data.total_energy[1]
    relative_error = abs.((data.total_energy .- E0) ./ E0)
    relative_error_plot = max.(relative_error, eps(Float64))

    stats = data.statistics

    # Panel 1: Local error
    p1 = plot(;
        title = "Local Energy Error",
        xlabel = "",
        ylabel = L"\left\lvert E_t - E_{t-1} \right\rvert",
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
        label = latexstring("\\max = ", _format_scientific(stats.local_error_max)),
    )
    hline!(
        p1,
        [stats.local_error_avg],
        linestyle = :dot,
        color = :gray,
        linewidth = 1,
        label = latexstring("\\mathrm{avg} = ", _format_scientific(stats.local_error_avg)),
    )

    # Panel 2: Relative energy error
    max_rel_err = maximum(relative_error)
    avg_rel_err = length(relative_error) > 1 ?
        sum(relative_error[2:end]) / (length(relative_error) - 1) : max_rel_err

    p2 = plot(;
        title = "Relative Energy Error",
        xlabel = L"t",
        ylabel = L"\left\lvert \Delta E / E_0 \right\rvert",
        yscale = :log10,
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(p2, data.t, relative_error_plot, label = "", linewidth = 1.5, color = :firebrick)
    hline!(
        p2,
        [max(max_rel_err, eps(Float64))],
        linestyle = :dash,
        color = :black,
        linewidth = 1,
        label = latexstring("\\max = ", _format_scientific(max_rel_err)),
    )
    hline!(
        p2,
        [max(avg_rel_err, eps(Float64))],
        linestyle = :dot,
        color = :gray,
        linewidth = 1,
        label = latexstring("\\mathrm{avg} = ", _format_scientific(avg_rel_err)),
    )

    plt = plot(p1, p2, layout = grid(2, 1, heights = [0.5, 0.5]), size = _multi_panel_size(2))
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
    component_labels = [L"F_x", L"F_y", L"F_z"]

    # Color palette for decomposition terms
    decomp_colors = [:black, :steelblue, :firebrick, :forestgreen]

    # =========================================================================
    # Panel 1: Force Magnitude with Statistics
    # =========================================================================
    p1 = plot(;
        title = "Force Magnitude — Pair ($i,$j)",
        xlabel = "",
        ylabel = L"\lVert \mathbf F_{ij} \rVert",
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
        label = latexstring("\\max = ", _format_scientific(stats.max)),
    )
    hline!(
        p1,
        [stats.min],
        linestyle = :dash,
        linewidth = 1,
        color = :steelblue,
        label = latexstring("\\min = ", _format_scientific(stats.min)),
    )
    hline!(
        p1,
        [stats.mean],
        linestyle = :dot,
        linewidth = 1,
        color = :gray,
        label = latexstring("\\mathrm{mean} = ", _format_scientific(stats.mean)),
    )

    # =========================================================================
    # Panel 2: Force Components (Fx, Fy, Fz)
    # =========================================================================
    p2 = plot(;
        title = "Force Components",
        xlabel = "",
        ylabel = L"F_k\ (\mathrm{component})",
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
        ylabel = L"F\ (\mathrm{signed})",
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

    plot!(p3, data.t, coulomb_signed, label = L"F_{\mathrm{C}} = q_i q_j/r^2", linewidth = 1.5, color = decomp_colors[1])
    plot!(p3, data.t, vv_signed, label = L"\mathbf v\cdot\mathbf v/c^2", linewidth = 1.5, color = decomp_colors[2])
    plot!(p3, data.t, ra_signed, label = L"\mathbf r\cdot\mathbf a/c^2", linewidth = 1.5, color = decomp_colors[3])
    plot!(p3, data.t, rv2_signed, label = L"-\tfrac{3}{2}(\hat{\mathbf r}\cdot\mathbf v)^2/c^2", linewidth = 1.5, color = decomp_colors[4])
    hline!(p3, [0.0], linestyle = :dash, color = :gray, linewidth = 0.5, label = "")

    # =========================================================================
    # Panel 4: Radial Form Decomposition (3 terms as signed scalars)
    # =========================================================================
    p4 = plot(;
        title = "Radial Form: F = Coulomb·(1 - ṙ²/(2c²) + r·r̈/c²)",
        xlabel = L"t",
        ylabel = L"F\ (\mathrm{signed})",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )

    # Reuse coulomb_signed and k_sign from above
    rdot2_signed = [sign(dot(data.radial_term_rdot2[t], data.coulomb[t])) * k_sign * norm(data.radial_term_rdot2[t]) for t = 1:n_times]
    rddot_signed = [sign(dot(data.radial_term_rddot[t], data.coulomb[t])) * k_sign * norm(data.radial_term_rddot[t]) for t = 1:n_times]

    plot!(p4, data.t, coulomb_signed, label = L"F_{\mathrm{C}} = q_i q_j/r^2", linewidth = 1.5, color = decomp_colors[1])
    plot!(p4, data.t, rdot2_signed, label = L"-\dot r^{\,2}/(2c^2)", linewidth = 1.5, color = decomp_colors[2])
    plot!(p4, data.t, rddot_signed, label = L"r\,\ddot r/c^2", linewidth = 1.5, color = decomp_colors[3])
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
        xlabel = L"r\ (\mathrm{separation})",
        ylabel = L"\dot r\ (\mathrm{radial\ velocity})",
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
        label = L"(r_0,\, \dot r_0)",
    )
    scatter!(
        plt,
        [ps.separation_distance[end]],
        [ps.radial_velocity[end]],
        marker = :square,
        markersize = 6,
        color = :firebrick,
        label = L"(r_f,\, \dot r_f)",
    )

    return plt
end

"""
    plot_momentum_errors(data::MomentumData) -> Plot

Plot conservation errors for linear and angular momentum as two stacked
log-scale panels.

- **Top panel** — linear drift `‖P(t) − P(0)‖` (steelblue).
- **Bottom panel** — angular drift `|L(t) − L(0)|` in 2D or `‖L − L₀‖` in
  3D (firebrick).

Each panel's legend reports the max absolute drift over the run, and —
when the initial magnitude is nonzero — the max relative drift
`max_t ‖Δ·‖ / ‖·₀‖`. Relative error is a scalar rescaling of the absolute
curve (the initial magnitude is constant), so on a log plot it is simply
a shifted copy; reporting it in the label avoids a redundant curve.

1D systems show only the linear panel.
"""
function WeberElectrodynamics.plot_momentum_errors(data::MomentumData)::Plots.Plot
    dims = data.dims
    nt = length(data.t)

    # --- Linear momentum drift: ‖P(t) − P(0)‖ ---
    P0 = data.linear_momentum_components[1, :]
    P0_mag = data.linear_momentum_magnitude[1]
    dP = Vector{Float64}(undef, nt)
    @inbounds for i = 1:nt
        s = 0.0
        for d = 1:dims
            s += (data.linear_momentum_components[i, d] - P0[d])^2
        end
        dP[i] = sqrt(s)
    end
    dP_max = maximum(dP)

    has_L = !(dims == 1 || isnothing(data.angular_momentum))

    linear_label = if P0_mag > 0
        latexstring(
            "\\mathrm{abs\\ max}\\ ", _format_scientific(dP_max),
            ",\\ \\mathrm{rel\\ max}\\ ", _format_scientific(dP_max / P0_mag),
        )
    else
        latexstring("\\mathrm{abs\\ max}\\ ", _format_scientific(dP_max))
    end

    p_lin = plot(;
        title = "Linear Momentum Drift ‖ΔP‖",
        xlabel = has_L ? "" : L"t",
        ylabel = L"\lVert \mathbf P(t) - \mathbf P(0) \rVert",
        yscale = :log10,
        legend = :topleft,
        PLOT_DEFAULTS...,
    )
    plot!(
        p_lin,
        data.t[2:end],
        max.(dP[2:end], eps(Float64)),
        label = linear_label,
        linewidth = 1.5,
        color = :steelblue,
    )

    if !has_L
        plot!(p_lin, size = _single_panel_size())
        return p_lin
    end

    # --- Angular momentum drift ---
    L0 = data.angular_momentum[1]
    if dims == 2
        dL = abs.(data.angular_momentum .- L0)
        L0_mag = abs(L0)
        L_sym = "|ΔLz|"
        ylabel_L = L"\lvert L_z(t) - L_z(0) \rvert"
    else  # dims == 3
        dL = Vector{Float64}(undef, nt)
        @inbounds for i = 1:nt
            s = 0.0
            for d = 1:3
                s += (data.angular_momentum[i][d] - L0[d])^2
            end
            dL[i] = sqrt(s)
        end
        L0_mag = sqrt(L0[1]^2 + L0[2]^2 + L0[3]^2)
        L_sym = "‖ΔL‖"
        ylabel_L = L"\lVert \mathbf L(t) - \mathbf L(0) \rVert"
    end
    dL_max = maximum(dL)

    angular_label = if L0_mag > 0
        latexstring(
            "\\mathrm{abs\\ max}\\ ", _format_scientific(dL_max),
            ",\\ \\mathrm{rel\\ max}\\ ", _format_scientific(dL_max / L0_mag),
        )
    else
        latexstring("\\mathrm{abs\\ max}\\ ", _format_scientific(dL_max))
    end

    p_ang = plot(;
        title = "Angular Momentum Drift $(L_sym)",
        xlabel = L"t",
        ylabel = ylabel_L,
        yscale = :log10,
        legend = :topleft,
        PLOT_DEFAULTS...,
    )
    plot!(
        p_ang,
        data.t[2:end],
        max.(dL[2:end], eps(Float64)),
        label = angular_label,
        linewidth = 1.5,
        color = :firebrick,
    )

    plt = plot(
        p_lin,
        p_ang,
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
        xlabel = L"n\ (\mathrm{time\ index})",
        ylabel = L"q_x",
        legend = :outertopright,
        size = _single_panel_size(),
        palette = :tab10,
        PLOT_DEFAULTS...,
    )

    for particle = 1:data.n_particles
        plot!(
            plt,
            data.trajectories[particle][:, 1],
            label = latexstring("\\mathrm{particle}\\ ", particle),
            linewidth = 1.5,
        )

        scatter!(
            plt,
            [1],
            [data.initial_positions[particle][1]],
            marker = :circle,
            markersize = 6,
            label = (particle == 1 ? L"\mathbf q_0" : ""),
        )
        scatter!(
            plt,
            [size(data.trajectories[particle], 1)],
            [data.final_positions[particle][1]],
            marker = :square,
            markersize = 6,
            label = (particle == 1 ? L"\mathbf q_f" : ""),
        )
    end

    return plt
end

function _plot_trajectories_2d(data::TrajectoryData)::Plots.Plot
    # Compute tight axis limits from all trajectory data
    all_x = vcat([data.trajectories[p][:, 1] for p in 1:data.n_particles]...)
    all_y = vcat([data.trajectories[p][:, 2] for p in 1:data.n_particles]...)
    xmin, xmax = extrema(all_x)
    ymin, ymax = extrema(all_y)
    dx = xmax - xmin
    dy = ymax - ymin
    margin = 0.1 * max(dx, dy, eps())
    xlims = (xmin - margin, xmax + margin)
    ylims = (ymin - margin, ymax + margin)

    plt = plot(;
        title = "2D Particle Trajectories",
        xlabel = L"q_x",
        ylabel = L"q_y",
        aspect_ratio = :equal,
        xlims = xlims,
        ylims = ylims,
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
            label = latexstring("\\mathrm{particle}\\ ", particle),
            linewidth = 1.5,
        )

        scatter!(
            plt,
            [data.initial_positions[particle][1]],
            [data.initial_positions[particle][2]],
            marker = :circle,
            markersize = 6,
            label = (particle == 1 ? L"\mathbf q_0" : ""),
        )
        scatter!(
            plt,
            [data.final_positions[particle][1]],
            [data.final_positions[particle][2]],
            marker = :square,
            markersize = 6,
            label = (particle == 1 ? L"\mathbf q_f" : ""),
        )
    end

    return plt
end

function _plot_trajectories_3d(data::TrajectoryData)::Plots.Plot
    plt = plot(;
        title = "3D Particle Trajectories",
        xlabel = L"q_x",
        ylabel = L"q_y",
        zlabel = L"q_z",
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
            label = latexstring("\\mathrm{particle}\\ ", particle),
            linewidth = 1.5,
        )

        scatter!(
            plt,
            [data.initial_positions[particle][1]],
            [data.initial_positions[particle][2]],
            [data.initial_positions[particle][3]],
            marker = :circle,
            markersize = 6,
            label = (particle == 1 ? L"\mathbf q_0" : ""),
        )
        scatter!(
            plt,
            [data.final_positions[particle][1]],
            [data.final_positions[particle][2]],
            [data.final_positions[particle][3]],
            marker = :square,
            markersize = 6,
            label = (particle == 1 ? L"\mathbf q_f" : ""),
        )
    end

    return plt
end


# =============================================================================
# Zöllner Electrogravitation Plots
# =============================================================================

"""
    plot_zollner_energy(data::EnergyData) -> Plot

Two-panel visualization of the Zöllner mismatch contribution to energy.

Panel 1: Total potential energy with the Zöllner extra contribution per pair overlaid.
Panel 2: Aggregate Zöllner residual (emergent gravitational energy) vs total energy.

When Zöllner is disabled (all κ = 1), both Zöllner traces are identically zero.
"""
function WeberElectrodynamics.plot_zollner_energy(data::EnergyData)::Plots.Plot
    colors = [:steelblue, :firebrick, :forestgreen, :darkorchid, :darkorange]

    # Panel 1: Total potential + Zöllner extra per pair
    p1 = plot(;
        title = "Potential Energy with Zöllner Contribution",
        xlabel = "",
        ylabel = L"U\ (\mathrm{energy})",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(
        p1, data.t, data.total_potential_energy,
        label = L"U_{\mathrm{total}}", linewidth = 2, color = :black,
    )
    for (pair_idx, ((i, j), pdata)) in enumerate(sort(collect(data.pair_energies), by = x -> x[1]))
        c = colors[mod1(pair_idx, length(colors))]
        kappa_str = @sprintf("%.4g", pdata.kappa)
        plot!(
            p1, data.t, pdata.total_pair_potential,
            label = latexstring("U_{", i, j, "}\\ \\kappa=", kappa_str),
            linewidth = 1.2, color = c, linestyle = :solid,
        )
        plot!(
            p1, data.t, pdata.zollner_extra_potential,
            label = latexstring("\\Delta U^{\\mathrm Z}_{", i, j, "}"),
            linewidth = 1.0, color = c, linestyle = :dash,
        )
    end

    # Panel 2: Zöllner residual vs total energy
    p2 = plot(;
        title = "Zöllner Gravitational Residual",
        xlabel = L"t",
        ylabel = L"E\ (\mathrm{energy})",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(
        p2, data.t, data.total_energy,
        label = L"H", linewidth = 2, color = :black,
    )
    plot!(
        p2, data.t, data.total_zollner_residual,
        label = L"\sum \Delta U^{\mathrm Z}\ (\mathrm{emergent\ gravity})",
        linewidth = 1.8, color = :firebrick, linestyle = :dash,
    )

    return plot(
        p1, p2;
        layout = grid(2, 1),
        size = _multi_panel_size(2),
    )
end

"""
    plot_zollner_force_residual(data::PairForceData) -> Plot

Two-panel visualization of the Zöllner extra force for one particle pair.

Panel 1: Total force magnitude alongside the Zöllner extra magnitude |（κ-1)·F_Coulomb|.
Panel 2: Ratio of extra Zöllner force to total force — shows relative size of mismatch.

Useful for verifying that `a` is small (ratio ≪ 1) while the emergent gravity is present.
"""
function WeberElectrodynamics.plot_zollner_force_residual(data::PairForceData)::Plots.Plot
    i, j = data.pair
    kappa_str = @sprintf("%.4g", data.kappa)

    # Panel 1: Total magnitude vs Zöllner extra
    p1 = plot(;
        title = "Force Magnitude — Pair ($i,$j), κ=$kappa_str",
        xlabel = "",
        ylabel = L"\lVert \mathbf F \rVert",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(
        p1, data.t, data.magnitude,
        label = L"\lVert \mathbf F_{\mathrm{total}} \rVert",
        linewidth = 2, color = :steelblue,
    )
    plot!(
        p1, data.t, data.zollner_extra_magnitude,
        label = L"\lVert (\kappa-1)\, \mathbf F_{\mathrm{C}} \rVert",
        linewidth = 1.5, color = :firebrick, linestyle = :dash,
    )

    # Panel 2: Ratio (κ-1)*F_coulomb / F_total
    ratio = ifelse.(
        data.magnitude .> eps(Float64),
        data.zollner_extra_magnitude ./ data.magnitude,
        zeros(length(data.magnitude)),
    )

    p2 = plot(;
        title = "Relative Zöllner Contribution",
        xlabel = L"t",
        ylabel = L"\lVert \Delta \mathbf F^{\mathrm Z} \rVert / \lVert \mathbf F_{\mathrm{total}} \rVert",
        legend = :outertopright,
        PLOT_DEFAULTS...,
    )
    plot!(
        p2, data.t, ratio,
        label = L"(\kappa-1)\ \mathrm{fraction}",
        linewidth = 1.5, color = :firebrick,
    )
    delta_kappa_str = @sprintf("%.4g", abs(data.kappa - 1.0))
    hline!(p2, [abs(data.kappa - 1.0)];
        label = latexstring("a = ", delta_kappa_str),
        linewidth = 1.0, color = :black, linestyle = :dot)

    return plot(
        p1, p2;
        layout = grid(2, 1),
        size = _multi_panel_size(2),
    )
end

"""
    plot_weber_vs_zollner(sol1, sol2; labels=["Weber", "Zöllner"]) -> Plot

Overlay trajectories from two solutions (same initial conditions, different κ) to
visualise the orbital divergence introduced by the Zöllner mismatch.

Both solutions must have the same number of particles and dimensions.
Solid lines = `sol1`, dashed lines = `sol2`.
"""
function WeberElectrodynamics.plot_weber_vs_zollner(
    sol1::WeberSolution,
    sol2::WeberSolution;
    labels::Vector{String} = ["Weber", "Zöllner"],
)::Plots.Plot
    n1 = sol1.prob.system.n_particles
    n2 = sol2.prob.system.n_particles
    dims = sol1.prob.system.dims
    @assert n1 == n2 "Both solutions must have the same number of particles"
    @assert dims == sol2.prob.system.dims "Both solutions must have the same dimensions"
    @assert dims in (2, 3) "Trajectory comparison only supported for 2D and 3D"

    traj1 = compute_trajectory_data(sol1, n1, dims)
    traj2 = compute_trajectory_data(sol2, n2, dims)

    colors = [:steelblue, :firebrick, :forestgreen, :darkorchid, :darkorange]

    plt = plot(;
        title = "Trajectory Comparison: $(labels[1]) vs $(labels[2])",
        xlabel = L"q_x",
        ylabel = L"q_y",
        legend = :outertopright,
        aspect_ratio = :equal,
        PLOT_DEFAULTS...,
        size = _square_size(),
    )

    if dims == 2
        for particle = 1:n1
            c = colors[mod1(particle, length(colors))]
            plot!(
                plt,
                traj1.trajectories[particle][:, 1],
                traj1.trajectories[particle][:, 2];
                label = latexstring("P_{", particle, "}\\ \\mathrm{", labels[1], "}"),
                linewidth = 1.5, color = c, linestyle = :solid,
            )
            plot!(
                plt,
                traj2.trajectories[particle][:, 1],
                traj2.trajectories[particle][:, 2];
                label = latexstring("P_{", particle, "}\\ \\mathrm{", labels[2], "}"),
                linewidth = 1.5, color = c, linestyle = :dash,
            )
        end
    else  # 3D
        for particle = 1:n1
            c = colors[mod1(particle, length(colors))]
            plot!(
                plt,
                traj1.trajectories[particle][:, 1],
                traj1.trajectories[particle][:, 2],
                traj1.trajectories[particle][:, 3];
                label = latexstring("P_{", particle, "}\\ \\mathrm{", labels[1], "}"),
                linewidth = 1.5, color = c, linestyle = :solid,
            )
            plot!(
                plt,
                traj2.trajectories[particle][:, 1],
                traj2.trajectories[particle][:, 2],
                traj2.trajectories[particle][:, 3];
                label = latexstring("P_{", particle, "}\\ \\mathrm{", labels[2], "}"),
                linewidth = 1.5, color = c, linestyle = :dash,
            )
        end
    end

    return plt
end

"""
    plot_zollner_phase_space(data1, data2; labels=["Weber","Zöllner"]) -> Plot

Overlaid (r, ṙ) phase portrait for the same pair from two different simulations.

Shows how the Zöllner mismatch shifts the phase-space orbit — e.g., orbits that are open
in standard Weber may curve inward under Zöllner gravity.
Solid = `data1`, dashed = `data2`.
"""
function WeberElectrodynamics.plot_zollner_phase_space(
    data1::PairForceData,
    data2::PairForceData;
    labels::Vector{String} = ["Weber", "Zöllner"],
)::Plots.Plot
    i, j = data1.pair
    plt = plot(;
        title = "Phase Space (r, ṙ) — Pair ($i,$j): $(labels[1]) vs $(labels[2])",
        xlabel = L"r\ (\mathrm{separation})",
        ylabel = L"\dot r\ (\mathrm{radial\ velocity})",
        legend = :outertopright,
        PLOT_DEFAULTS...,
        size = _square_size(),
    )

    plot!(
        plt,
        data1.phase_space.separation_distance,
        data1.phase_space.radial_velocity;
        label = labels[1],
        linewidth = 1.5, color = :steelblue, linestyle = :solid,
    )
    scatter!(
        plt,
        [data1.phase_space.separation_distance[1]],
        [data1.phase_space.radial_velocity[1]];
        marker = :circle, markersize = 6, color = :steelblue, label = "",
    )

    plot!(
        plt,
        data2.phase_space.separation_distance,
        data2.phase_space.radial_velocity;
        label = labels[2],
        linewidth = 1.5, color = :firebrick, linestyle = :dash,
    )
    scatter!(
        plt,
        [data2.phase_space.separation_distance[1]],
        [data2.phase_space.radial_velocity[1]];
        marker = :circle, markersize = 6, color = :firebrick, label = "",
    )

    return plt
end

end # module
