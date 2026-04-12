module WeberElectrodynamicsPythonPlotExt

using WeberElectrodynamics
using WeberElectrodynamics: @sprintf, norm, dot
using WeberElectrodynamics: TrajectoryData, EnergyData, PairForceData, MomentumData,
    WeberSolution, compute_trajectory_data
using PythonPlot

const plt = PythonPlot
const mpl = PythonPlot.matplotlib

# =============================================================================
# Publication figure sizes (inches — matplotlib's native unit)
# =============================================================================

const SINGLE_COLUMN_WIDTH_IN = 3.5   # inches
const DOUBLE_COLUMN_WIDTH_IN = 7.0   # inches
const SINGLE_PANEL_ASPECT = 0.75     # height/width for single-panel layouts
const MULTI_PANEL_ROW = 0.5          # height per row in multi-panel layouts

_single_panel_size() = (SINGLE_COLUMN_WIDTH_IN, SINGLE_COLUMN_WIDTH_IN * SINGLE_PANEL_ASPECT)
_multi_panel_size(n::Int) = (SINGLE_COLUMN_WIDTH_IN, SINGLE_COLUMN_WIDTH_IN * MULTI_PANEL_ROW * n)
_square_size() = (SINGLE_COLUMN_WIDTH_IN, SINGLE_COLUMN_WIDTH_IN)

# =============================================================================
# Matplotlib style
# =============================================================================

const _STYLE_APPLIED = Ref(false)

"""
    _apply_style!(; usetex=false, force=false)

Set matplotlib rcParams to a publication-quality default. Idempotent — called
lazily by each plot function on first use. Pass `force=true` to reapply.
"""
function _apply_style!(; usetex::Bool = false, force::Bool = false)
    if _STYLE_APPLIED[] && !force
        return nothing
    end
    rc = mpl.rcParams
    rc["figure.dpi"] = 300
    rc["savefig.dpi"] = 300
    rc["figure.autolayout"] = true

    rc["font.family"] = "serif"
    rc["font.serif"] = ["DejaVu Serif", "Bitstream Vera Serif", "Computer Modern Roman"]
    rc["mathtext.fontset"] = "cm"
    rc["text.usetex"] = usetex
    rc["font.size"] = 10
    rc["axes.titlesize"] = 11
    rc["axes.labelsize"] = 10
    rc["xtick.labelsize"] = 9
    rc["ytick.labelsize"] = 9
    rc["legend.fontsize"] = 8
    rc["figure.titlesize"] = 12

    rc["axes.linewidth"] = 0.8
    rc["axes.spines.top"] = false
    rc["axes.spines.right"] = false
    rc["axes.grid"] = true
    rc["grid.alpha"] = 0.3
    rc["grid.linestyle"] = "--"
    rc["grid.linewidth"] = 0.5

    rc["xtick.direction"] = "in"
    rc["ytick.direction"] = "in"
    rc["xtick.major.size"] = 4
    rc["ytick.major.size"] = 4
    rc["xtick.minor.size"] = 2
    rc["ytick.minor.size"] = 2
    rc["xtick.minor.visible"] = true
    rc["ytick.minor.visible"] = true

    rc["lines.linewidth"] = 1.5
    rc["lines.markersize"] = 4

    rc["legend.frameon"] = false
    rc["legend.borderaxespad"] = 0.4

    _STYLE_APPLIED[] = true
    return nothing
end

"""
    set_matplotlib_style!(; usetex=false)

Reapply the publication rcParams, optionally enabling LaTeX rendering
(requires a working TeX installation on the host).
"""
set_matplotlib_style!(; usetex::Bool = false) = _apply_style!(; usetex = usetex, force = true)

# =============================================================================
# Helper functions
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

"""Place legend outside the axes on the top-right (matches Plots' :outertopright)."""
function _outer_legend!(ax)
    ax.legend(; loc = "upper left", bbox_to_anchor = (1.02, 1.0), borderaxespad = 0.0)
end

# Color palette (matches Plots extension symbol palette)
const CLR_BLACK = "black"
const CLR_BLUE = "steelblue"
const CLR_RED = "firebrick"
const CLR_GREEN = "forestgreen"
const CLR_GRAY = "gray"
const CLR_PURPLE = "darkorchid"
const CLR_ORANGE = "darkorange"
const ZOLLNER_PALETTE = [CLR_BLUE, CLR_RED, CLR_GREEN, CLR_PURPLE, CLR_ORANGE]

# =============================================================================
# plot_trajectories
# =============================================================================

function WeberElectrodynamics.plot_trajectories(data::TrajectoryData)
    _apply_style!()
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

function _plot_trajectories_1d(data::TrajectoryData)
    fig, ax = plt.subplots(; figsize = _single_panel_size())
    cmap = mpl.pyplot.get_cmap("tab10")

    for particle = 1:data.n_particles
        color = cmap(mod(particle - 1, 10))
        xs = 1:size(data.trajectories[particle], 1)
        ax.plot(xs, data.trajectories[particle][:, 1];
            label = "Particle $particle", linewidth = 1.5, color = color)

        initial_label = particle == 1 ? "Initial" : nothing
        final_label = particle == 1 ? "Final" : nothing
        ax.scatter([1], [data.initial_positions[particle][1]];
            marker = "o", s = 36, color = color,
            edgecolors = "black", linewidths = 0.5, zorder = 5,
            label = initial_label)
        ax.scatter([length(xs)], [data.final_positions[particle][1]];
            marker = "s", s = 36, color = color,
            edgecolors = "black", linewidths = 0.5, zorder = 5,
            label = final_label)
    end

    ax.set_title("1D Particle Trajectories")
    ax.set_xlabel("Time Index")
    ax.set_ylabel(raw"$x$")
    _outer_legend!(ax)
    fig.tight_layout()
    return fig
end

function _plot_trajectories_2d(data::TrajectoryData)
    all_x = vcat([data.trajectories[p][:, 1] for p = 1:data.n_particles]...)
    all_y = vcat([data.trajectories[p][:, 2] for p = 1:data.n_particles]...)
    xmin, xmax = extrema(all_x)
    ymin, ymax = extrema(all_y)
    dx = xmax - xmin
    dy = ymax - ymin
    margin = 0.1 * max(dx, dy, eps())

    fig, ax = plt.subplots(; figsize = _square_size())
    cmap = mpl.pyplot.get_cmap("tab10")

    for particle = 1:data.n_particles
        color = cmap(mod(particle - 1, 10))
        ax.plot(
            data.trajectories[particle][:, 1],
            data.trajectories[particle][:, 2];
            label = "Particle $particle", linewidth = 1.5, color = color,
        )
        initial_label = particle == 1 ? "Initial" : nothing
        final_label = particle == 1 ? "Final" : nothing
        ax.scatter(
            [data.initial_positions[particle][1]],
            [data.initial_positions[particle][2]];
            marker = "o", s = 36, color = color,
            edgecolors = "black", linewidths = 0.5, zorder = 5,
            label = initial_label,
        )
        ax.scatter(
            [data.final_positions[particle][1]],
            [data.final_positions[particle][2]];
            marker = "s", s = 36, color = color,
            edgecolors = "black", linewidths = 0.5, zorder = 5,
            label = final_label,
        )
    end

    ax.set_title("2D Particle Trajectories")
    ax.set_xlabel(raw"$x$")
    ax.set_ylabel(raw"$y$")
    ax.set_xlim(xmin - margin, xmax + margin)
    ax.set_ylim(ymin - margin, ymax + margin)
    ax.set_aspect("equal")
    _outer_legend!(ax)
    fig.tight_layout()
    return fig
end

function _plot_trajectories_3d(data::TrajectoryData)
    w, _ = _single_panel_size()
    fig = plt.figure(; figsize = (w, w * 0.85))
    ax = fig.add_subplot(111; projection = "3d")
    cmap = mpl.pyplot.get_cmap("tab10")

    for particle = 1:data.n_particles
        color = cmap(mod(particle - 1, 10))
        ax.plot(
            data.trajectories[particle][:, 1],
            data.trajectories[particle][:, 2],
            data.trajectories[particle][:, 3];
            label = "Particle $particle", linewidth = 1.5, color = color,
        )
        initial_label = particle == 1 ? "Initial" : nothing
        final_label = particle == 1 ? "Final" : nothing
        ax.scatter(
            [data.initial_positions[particle][1]],
            [data.initial_positions[particle][2]],
            [data.initial_positions[particle][3]];
            marker = "o", s = 36, color = color,
            edgecolors = "black", linewidths = 0.5, label = initial_label,
        )
        ax.scatter(
            [data.final_positions[particle][1]],
            [data.final_positions[particle][2]],
            [data.final_positions[particle][3]];
            marker = "s", s = 36, color = color,
            edgecolors = "black", linewidths = 0.5, label = final_label,
        )
    end

    ax.set_title("3D Particle Trajectories")
    ax.set_xlabel(raw"$x$")
    ax.set_ylabel(raw"$y$")
    ax.set_zlabel(raw"$z$")
    ax.legend(; loc = "upper left", bbox_to_anchor = (1.05, 1.0))
    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_energy
# =============================================================================

function WeberElectrodynamics.plot_energy(data::EnergyData)
    _apply_style!()
    E0 = data.total_energy[1]
    relative_error = abs.((data.total_energy .- E0) ./ E0)
    relative_error = max.(relative_error, eps(Float64))

    fig, axs = plt.subplots(2, 1;
        figsize = _multi_panel_size(2),
        gridspec_kw = Dict("height_ratios" => [0.6, 0.4]))
    ax1 = axs[0]
    ax2 = axs[1]

    ax1.plot(data.t, data.total_energy; label = "Total E", linewidth = 2, color = CLR_BLACK)
    ax1.plot(data.t, data.kinetic_energy; label = "Kinetic T", linewidth = 1.5, color = CLR_BLUE)
    ax1.plot(data.t, data.total_potential_energy; label = "Potential U", linewidth = 1.5, color = CLR_RED)
    ax1.set_title("Energy Conservation")
    ax1.set_ylabel("Energy")
    _outer_legend!(ax1)

    max_err = maximum(relative_error)
    avg_err = length(relative_error) > 1 ?
        sum(relative_error[2:end]) / (length(relative_error) - 1) : max_err
    max_idx = argmax(relative_error)

    ax2.plot(data.t, relative_error; linewidth = 1.5, color = CLR_BLACK)
    ax2.set_yscale("log")
    ax2.axhline(max_err; linestyle = "--", linewidth = 1, color = CLR_RED,
        label = "max = $(_format_scientific(max_err))")
    ax2.axhline(avg_err; linestyle = ":", linewidth = 1, color = CLR_GRAY,
        label = "avg = $(_format_scientific(avg_err))")
    ax2.axvline(data.t[max_idx]; linestyle = "--", color = CLR_GRAY, alpha = 0.5)
    ax2.set_title("Relative Energy Error")
    ax2.set_xlabel(raw"Time $t$")
    ax2.set_ylabel(raw"$|\Delta E / E_0|$")
    _outer_legend!(ax2)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_pair_energy
# =============================================================================

function WeberElectrodynamics.plot_pair_energy(data::EnergyData, pair::Tuple{Int,Int})
    _apply_style!()
    if !haskey(data.pair_energies, pair)
        available = sort(collect(keys(data.pair_energies)))
        throw(ArgumentError("Pair $pair not found. Available pairs: $available"))
    end
    pair_data = data.pair_energies[pair]
    i, j = pair

    fig, axs = plt.subplots(2, 1;
        figsize = _multi_panel_size(2),
        gridspec_kw = Dict("height_ratios" => [0.6, 0.4]))
    ax1 = axs[0]
    ax2 = axs[1]

    ax1.plot(data.t, pair_data.coulomb_term; label = raw"Coulomb $q_i q_j / r$",
        linewidth = 1.5, color = CLR_RED)
    ax1.plot(data.t, pair_data.velocity_term; label = "Velocity term",
        linewidth = 1.5, color = CLR_BLUE)
    ax1.plot(data.t, pair_data.total_pair_potential; label = "Total",
        linewidth = 2, color = CLR_BLACK)
    ax1.set_title("Pair ($i,$j) Weber Potential")
    ax1.set_ylabel("Energy")
    _outer_legend!(ax1)

    rdot = pair_data.radial_velocity
    rdot_pos = [v >= 0 ? v : 0.0 for v in rdot]
    rdot_neg = [v < 0 ? v : 0.0 for v in rdot]

    ax2.fill_between(data.t, rdot_pos, 0; alpha = 0.3, color = CLR_RED, linewidth = 0,
        label = raw"Separating ($\dot r > 0$)")
    ax2.fill_between(data.t, rdot_neg, 0; alpha = 0.3, color = CLR_BLUE, linewidth = 0,
        label = raw"Approaching ($\dot r < 0$)")
    ax2.plot(data.t, rdot; label = raw"$\dot r$", linewidth = 1.5, color = CLR_BLACK)
    ax2.axhline(0.0; linestyle = "--", color = CLR_GRAY, linewidth = 1)
    ax2.set_title("Radial Dynamics")
    ax2.set_xlabel(raw"Time $t$")
    ax2.set_ylabel(raw"Radial Velocity $\dot r$")
    _outer_legend!(ax2)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_energy_errors
# =============================================================================

function WeberElectrodynamics.plot_energy_errors(data::EnergyData)
    _apply_style!()
    local_errors = _compute_local_errors(data.total_energy)
    local_errors_plot = max.(local_errors, eps(Float64))
    global_percent = _compute_global_error_percent(data.total_energy)
    global_percent_plot = max.(global_percent, eps(Float64))
    h_error = max.(data.hamiltonian_validation_error, eps(Float64))
    stats = data.statistics

    fig, axs = plt.subplots(3, 1;
        figsize = _multi_panel_size(3),
        gridspec_kw = Dict("height_ratios" => [0.35, 0.35, 0.30]))
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]

    ax1.plot(data.t[2:end], local_errors_plot; linewidth = 1, color = CLR_BLUE)
    ax1.set_yscale("log")
    ax1.axhline(stats.local_error_max; linestyle = "--", color = CLR_RED, linewidth = 1,
        label = "max = $(_format_scientific(stats.local_error_max))")
    ax1.axhline(stats.local_error_avg; linestyle = ":", color = CLR_GRAY, linewidth = 1,
        label = "avg = $(_format_scientific(stats.local_error_avg))")
    ax1.set_title("Local Energy Error")
    ax1.set_ylabel(raw"$|E_t - E_{t-1}|$")
    _outer_legend!(ax1)

    ax2.plot(data.t, global_percent_plot; linewidth = 1.5, color = CLR_RED)
    ax2.set_yscale("log")
    ax2.axhline(stats.global_error_percent_max; linestyle = "--", color = CLR_BLACK, linewidth = 1,
        label = "max = $(_format_scientific(stats.global_error_percent_max))%")
    ax2.set_title("Global Energy Drift")
    ax2.set_ylabel("Global Error (%)")
    _outer_legend!(ax2)

    max_h_err = maximum(data.hamiltonian_validation_error)
    ax3.plot(data.t, h_error; linewidth = 1, color = CLR_BLACK)
    ax3.set_yscale("log")
    ax3.axhline(max(max_h_err, eps(Float64)); linestyle = "--", color = CLR_GRAY, linewidth = 1,
        label = "max = $(_format_scientific(max_h_err))")
    ax3.set_title("Hamiltonian Validation")
    ax3.set_xlabel(raw"Time $t$")
    ax3.set_ylabel(raw"$|H_{\mathrm{computed}} - H_{\mathrm{compiled}}|$")
    _outer_legend!(ax3)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_pair_forces
# =============================================================================

function WeberElectrodynamics.plot_pair_forces(data::PairForceData)
    _apply_style!()
    n_times = length(data.t)
    dims = data.dims
    i, j = data.pair

    component_colors = [CLR_BLUE, CLR_RED, CLR_GREEN]
    component_labels = ["Fx", "Fy", "Fz"]
    decomp_colors = [CLR_BLACK, CLR_BLUE, CLR_RED, CLR_GREEN]

    w = SINGLE_COLUMN_WIDTH_IN
    fig, axs = plt.subplots(4, 1; figsize = (w, w * 2.0))
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    ax4 = axs[3]

    # Panel 1: force magnitude
    ax1.plot(data.t, data.magnitude; linewidth = 1.5, color = CLR_BLACK)
    stats = data.stats
    ax1.axhline(stats.max; linestyle = "--", linewidth = 1, color = CLR_RED,
        label = "max = $(_format_scientific(stats.max))")
    ax1.axhline(stats.min; linestyle = "--", linewidth = 1, color = CLR_BLUE,
        label = "min = $(_format_scientific(stats.min))")
    ax1.axhline(stats.mean; linestyle = ":", linewidth = 1, color = CLR_GRAY,
        label = "mean = $(_format_scientific(stats.mean))")
    ax1.set_title("Force Magnitude |F| — Pair ($i,$j)")
    ax1.set_ylabel(raw"$|F|$")
    _outer_legend!(ax1)

    # Panel 2: force components
    for d = 1:dims
        component = [data.force[t][d] for t = 1:n_times]
        ax2.plot(data.t, component;
            label = component_labels[d], linewidth = 1.5, color = component_colors[d])
    end
    ax2.axhline(0.0; linestyle = "--", color = CLR_GRAY, linewidth = 0.5)
    ax2.set_title("Force Components")
    ax2.set_ylabel("Force")
    _outer_legend!(ax2)

    # Panel 3: vector form decomposition
    k_sign = sign(data.charge_product)
    coulomb_signed = [k_sign * norm(data.coulomb[t]) for t = 1:n_times]
    vv_signed = [sign(dot(data.vector_term_vv[t], data.coulomb[t])) * k_sign *
                 norm(data.vector_term_vv[t]) for t = 1:n_times]
    ra_signed = [sign(dot(data.vector_term_ra[t], data.coulomb[t])) * k_sign *
                 norm(data.vector_term_ra[t]) for t = 1:n_times]
    rv2_signed = [sign(dot(data.vector_term_rv2[t], data.coulomb[t])) * k_sign *
                  norm(data.vector_term_rv2[t]) for t = 1:n_times]

    ax3.plot(data.t, coulomb_signed; label = "Coulomb", linewidth = 1.5, color = decomp_colors[1])
    ax3.plot(data.t, vv_signed; label = raw"$v \cdot v / c^2$",
        linewidth = 1.5, color = decomp_colors[2])
    ax3.plot(data.t, ra_signed; label = raw"$r \cdot a / c^2$",
        linewidth = 1.5, color = decomp_colors[3])
    ax3.plot(data.t, rv2_signed; label = raw"$-1.5(\hat r \cdot v)^2 / c^2$",
        linewidth = 1.5, color = decomp_colors[4])
    ax3.axhline(0.0; linestyle = "--", color = CLR_GRAY, linewidth = 0.5)
    ax3.set_title(raw"Vector Form: $F = F_C (1 + (v\cdot v + r\cdot a - 1.5(\hat r\cdot v)^2)/c^2)$")
    ax3.set_ylabel("Signed Force")
    _outer_legend!(ax3)

    # Panel 4: radial form decomposition
    rdot2_signed = [sign(dot(data.radial_term_rdot2[t], data.coulomb[t])) * k_sign *
                    norm(data.radial_term_rdot2[t]) for t = 1:n_times]
    rddot_signed = [sign(dot(data.radial_term_rddot[t], data.coulomb[t])) * k_sign *
                    norm(data.radial_term_rddot[t]) for t = 1:n_times]

    ax4.plot(data.t, coulomb_signed; label = "Coulomb", linewidth = 1.5, color = decomp_colors[1])
    ax4.plot(data.t, rdot2_signed; label = raw"$-\dot r^2/(2c^2)$",
        linewidth = 1.5, color = decomp_colors[2])
    ax4.plot(data.t, rddot_signed; label = raw"$r \cdot \ddot r / c^2$",
        linewidth = 1.5, color = decomp_colors[3])
    ax4.axhline(0.0; linestyle = "--", color = CLR_GRAY, linewidth = 0.5)
    ax4.set_title(raw"Radial Form: $F = F_C (1 - \dot r^2/(2c^2) + r\ddot r/c^2)$")
    ax4.set_xlabel(raw"Time $t$")
    ax4.set_ylabel("Signed Force")
    _outer_legend!(ax4)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_phase_space
# =============================================================================

function WeberElectrodynamics.plot_phase_space(data::PairForceData)
    _apply_style!()
    ps = data.phase_space
    i, j = data.pair

    fig, ax = plt.subplots(; figsize = _square_size())
    ax.plot(ps.separation_distance, ps.radial_velocity;
        linewidth = 1.5, color = CLR_BLACK)
    ax.scatter([ps.separation_distance[1]], [ps.radial_velocity[1]];
        marker = "o", s = 50, color = CLR_BLUE, edgecolors = "black",
        linewidths = 0.5, zorder = 5, label = "Initial")
    ax.scatter([ps.separation_distance[end]], [ps.radial_velocity[end]];
        marker = "s", s = 50, color = CLR_RED, edgecolors = "black",
        linewidths = 0.5, zorder = 5, label = "Final")

    ax.set_title("Phase Space Portrait — Pair ($i,$j)")
    ax.set_xlabel(raw"Separation $r$")
    ax.set_ylabel(raw"Radial Velocity $\dot r$")
    _outer_legend!(ax)
    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_momentum
# =============================================================================

function WeberElectrodynamics.plot_momentum(data::MomentumData)
    _apply_style!()
    dims = data.dims
    component_colors = [CLR_BLUE, CLR_RED, CLR_GREEN]
    component_labels = ["Px", "Py", "Pz"]

    one_panel = dims == 1 || isnothing(data.angular_momentum)

    if one_panel
        fig, ax1 = plt.subplots(; figsize = _single_panel_size())
    else
        fig, axs = plt.subplots(2, 1; figsize = _multi_panel_size(2))
        ax1 = axs[0]
        ax2 = axs[1]
    end

    for d = 1:dims
        ax1.plot(data.t, data.linear_momentum_components[:, d];
            label = component_labels[d], linewidth = 1.5, color = component_colors[d])
    end
    ax1.plot(data.t, data.linear_momentum_magnitude;
        label = raw"$|P|$", linewidth = 2, color = CLR_BLACK, linestyle = "--")
    ax1.axhline(0.0; linestyle = ":", color = CLR_GRAY, linewidth = 0.5)
    ax1.set_title("Total Linear Momentum")
    ax1.set_ylabel(raw"Momentum $P$")
    if one_panel
        ax1.set_xlabel(raw"Time $t$")
    end
    _outer_legend!(ax1)

    if one_panel
        fig.tight_layout()
        return fig
    end

    if dims == 2
        ax2.plot(data.t, data.angular_momentum;
            label = raw"$L_z$", linewidth = 1.5, color = CLR_BLACK)
        ax2.axhline(0.0; linestyle = ":", color = CLR_GRAY, linewidth = 0.5)
        ax2.set_title("Total Angular Momentum (z-component)")
        ax2.set_ylabel(raw"Angular Momentum $L_z$")
    else  # 3D
        L_labels = ["Lx", "Ly", "Lz"]
        for d = 1:3
            L_component = [data.angular_momentum[t][d] for t = 1:length(data.t)]
            ax2.plot(data.t, L_component;
                label = L_labels[d], linewidth = 1.5, color = component_colors[d])
        end
        ax2.plot(data.t, data.angular_momentum_magnitude;
            label = raw"$|L|$", linewidth = 2, color = CLR_BLACK, linestyle = "--")
        ax2.axhline(0.0; linestyle = ":", color = CLR_GRAY, linewidth = 0.5)
        ax2.set_title("Total Angular Momentum")
        ax2.set_ylabel(raw"Angular Momentum $L$")
    end
    ax2.set_xlabel(raw"Time $t$")
    _outer_legend!(ax2)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_zollner_energy
# =============================================================================

function WeberElectrodynamics.plot_zollner_energy(data::EnergyData)
    _apply_style!()

    fig, axs = plt.subplots(2, 1; figsize = _multi_panel_size(2))
    ax1 = axs[0]
    ax2 = axs[1]

    ax1.plot(data.t, data.total_potential_energy;
        label = raw"Total $V$", linewidth = 2, color = CLR_BLACK)
    for (pair_idx, ((i, j), pdata)) in enumerate(sort(collect(data.pair_energies); by = x -> x[1]))
        c = ZOLLNER_PALETTE[mod1(pair_idx, length(ZOLLNER_PALETTE))]
        kappa_str = @sprintf("%.4g", pdata.kappa)
        ax1.plot(data.t, pdata.total_pair_potential;
            label = "V($i,$j) κ=$kappa_str", linewidth = 1.2, color = c, linestyle = "-")
        ax1.plot(data.t, pdata.zollner_extra_potential;
            label = "ΔV_Z($i,$j)", linewidth = 1.0, color = c, linestyle = "--")
    end
    ax1.set_title("Potential Energy with Zöllner Contribution")
    ax1.set_ylabel("Energy")
    _outer_legend!(ax1)

    ax2.plot(data.t, data.total_energy;
        label = raw"Total $E$", linewidth = 2, color = CLR_BLACK)
    ax2.plot(data.t, data.total_zollner_residual;
        label = "Σ ΔV_Zöllner (emergent gravity)", linewidth = 1.8,
        color = CLR_RED, linestyle = "--")
    ax2.set_title("Zöllner Gravitational Residual")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Energy")
    _outer_legend!(ax2)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_zollner_force_residual
# =============================================================================

function WeberElectrodynamics.plot_zollner_force_residual(data::PairForceData)
    _apply_style!()
    i, j = data.pair
    kappa_str = @sprintf("%.4g", data.kappa)

    fig, axs = plt.subplots(2, 1; figsize = _multi_panel_size(2))
    ax1 = axs[0]
    ax2 = axs[1]

    ax1.plot(data.t, data.magnitude;
        label = raw"$|F_{\mathrm{total}}|$", linewidth = 2, color = CLR_BLUE)
    ax1.plot(data.t, data.zollner_extra_magnitude;
        label = raw"$|(\kappa - 1)\cdot F_C|$", linewidth = 1.5,
        color = CLR_RED, linestyle = "--")
    ax1.set_title("Force Magnitude — Pair ($i,$j), κ=$kappa_str")
    ax1.set_ylabel(raw"$|F|$")
    _outer_legend!(ax1)

    ratio = ifelse.(
        data.magnitude .> eps(Float64),
        data.zollner_extra_magnitude ./ data.magnitude,
        zeros(length(data.magnitude)),
    )
    ax2.plot(data.t, ratio;
        label = "(κ−1) fraction", linewidth = 1.5, color = CLR_RED)
    delta_kappa_str = @sprintf("%.4g", abs(data.kappa - 1.0))
    ax2.axhline(abs(data.kappa - 1.0);
        label = "a = $delta_kappa_str", linewidth = 1.0, color = CLR_BLACK, linestyle = ":")
    ax2.set_title("Relative Zöllner Contribution")
    ax2.set_xlabel("Time")
    ax2.set_ylabel(raw"$|\Delta F_{\mathrm{Zöllner}}| / |F_{\mathrm{total}}|$")
    _outer_legend!(ax2)

    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_weber_vs_zollner
# =============================================================================

function WeberElectrodynamics.plot_weber_vs_zollner(
    sol1::WeberSolution,
    sol2::WeberSolution;
    labels::Vector{String} = ["Weber", "Zöllner"],
)
    _apply_style!()
    n1 = sol1.prob.system.n_particles
    n2 = sol2.prob.system.n_particles
    dims = sol1.prob.system.dims
    @assert n1 == n2 "Both solutions must have the same number of particles"
    @assert dims == sol2.prob.system.dims "Both solutions must have the same dimensions"
    @assert dims in (2, 3) "Trajectory comparison only supported for 2D and 3D"

    traj1 = compute_trajectory_data(sol1, n1, dims)
    traj2 = compute_trajectory_data(sol2, n2, dims)

    cmap = mpl.pyplot.get_cmap("tab10")

    if dims == 2
        fig, ax = plt.subplots(; figsize = _square_size())
        for particle = 1:n1
            c = cmap(mod(particle - 1, 10))
            ax.plot(
                traj1.trajectories[particle][:, 1],
                traj1.trajectories[particle][:, 2];
                label = "P$particle $(labels[1])",
                linewidth = 1.5, color = c, linestyle = "-",
            )
            ax.plot(
                traj2.trajectories[particle][:, 1],
                traj2.trajectories[particle][:, 2];
                label = "P$particle $(labels[2])",
                linewidth = 1.5, color = c, linestyle = "--",
            )
        end
        ax.set_xlabel(raw"$x$")
        ax.set_ylabel(raw"$y$")
        ax.set_aspect("equal")
    else  # 3D
        w = SINGLE_COLUMN_WIDTH_IN
        fig = plt.figure(; figsize = (w, w))
        ax = fig.add_subplot(111; projection = "3d")
        for particle = 1:n1
            c = cmap(mod(particle - 1, 10))
            ax.plot(
                traj1.trajectories[particle][:, 1],
                traj1.trajectories[particle][:, 2],
                traj1.trajectories[particle][:, 3];
                label = "P$particle $(labels[1])",
                linewidth = 1.5, color = c, linestyle = "-",
            )
            ax.plot(
                traj2.trajectories[particle][:, 1],
                traj2.trajectories[particle][:, 2],
                traj2.trajectories[particle][:, 3];
                label = "P$particle $(labels[2])",
                linewidth = 1.5, color = c, linestyle = "--",
            )
        end
        ax.set_xlabel(raw"$x$")
        ax.set_ylabel(raw"$y$")
        ax.set_zlabel(raw"$z$")
    end

    ax.set_title("Trajectory Comparison: $(labels[1]) vs $(labels[2])")
    _outer_legend!(ax)
    fig.tight_layout()
    return fig
end

# =============================================================================
# plot_zollner_phase_space
# =============================================================================

function WeberElectrodynamics.plot_zollner_phase_space(
    data1::PairForceData,
    data2::PairForceData;
    labels::Vector{String} = ["Weber", "Zöllner"],
)
    _apply_style!()
    i, j = data1.pair

    fig, ax = plt.subplots(; figsize = _square_size())

    ax.plot(
        data1.phase_space.separation_distance,
        data1.phase_space.radial_velocity;
        label = labels[1], linewidth = 1.5, color = CLR_BLUE, linestyle = "-",
    )
    ax.scatter(
        [data1.phase_space.separation_distance[1]],
        [data1.phase_space.radial_velocity[1]];
        marker = "o", s = 50, color = CLR_BLUE,
        edgecolors = "black", linewidths = 0.5, zorder = 5,
    )

    ax.plot(
        data2.phase_space.separation_distance,
        data2.phase_space.radial_velocity;
        label = labels[2], linewidth = 1.5, color = CLR_RED, linestyle = "--",
    )
    ax.scatter(
        [data2.phase_space.separation_distance[1]],
        [data2.phase_space.radial_velocity[1]];
        marker = "o", s = 50, color = CLR_RED,
        edgecolors = "black", linewidths = 0.5, zorder = 5,
    )

    ax.set_title("Phase Space (r, ṙ) — Pair ($i,$j): $(labels[1]) vs $(labels[2])")
    ax.set_xlabel(raw"$r$ (separation)")
    ax.set_ylabel(raw"$\dot r$ (radial velocity)")
    _outer_legend!(ax)
    fig.tight_layout()
    return fig
end

end # module
