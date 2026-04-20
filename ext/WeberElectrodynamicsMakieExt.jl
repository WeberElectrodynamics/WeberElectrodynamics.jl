module WeberElectrodynamicsMakieExt

using WeberElectrodynamics
using WeberElectrodynamics: @sprintf, norm, dot,
    HamiltonianProblem, HamiltonianIntegrator, HamiltonianSystem, HamiltonianSolution,
    SymmetricProjectionIntegrator, HamiltonianAlgorithm,
    _pair_index, compute_total_kinetic_energy, compute_pair_weber_components
using CommonSolve: init, step!
using Makie

# =============================================================================
# Data Source Abstraction
# =============================================================================

abstract type AnimationDataSource end

mutable struct StreamingSource <: AnimationDataSource
    integrator::HamiltonianIntegrator
    total_steps::Int
end

mutable struct ReplaySource <: AnimationDataSource
    sol::HamiltonianSolution
    stride::Int
    cursor::Int
end

# =============================================================================
# Rolling Buffer
# =============================================================================

mutable struct RollingBuffer
    capacity::Int
    count::Int
    cursor::Int

    t::Vector{Float64}

    # Per-particle positions: positions[particle][dim, slot]
    positions::Vector{Matrix{Float64}}

    # Energy trace (used for the live error display)
    total_energy::Vector{Float64}

    # Per-pair phase space
    pair_separation::Dict{Tuple{Int,Int},Vector{Float64}}
    pair_radial_velocity::Dict{Tuple{Int,Int},Vector{Float64}}

    # Per-particle phase space (q, p components)
    particle_q::Dict{Int,Matrix{Float64}}
    particle_p::Dict{Int,Matrix{Float64}}
end

function RollingBuffer(capacity::Int, prob::HamiltonianProblem)
    n = prob.system.n_particles
    dims = prob.system.dims

    positions = [Matrix{Float64}(undef, dims, capacity) for _ in 1:n]

    pair_sep = Dict{Tuple{Int,Int},Vector{Float64}}()
    pair_rdot = Dict{Tuple{Int,Int},Vector{Float64}}()
    for i in 1:n, j in (i+1):n
        pair_sep[(i, j)] = Vector{Float64}(undef, capacity)
        pair_rdot[(i, j)] = Vector{Float64}(undef, capacity)
    end

    particle_q = Dict(i => Matrix{Float64}(undef, dims, capacity) for i in 1:n)
    particle_p = Dict(i => Matrix{Float64}(undef, dims, capacity) for i in 1:n)

    RollingBuffer(
        capacity, 0, 1,
        Vector{Float64}(undef, capacity),
        positions,
        Vector{Float64}(undef, capacity),
        pair_sep, pair_rdot,
        particle_q, particle_p,
    )
end

function reset!(buf::RollingBuffer)
    buf.count = 0
    buf.cursor = 1
    return nothing
end

# =============================================================================
# Per-Step Computation
# =============================================================================

function _compute_step_energy(q::AbstractVector{Float64}, p::AbstractVector{Float64}, prob::HamiltonianProblem)
    masses = prob.masses
    charges = prob.charges
    c = prob.c
    dims = prob.system.dims
    n = prob.system.n_particles
    kappas = prob.kappas

    KE = compute_total_kinetic_energy(p, masses, dims)
    PE = 0.0
    for i in 1:n, j in (i+1):n
        kappa_ij = kappas[_pair_index(i, j, n)]
        coulomb, velocity, _, _ = compute_pair_weber_components(
            q, p, i, j, masses, charges, c, dims, kappa_ij,
        )
        PE += coulomb + velocity
    end
    return KE + PE
end

function _compute_step_pair_phase(
    q::AbstractVector{Float64}, p::AbstractVector{Float64},
    i::Int, j::Int, masses::AbstractVector{Float64}, dims::Int,
)
    qi_start = (i - 1) * dims
    qj_start = (j - 1) * dims
    mi, mj = masses[i], masses[j]

    r_sq = 0.0
    r_dot_v = 0.0
    @inbounds for d in 1:dims
        dq = q[qi_start+d] - q[qj_start+d]
        dv = p[qi_start+d] / mi - p[qj_start+d] / mj
        r_sq += dq^2
        r_dot_v += dq * dv
    end
    r = sqrt(r_sq)
    rdot = r_dot_v / r
    return r, rdot
end

# =============================================================================
# Buffer Push
# =============================================================================

function push_step!(buf::RollingBuffer, t::Float64, q::AbstractVector{Float64},
                    p::AbstractVector{Float64}, prob::HamiltonianProblem)
    idx = buf.cursor
    dims = prob.system.dims
    n = prob.system.n_particles

    buf.t[idx] = t

    # Positions per particle
    @inbounds for particle in 1:n
        base = (particle - 1) * dims
        for d in 1:dims
            buf.positions[particle][d, idx] = q[base+d]
        end
    end

    # Total energy (only used for live error display)
    buf.total_energy[idx] = _compute_step_energy(q, p, prob)

    # Pair phase space
    @inbounds for i in 1:n, j in (i+1):n
        r, rdot = _compute_step_pair_phase(q, p, i, j, prob.masses, dims)
        buf.pair_separation[(i, j)][idx] = r
        buf.pair_radial_velocity[(i, j)][idx] = rdot
    end

    # Per-particle phase space
    @inbounds for particle in 1:n
        base = (particle - 1) * dims
        for d in 1:dims
            buf.particle_q[particle][d, idx] = q[base+d]
            buf.particle_p[particle][d, idx] = p[base+d]
        end
    end

    # Advance cursor
    buf.cursor = (idx % buf.capacity) + 1
    buf.count = min(buf.count + 1, buf.capacity)

    return nothing
end

# =============================================================================
# Buffer Linearization
# =============================================================================

function _linearize(arr::Vector{Float64}, buf::RollingBuffer)
    if buf.count == 0
        return Float64[]
    end
    if buf.count < buf.capacity
        return arr[1:buf.count]
    end
    oldest = buf.cursor  # next write position = oldest entry
    return vcat(@view(arr[oldest:buf.capacity]), @view(arr[1:oldest-1]))
end

function _linearize_col(mat::Matrix{Float64}, row::Int, buf::RollingBuffer)
    if buf.count == 0
        return Float64[]
    end
    if buf.count < buf.capacity
        return vec(mat[row, 1:buf.count])
    end
    oldest = buf.cursor
    return vcat(@view(mat[row, oldest:buf.capacity]), @view(mat[row, 1:oldest-1]))
end

# Log-linear speed range: 1,2,...,9, 10,20,...,90, 100,200,...,900, 1000
function _log_linear_range(max_val::Int)
    values = Int[]
    decade = 1
    while decade <= max_val
        for k in 1:9
            v = k * decade
            v <= max_val && push!(values, v)
        end
        decade *= 10
    end
    return values
end

# =============================================================================
# Animation State
# =============================================================================

mutable struct AnimationState
    source::AnimationDataSource
    prob::HamiltonianProblem
    buffer::RollingBuffer
    is_playing::Observable{Bool}
    timer::Union{Nothing,Timer}
    tail_length::Observable{Int}
    compute_batch::Observable{Int}
    buffer_size::Int

    # Reference value for energy-error display
    E0::Float64

    # Phase space selection
    phase_selection::Observable{String}
    phase_component::Observable{Int}

    # For streaming reset
    alg::SymmetricProjectionIntegrator
    extended_prob::Union{Nothing,HamiltonianProblem}

    # Axes for limit reset (mixed Axis / Axis3)
    axes::Vector{Any}
end

# =============================================================================
# Energy-Error Computation
# =============================================================================

function _format_energy_error(state::AnimationState)
    buf = state.buffer
    if buf.count < 2
        return "ΔE/E₀ (max) = --"
    end
    E_lin = _linearize(buf.total_energy, buf)
    E0 = state.E0
    if abs(E0) < eps()
        err = maximum(abs.(E_lin .- E0))
        return @sprintf("|ΔE| (max) = %.3e", err)
    end
    err_pct = maximum(abs.(E_lin .- E0)) / abs(E0) * 100
    return @sprintf("ΔE/E₀ (max) = %.3e %%", err_pct)
end

# =============================================================================
# Plot Observables
# =============================================================================

struct PlotObservables
    # Trajectory (dims-aware)
    traj_x::Vector{Observable{Vector{Float64}}}
    traj_y::Vector{Observable{Vector{Float64}}}
    traj_z::Vector{Observable{Vector{Float64}}}     # empty when dims==2
    marker_2d::Vector{Observable{Point2f}}           # populated when dims==2
    marker_3d::Vector{Observable{Point3f}}           # populated when dims==3

    # Phase space
    phase_x::Observable{Vector{Float64}}
    phase_y::Observable{Vector{Float64}}
    phase_marker::Observable{Point2f}

    # Live error display
    energy_error_text::Observable{String}

    # Info display
    time_text::Observable{String}
    step_text::Observable{String}
end

function _create_observables(prob::HamiltonianProblem)
    n = prob.system.n_particles
    dims = prob.system.dims

    traj_x = [Observable(Float64[]) for _ in 1:n]
    traj_y = [Observable(Float64[]) for _ in 1:n]
    traj_z = [Observable(Float64[]) for _ in 1:n]  # always allocated, used only when dims==3

    if dims == 3
        marker_2d = Observable{Point2f}[]
        marker_3d = [Observable(Point3f(0, 0, 0)) for _ in 1:n]
    else
        marker_2d = [Observable(Point2f(0, 0)) for _ in 1:n]
        marker_3d = Observable{Point3f}[]
    end

    phase_x = Observable(Float64[])
    phase_y = Observable(Float64[])
    phase_marker = Observable(Point2f(0, 0))

    energy_error_text = Observable("ΔE/E₀ (max) = --")
    time_text = Observable("t = 0.0000")
    step_text = Observable("step = 0")

    PlotObservables(
        traj_x, traj_y, traj_z, marker_2d, marker_3d,
        phase_x, phase_y, phase_marker,
        energy_error_text,
        time_text, step_text,
    )
end

# =============================================================================
# Observable Update
# =============================================================================

function _update_observables!(state::AnimationState, obs::PlotObservables)
    buf = state.buffer
    prob = state.prob
    n = prob.system.n_particles
    dims = prob.system.dims

    # Trajectory tails (last tail_length entries)
    tail_len = min(state.tail_length[], buf.count)
    for particle in 1:n
        x_full = _linearize_col(buf.positions[particle], 1, buf)
        y_full = _linearize_col(buf.positions[particle], 2, buf)
        start_idx = max(1, length(x_full) - tail_len + 1)
        x_tail = x_full[start_idx:end]
        y_tail = y_full[start_idx:end]
        obs.traj_x[particle][] = x_tail
        obs.traj_y[particle][] = y_tail

        if dims == 3
            z_full = _linearize_col(buf.positions[particle], 3, buf)
            z_tail = z_full[start_idx:end]
            obs.traj_z[particle][] = z_tail
            if !isempty(x_full)
                obs.marker_3d[particle][] = Point3f(x_full[end], y_full[end], z_full[end])
            end
        else
            if !isempty(x_full)
                obs.marker_2d[particle][] = Point2f(x_full[end], y_full[end])
            end
        end
    end

    # Phase space
    _update_phase_space!(state, obs, buf)

    # Energy-error label
    obs.energy_error_text[] = _format_energy_error(state)

    # Info
    if state.source isa StreamingSource
        integ = state.source.integrator
        obs.time_text[] = @sprintf("t = %.4f", integ.t)
        obs.step_text[] = @sprintf("step = %d", state.source.total_steps)
    elseif state.source isa ReplaySource
        src = state.source
        t_val = buf.count > 0 ? buf.t[buf.cursor == 1 ? buf.capacity : buf.cursor - 1] : 0.0
        obs.time_text[] = @sprintf("t = %.4f", t_val)
        obs.step_text[] = @sprintf("frame = %d / %d", src.cursor - 1, length(src.sol.t))
    end

    return nothing
end

function _update_phase_space!(state::AnimationState, obs::PlotObservables, buf::RollingBuffer)
    sel = state.phase_selection[]

    if startswith(sel, "Pair")
        m = match(r"Pair \((\d+),(\d+)\)", sel)
        if isnothing(m)
            return
        end
        i, j = parse(Int, m[1]), parse(Int, m[2])
        key = i < j ? (i, j) : (j, i)
        if haskey(buf.pair_separation, key)
            obs.phase_x[] = _linearize(buf.pair_separation[key], buf)
            obs.phase_y[] = _linearize(buf.pair_radial_velocity[key], buf)
        end
    elseif startswith(sel, "Particle")
        m = match(r"Particle (\d+)", sel)
        if isnothing(m)
            return
        end
        particle = parse(Int, m[1])
        comp = state.phase_component[]
        if haskey(buf.particle_q, particle)
            obs.phase_x[] = _linearize_col(buf.particle_q[particle], comp, buf)
            obs.phase_y[] = _linearize_col(buf.particle_p[particle], comp, buf)
        end
    end

    px = obs.phase_x[]
    py = obs.phase_y[]
    if !isempty(px) && !isempty(py)
        obs.phase_marker[] = Point2f(px[end], py[end])
    end

    return nothing
end

# =============================================================================
# Auto-Scaling
# =============================================================================

function _get_xlims(ax::Axis)
    fl = ax.finallimits[]
    return (fl.origin[1], fl.origin[1] + fl.widths[1])
end

function _get_ylims(ax::Axis)
    fl = ax.finallimits[]
    return (fl.origin[2], fl.origin[2] + fl.widths[2])
end

function _autoscale!(ax::Axis, x_obs::Observable, y_obs::Observable;
                     padding::Float64 = 0.05, track_x::Bool = false)
    on(y_obs) do yd
        xd = x_obs[]
        if isempty(xd) || isempty(yd) || length(xd) < 2
            return
        end
        xmin, xmax = extrema(xd)
        ymin, ymax = extrema(yd)
        dx = xmax - xmin
        dy = ymax - ymin
        if dx < eps(Float64)
            dx = 1.0
        end
        if dy < eps(Float64)
            dy = max(abs(ymin) * 0.1, eps(Float64))
        end
        new_xl = xmin - padding * dx
        new_xr = xmax + padding * dx
        new_yb = ymin - padding * dy
        new_yt = ymax + padding * dy

        if !track_x
            cur_xl, cur_xr = _get_xlims(ax)
            new_xl = min(new_xl, cur_xl)
            new_xr = max(new_xr, cur_xr)
        end

        cur_yb, cur_yt = _get_ylims(ax)
        new_yb = min(new_yb, cur_yb)
        new_yt = max(new_yt, cur_yt)

        xlims!(ax, new_xl, new_xr)
        ylims!(ax, new_yb, new_yt)
    end
    return nothing
end

function _autoscale_trajectory_2d!(ax::Axis,
                                   traj_xs::Vector{Observable{Vector{Float64}}},
                                   traj_ys::Vector{Observable{Vector{Float64}}};
                                   padding::Float64 = 0.05)
    on(traj_xs[1]) do _
        all_x = Float64[]
        all_y = Float64[]
        for (xo, yo) in zip(traj_xs, traj_ys)
            append!(all_x, xo[])
            append!(all_y, yo[])
        end
        if isempty(all_x) || isempty(all_y)
            return
        end
        xmin, xmax = extrema(all_x)
        ymin, ymax = extrema(all_y)
        span = max(xmax - xmin, ymax - ymin)
        if span < eps(Float64)
            span = 1.0
        end
        cx = (xmin + xmax) / 2
        cy = (ymin + ymax) / 2
        half = span / 2 * (1 + padding)

        new_xl = cx - half
        new_xr = cx + half
        new_yb = cy - half
        new_yt = cy + half

        cur_xl, cur_xr = _get_xlims(ax)
        cur_yb, cur_yt = _get_ylims(ax)
        new_xl = min(new_xl, cur_xl)
        new_xr = max(new_xr, cur_xr)
        new_yb = min(new_yb, cur_yb)
        new_yt = max(new_yt, cur_yt)

        xlims!(ax, new_xl, new_xr)
        ylims!(ax, new_yb, new_yt)
    end
    return nothing
end

function _get_lims_3d(ax::Axis3)
    fl = ax.finallimits[]
    o = fl.origin
    w = fl.widths
    return (o[1], o[1] + w[1], o[2], o[2] + w[2], o[3], o[3] + w[3])
end

function _autoscale_trajectory_3d!(ax::Axis3,
                                   traj_xs::Vector{Observable{Vector{Float64}}},
                                   traj_ys::Vector{Observable{Vector{Float64}}},
                                   traj_zs::Vector{Observable{Vector{Float64}}};
                                   padding::Float64 = 0.05)
    on(traj_xs[1]) do _
        all_x = Float64[]
        all_y = Float64[]
        all_z = Float64[]
        for (xo, yo, zo) in zip(traj_xs, traj_ys, traj_zs)
            append!(all_x, xo[])
            append!(all_y, yo[])
            append!(all_z, zo[])
        end
        if isempty(all_x) || isempty(all_y) || isempty(all_z)
            return
        end
        xmin, xmax = extrema(all_x)
        ymin, ymax = extrema(all_y)
        zmin, zmax = extrema(all_z)
        span = max(xmax - xmin, ymax - ymin, zmax - zmin)
        if span < eps(Float64)
            span = 1.0
        end
        cx = (xmin + xmax) / 2
        cy = (ymin + ymax) / 2
        cz = (zmin + zmax) / 2
        half = span / 2 * (1 + padding)

        new_xl = cx - half
        new_xr = cx + half
        new_yb = cy - half
        new_yt = cy + half
        new_zl = cz - half
        new_zr = cz + half

        cur_xl, cur_xr, cur_yb, cur_yt, cur_zl, cur_zr = _get_lims_3d(ax)
        new_xl = min(new_xl, cur_xl)
        new_xr = max(new_xr, cur_xr)
        new_yb = min(new_yb, cur_yb)
        new_yt = max(new_yt, cur_yt)
        new_zl = min(new_zl, cur_zl)
        new_zr = max(new_zr, cur_zr)

        limits!(ax, new_xl, new_xr, new_yb, new_yt, new_zl, new_zr)
    end
    return nothing
end

# =============================================================================
# Phase Space Helpers
# =============================================================================

function _phase_title(state::AnimationState)
    sel = state.phase_selection[]
    if startswith(sel, "Pair")
        return "Phase Space (r, dr/dt) — $sel"
    else
        comp = state.phase_component[]
        coord = ["x", "y", "z"][comp]
        return "Phase Space ($coord, p_$coord) — $sel"
    end
end

function _phase_xlabel(state::AnimationState)
    sel = state.phase_selection[]
    if startswith(sel, "Pair")
        return "r"
    else
        comp = state.phase_component[]
        return ["x", "y", "z"][comp]
    end
end

function _phase_ylabel(state::AnimationState)
    sel = state.phase_selection[]
    if startswith(sel, "Pair")
        return "dr/dt"
    else
        comp = state.phase_component[]
        return "p_" * ["x", "y", "z"][comp]
    end
end

# =============================================================================
# Figure Construction
# =============================================================================

_particle_color(i::Integer) = (cs = Makie.wong_colors(); cs[mod1(i, length(cs))])

function _alpha_gradient(col, n::Integer; min_alpha = 0.15)
    n <= 0 && return RGBAf[]
    base = Makie.to_color(col)
    n == 1 && return [RGBAf(base.r, base.g, base.b, 1.0)]
    return [RGBAf(base.r, base.g, base.b,
                  min_alpha + (1.0 - min_alpha) * (k - 1) / (n - 1)) for k in 1:n]
end

function _weber_theme()
    base = theme_latexfonts()
    overlay = Theme(
        fontsize = 14,
        backgroundcolor = :white,
        palette = (color = Makie.wong_colors(),),
        Axis = (
            xgridcolor = (:black, 0.08), ygridcolor = (:black, 0.08),
            xminorgridvisible = false, yminorgridvisible = false,
            topspinevisible = false, rightspinevisible = false,
            xtickalign = 1, ytickalign = 1,
            titlesize = 15, xlabelsize = 13, ylabelsize = 13,
        ),
        Axis3 = (
            xgridcolor = (:black, 0.12),
            ygridcolor = (:black, 0.12),
            zgridcolor = (:black, 0.12),
            titlesize = 15, xlabelsize = 13, ylabelsize = 13, zlabelsize = 13,
            protrusions = 40,
        ),
        Label = (fontsize = 13,),
    )
    return merge(overlay, base)
end

function _build_figure(state::AnimationState, obs::PlotObservables;
                       figure_size::Tuple{Int,Int} = (1200, 800))
    return with_theme(_weber_theme()) do
        _build_figure_impl(state, obs, figure_size)
    end
end

function _build_figure_impl(state::AnimationState, obs::PlotObservables,
                            figure_size::Tuple{Int,Int})
    fig = Figure(; size = figure_size, px_per_unit = 2)

    prob = state.prob
    n = prob.system.n_particles
    dims = prob.system.dims

    # =========================================================================
    # Big trajectory panel (Axis3 for 3D, Axis for 2D)
    # =========================================================================
    particle_cols = [_particle_color(i) for i in 1:n]
    legend_entries = [[LineElement(color = particle_cols[i], linewidth = 2.2),
                       MarkerElement(marker = :circle, color = particle_cols[i],
                                     markersize = 10, strokewidth = 0.6,
                                     strokecolor = :black)] for i in 1:n]
    legend_labels = ["P$i" for i in 1:n]

    if dims == 3
        ax_traj = Axis3(fig[1:3, 1:3];
            title = "Particle Trajectories",
            xlabel = "x", ylabel = "y", zlabel = "z",
            aspect = :data,
            perspectiveness = 0.3,
            viewmode = :fit,
        )
        for particle in 1:n
            col = particle_cols[particle]
            xobs = obs.traj_x[particle]
            yobs = obs.traj_y[particle]
            zobs = obs.traj_z[particle]
            pts = lift(xobs, yobs, zobs) do x, y, z
                m = min(length(x), length(y), length(z))
                Point3f[Point3f(x[k], y[k], z[k]) for k in 1:m]
            end
            trail_colors = lift(pts) do p
                _alpha_gradient(col, length(p))
            end
            lines!(ax_traj, pts; color = trail_colors, linewidth = 2.2,
                linecap = :round, joinstyle = :round)
            scatter!(ax_traj, obs.marker_3d[particle];
                color = col, markersize = 14,
                strokewidth = 0.6, strokecolor = :black, fxaa = true)
        end
        axislegend(ax_traj, legend_entries, legend_labels;
            position = :rt, framevisible = true,
            framecolor = (:black, 0.15), labelsize = 11)
    else
        ax_traj = Axis(fig[1:3, 1:3];
            title = "Particle Trajectories",
            xlabel = "x", ylabel = "y",
            aspect = DataAspect(),
        )
        for particle in 1:n
            col = particle_cols[particle]
            xobs = obs.traj_x[particle]
            yobs = obs.traj_y[particle]
            trail_colors = lift(xobs, yobs) do x, y
                _alpha_gradient(col, min(length(x), length(y)))
            end
            lines!(ax_traj, xobs, yobs; color = trail_colors, linewidth = 2.2,
                linecap = :round, joinstyle = :round)
            scatter!(ax_traj, obs.marker_2d[particle];
                color = col, markersize = 14,
                strokewidth = 0.6, strokecolor = :black)
        end
        axislegend(ax_traj, legend_entries, legend_labels;
            position = :rt, framevisible = true,
            framecolor = (:black, 0.15), labelsize = 11)
    end

    # =========================================================================
    # Right sidebar: phase menu, phase plot, info column
    # =========================================================================
    pairs = [(i, j) for i in 1:n for j in (i+1):n]
    pair_labels = ["Pair ($i,$j)" for (i, j) in pairs]
    particle_labels = ["Particle $i" for i in 1:n]
    all_labels = vcat(pair_labels, particle_labels)

    # Phase menu (top of sidebar)
    menu_grid = fig[1, 4] = GridLayout()
    Label(menu_grid[1, 1], "Phase:"; fontsize = 11, halign = :right)
    menu = Menu(menu_grid[1, 2]; options = all_labels,
                default = state.phase_selection[], width = 160)

    # Phase space panel
    ax_phase = Axis(fig[2, 4];
        title = _phase_title(state),
        xlabel = _phase_xlabel(state),
        ylabel = _phase_ylabel(state),
    )
    lines!(ax_phase, obs.phase_x, obs.phase_y;
        color = (:black, 0.65), linewidth = 1.6, linecap = :round)
    scatter!(ax_phase, obs.phase_marker;
        color = :black, markersize = 11,
        strokewidth = 0.6, strokecolor = :white)

    # Info column: energy error + time + step
    info_grid = fig[3, 4] = GridLayout()
    Label(info_grid[1, 1], obs.energy_error_text;
        fontsize = 14, halign = :left, font = :bold)
    Label(info_grid[2, 1], obs.time_text;
        fontsize = 13, halign = :left)
    Label(info_grid[3, 1], obs.step_text;
        fontsize = 13, halign = :left)

    # Sidebar column width
    colsize!(fig.layout, 4, Relative(0.25))

    # =========================================================================
    # Bottom controls row
    # =========================================================================
    controls_grid = fig[4, 1:4] = GridLayout()
    rowsize!(fig.layout, 4, Fixed(60))

    play_label = @lift($(state.is_playing) ? "|| Pause" : "> Play")
    play_btn = Button(controls_grid[1, 1]; label = play_label, width = 90)
    reset_btn = Button(controls_grid[1, 2]; label = "Reset", width = 70)

    Label(controls_grid[1, 3], "Trail:"; fontsize = 11, halign = :right)
    trail_slider = Slider(controls_grid[1, 4]; range = 10:10:state.buffer_size,
                          startvalue = state.tail_length[])

    Label(controls_grid[1, 5], "Speed:"; fontsize = 11, halign = :right)
    speed_range = _log_linear_range(1000)
    speed_slider = Slider(controls_grid[1, 6]; range = speed_range,
                          startvalue = state.compute_batch[])

    on(play_btn.clicks) do _
        state.is_playing[] = !state.is_playing[]
    end

    on(reset_btn.clicks) do _
        _reset_animation!(state, obs)
    end

    on(trail_slider.value) do val
        state.tail_length[] = val
    end

    on(speed_slider.value) do val
        state.compute_batch[] = val
    end

    on(menu.selection) do sel
        state.phase_selection[] = sel
        ax_phase.title[] = _phase_title(state)
        ax_phase.xlabel[] = _phase_xlabel(state)
        ax_phase.ylabel[] = _phase_ylabel(state)
        _update_phase_space!(state, obs, state.buffer)
    end

    # =========================================================================
    # Auto-scaling setup
    # =========================================================================
    if dims == 3
        _autoscale_trajectory_3d!(ax_traj, obs.traj_x, obs.traj_y, obs.traj_z)
    else
        _autoscale_trajectory_2d!(ax_traj, obs.traj_x, obs.traj_y)
    end
    _autoscale!(ax_phase, obs.phase_x, obs.phase_y)

    state.axes = Any[ax_traj, ax_phase]

    return fig
end

# =============================================================================
# Reset
# =============================================================================

function _reset_animation!(state::AnimationState, obs::PlotObservables)
    state.is_playing[] = false
    reset!(state.buffer)

    for ax in state.axes
        autolimits!(ax)
    end

    if state.source isa StreamingSource
        prob_to_init = isnothing(state.extended_prob) ? state.prob : state.extended_prob
        state.source.integrator = init(prob_to_init, state.alg)
        state.source.total_steps = 0
        push_step!(state.buffer, state.source.integrator.t,
                   state.source.integrator.q, state.source.integrator.p, state.prob)
    elseif state.source isa ReplaySource
        state.source.cursor = 1
        sol = state.source.sol
        push_step!(state.buffer, sol.t[1], sol.q[1], sol.p[1], state.prob)
    end

    _update_observables!(state, obs)
    return nothing
end

# =============================================================================
# Integrator Recycling (Streaming)
# =============================================================================

function _recycle_integrator!(integ::HamiltonianIntegrator)
    span = integ.t_end - integ.t_history[1]
    integ.step_count = 0
    integ.t_end = integ.t + span
    integ.t_history[1] = integ.t
    integ.q_history[1] .= integ.q
    integ.p_history[1] .= integ.p
    fill!(integ.buffers.μ, 0.0)
    return nothing
end

# =============================================================================
# Compute Loop (Timer-Based)
# =============================================================================

function _advance_source!(state::AnimationState)
    if state.source isa StreamingSource
        integ = state.source.integrator
        more = step!(integ)
        if !more
            _recycle_integrator!(integ)
            more = step!(integ)
            !more && return false
        end
        state.source.total_steps += 1
        push_step!(state.buffer, integ.t, integ.q, integ.p, state.prob)
        return true
    elseif state.source isa ReplaySource
        src = state.source
        sol = src.sol
        next_cursor = src.cursor + src.stride
        if next_cursor > length(sol.t)
            return false
        end
        src.cursor = next_cursor
        push_step!(state.buffer, sol.t[next_cursor], sol.q[next_cursor], sol.p[next_cursor], state.prob)
        return true
    end
    return false
end

function _start_animation!(state::AnimationState, obs::PlotObservables)
    state.timer = Timer(0.0; interval = 1 / 60) do _
        if !state.is_playing[]
            return
        end

        batch = state.compute_batch[]
        for _ in 1:batch
            more = _advance_source!(state)
            if !more
                state.is_playing[] = false
                return
            end
        end

        _update_observables!(state, obs)
    end
    return nothing
end

# =============================================================================
# Extended Problem Construction (for Streaming)
# =============================================================================

function _make_extended_problem(prob::HamiltonianProblem)
    max_time = prob.tspan[1] + 1000 * prob.dt
    return HamiltonianProblem(
        prob.system,
        (prob.tspan[1], max_time),
        prob.q_initial,
        prob.p_initial;
        masses = prob.masses,
        charges = prob.charges,
        c = prob.c,
        dt = prob.dt,
        convergence_tolerance = prob.convergence_tolerance,
        maximum_iterations = prob.maximum_iterations,
        regularization = prob.regularization,
        zollner = prob.zollner,
    )
end

# =============================================================================
# Public API: Streaming Mode
# =============================================================================

"""
    animate_weber(prob::HamiltonianProblem; kwargs...) -> screen

Launch an interactive real-time animation of a Weber simulation.

Streams simulation data into a rolling buffer and displays a live dashboard
with one large trajectory panel (`Axis3` when `prob.system.dims == 3`,
otherwise `Axis`), a phase-space panel with selector dropdown, and a
live energy-error readout. Requires a windowed Makie backend (GLMakie or
WGLMakie recommended). Requires at least 2D (`prob.system.dims ≥ 2`).

# Keywords
- `buffer_size=2000`: Maximum timesteps retained in the rolling buffer.
- `tail_length=200`: Length of visible trajectory tail (must be ≤ `buffer_size`).
- `compute_batch=1`: Integration steps computed per animation frame.
- `initial_pair=(1,2)`: Default pair shown in the phase-space panel.
- `phase_mode=:pair`: Phase-space mode; `:pair` for pair separation portrait,
  `:particle` for single-particle phase space.
- `initial_particle=1`: Particle index for `:particle` phase-space mode.
- `initial_component=1`: Component index for particle phase-space mode.
- `figure_size=(1200,800)`: Window size in pixels `(width, height)`.
- `alg=SymmetricProjectionIntegrator()`: Integrator algorithm.

# Returns
- A Makie screen handle; close the window to stop the simulation.
"""
function WeberElectrodynamics.animate_weber(
    prob::HamiltonianProblem;
    buffer_size::Int = 2000,
    tail_length::Int = 200,
    compute_batch::Int = 1,
    initial_pair::Tuple{Int,Int} = (1, min(2, prob.system.n_particles)),
    phase_mode::Symbol = :pair,
    initial_particle::Int = 1,
    initial_component::Int = 1,
    figure_size::Tuple{Int,Int} = (1200, 800),
    alg::SymmetricProjectionIntegrator = SymmetricProjectionIntegrator(),
)
    @assert prob.system.dims >= 2 "Animation viewer requires at least 2D (got $(prob.system.dims)D)"
    @assert buffer_size > 0 "buffer_size must be positive"
    @assert tail_length > 0 "tail_length must be positive"
    @assert compute_batch > 0 "compute_batch must be positive"
    @assert tail_length <= buffer_size "tail_length must be <= buffer_size"

    extended_prob = _make_extended_problem(prob)
    integrator = init(extended_prob, alg)

    q0, p0 = prob.q_initial, prob.p_initial
    E0 = _compute_step_energy(q0, p0, prob)

    phase_sel = phase_mode == :pair ?
        "Pair ($(initial_pair[1]),$(initial_pair[2]))" :
        "Particle $initial_particle"

    buffer = RollingBuffer(buffer_size, prob)
    source = StreamingSource(integrator, 0)

    state = AnimationState(
        source, prob, buffer,
        Observable(false), nothing,
        Observable(tail_length), Observable(compute_batch),
        buffer_size,
        E0,
        Observable(phase_sel), Observable(initial_component),
        alg, extended_prob,
        Any[],
    )

    push_step!(buffer, integrator.t, integrator.q, integrator.p, prob)

    obs = _create_observables(prob)
    _update_observables!(state, obs)

    fig = _build_figure(state, obs; figure_size = figure_size)
    _start_animation!(state, obs)

    Makie.inline!(false)
    screen = display(fig)

    return screen
end

# =============================================================================
# Public API: Replay Mode
# =============================================================================

"""
    animate_weber(sol::HamiltonianSolution; kwargs...) -> screen

Replay a completed `HamiltonianSolution` as an interactive animated dashboard.

Identical dashboard layout to the streaming form (single dominant trajectory
panel, phase-space sidebar, live energy-error readout), but replays
pre-computed trajectory data. Requires at least 2D
(`sol.prob.system.dims ≥ 2`).

# Keywords
- `buffer_size=2000`: Rolling display buffer size.
- `tail_length=200`: Visible trajectory tail length (must be ≤ `buffer_size`).
- `compute_batch=1`: Replay steps advanced per animation frame.
- `initial_pair=(1,2)`: Default pair for the phase-space panel.
- `phase_mode=:pair`: `:pair` or `:particle` phase-space mode.
- `initial_particle=1`, `initial_component=1`: Phase-space selection for
  `:particle` mode.
- `figure_size=(1200,800)`: Window size in pixels `(width, height)`.
- `stride=1`: Skip every `stride`-th stored timestep during replay.

# Returns
- A Makie screen handle; close the window to stop replay.
"""
function WeberElectrodynamics.animate_weber(
    sol::HamiltonianSolution;
    buffer_size::Int = 2000,
    tail_length::Int = 200,
    compute_batch::Int = 1,
    initial_pair::Tuple{Int,Int} = (1, min(2, sol.prob.system.n_particles)),
    phase_mode::Symbol = :pair,
    initial_particle::Int = 1,
    initial_component::Int = 1,
    figure_size::Tuple{Int,Int} = (1200, 800),
    stride::Int = 1,
)
    prob = sol.prob
    @assert prob.system.dims >= 2 "Animation viewer requires at least 2D (got $(prob.system.dims)D)"
    @assert buffer_size > 0 "buffer_size must be positive"
    @assert tail_length > 0 "tail_length must be positive"
    @assert compute_batch > 0 "compute_batch must be positive"
    @assert tail_length <= buffer_size "tail_length must be <= buffer_size"
    @assert stride > 0 "stride must be positive"

    q0, p0 = sol.q[1], sol.p[1]
    E0 = _compute_step_energy(q0, p0, prob)

    phase_sel = phase_mode == :pair ?
        "Pair ($(initial_pair[1]),$(initial_pair[2]))" :
        "Particle $initial_particle"

    buffer = RollingBuffer(buffer_size, prob)
    source = ReplaySource(sol, stride, 1)

    state = AnimationState(
        source, prob, buffer,
        Observable(false), nothing,
        Observable(tail_length), Observable(compute_batch),
        buffer_size,
        E0,
        Observable(phase_sel), Observable(initial_component),
        SymmetricProjectionIntegrator(), nothing,
        Any[],
    )

    push_step!(buffer, sol.t[1], sol.q[1], sol.p[1], prob)

    obs = _create_observables(prob)
    _update_observables!(state, obs)

    fig = _build_figure(state, obs; figure_size = figure_size)
    _start_animation!(state, obs)

    Makie.inline!(false)
    screen = display(fig)

    return screen
end

end # module
