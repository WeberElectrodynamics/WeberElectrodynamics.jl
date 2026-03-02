module WeberElectrodynamicsMakieExt

using WeberElectrodynamics
using WeberElectrodynamics: @sprintf, norm, dot,
    WeberProblem, WeberIntegrator, WeberSystem, WeberSolution,
    SymmetricProjectionIntegrator, WeberAlgorithm,
    _pair_index, compute_total_kinetic_energy, compute_pair_weber_components
using CommonSolve: init, step!
using Makie

# =============================================================================
# Data Source Abstraction
# =============================================================================

abstract type AnimationDataSource end

mutable struct StreamingSource <: AnimationDataSource
    integrator::WeberIntegrator
    total_steps::Int
end

mutable struct ReplaySource <: AnimationDataSource
    sol::WeberSolution
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

    # Energy traces
    total_energy::Vector{Float64}
    kinetic_energy::Vector{Float64}
    potential_energy::Vector{Float64}

    # Linear momentum magnitude
    linear_momentum_mag::Vector{Float64}

    # Angular momentum (scalar for 2D)
    angular_momentum::Vector{Float64}

    # Per-pair phase space
    pair_separation::Dict{Tuple{Int,Int},Vector{Float64}}
    pair_radial_velocity::Dict{Tuple{Int,Int},Vector{Float64}}

    # Per-particle phase space (q, p components)
    particle_q::Dict{Int,Matrix{Float64}}
    particle_p::Dict{Int,Matrix{Float64}}
end

function RollingBuffer(capacity::Int, prob::WeberProblem)
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
        Vector{Float64}(undef, capacity),
        Vector{Float64}(undef, capacity),
        Vector{Float64}(undef, capacity),
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

function _compute_step_energy(q::AbstractVector{Float64}, p::AbstractVector{Float64}, prob::WeberProblem)
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
    return KE, PE, KE + PE
end

function _compute_step_linear_momentum(p::AbstractVector{Float64}, dims::Int, n::Int)
    P = zeros(dims)
    @inbounds for i in 1:n
        base = (i - 1) * dims
        for d in 1:dims
            P[d] += p[base+d]
        end
    end
    return norm(P)
end

function _compute_step_angular_momentum_2d(q::AbstractVector{Float64}, p::AbstractVector{Float64}, n::Int)
    L = 0.0
    @inbounds for i in 1:n
        x_idx = (i - 1) * 2 + 1
        y_idx = (i - 1) * 2 + 2
        L += q[x_idx] * p[y_idx] - q[y_idx] * p[x_idx]
    end
    return L
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
                    p::AbstractVector{Float64}, prob::WeberProblem)
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

    # Energy
    KE, PE, E = _compute_step_energy(q, p, prob)
    buf.kinetic_energy[idx] = KE
    buf.potential_energy[idx] = PE
    buf.total_energy[idx] = E

    # Linear momentum
    buf.linear_momentum_mag[idx] = _compute_step_linear_momentum(p, dims, n)

    # Angular momentum (2D)
    if dims >= 2
        buf.angular_momentum[idx] = _compute_step_angular_momentum_2d(q, p, n)
    else
        buf.angular_momentum[idx] = 0.0
    end

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

function _trim_to_window(data::Vector{Float64}, window::Int)
    n = length(data)
    n <= window ? data : data[end-window+1:end]
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
    prob::WeberProblem
    buffer::RollingBuffer
    is_playing::Observable{Bool}
    timer::Union{Nothing,Timer}
    tail_length::Observable{Int}
    display_window::Observable{Int}
    compute_batch::Observable{Int}
    buffer_size::Int

    # Reference values for error computation
    E0::Float64
    P0_mag::Float64
    L0::Float64

    # Phase space selection
    phase_selection::Observable{String}
    phase_component::Observable{Int}

    # For streaming reset
    alg::SymmetricProjectionIntegrator
    extended_prob::Union{Nothing,WeberProblem}

    # Axes for limit reset
    axes::Vector{Axis}
end

# =============================================================================
# Error Computation
# =============================================================================

function _compute_errors(state::AnimationState)
    buf = state.buffer
    if buf.count < 2
        return fill("--", 6)
    end

    # Current and previous indices (in circular buffer)
    curr_idx = buf.cursor == 1 ? buf.capacity : buf.cursor - 1
    prev_idx = curr_idx == 1 ? buf.capacity : curr_idx - 1
    if buf.count < buf.capacity && curr_idx <= 1
        return fill("--", 6)
    end

    E_curr = buf.total_energy[curr_idx]
    E_prev = buf.total_energy[prev_idx]
    P_curr = buf.linear_momentum_mag[curr_idx]
    P_prev = buf.linear_momentum_mag[prev_idx]
    L_curr = buf.angular_momentum[curr_idx]
    L_prev = buf.angular_momentum[prev_idx]

    E0 = state.E0
    P0 = state.P0_mag
    L0 = state.L0

    e_local = @sprintf("%.2e", abs(E_curr - E_prev))
    e_global = abs(E0) > eps() ? @sprintf("%.2e%%", abs((E_curr - E0) / E0) * 100) : @sprintf("%.2e", abs(E_curr - E0))
    p_local = @sprintf("%.2e", abs(P_curr - P_prev))
    p_global = @sprintf("%.2e", abs(P_curr - P0))
    l_local = @sprintf("%.2e", abs(L_curr - L_prev))
    l_global = abs(L0) > eps() ? @sprintf("%.2e%%", abs((L_curr - L0) / L0) * 100) : @sprintf("%.2e", abs(L_curr - L0))

    return [e_local, e_global, p_local, p_global, l_local, l_global]
end

# =============================================================================
# Observable Update
# =============================================================================

struct PlotObservables
    # Trajectory panel
    traj_x::Vector{Observable{Vector{Float64}}}
    traj_y::Vector{Observable{Vector{Float64}}}
    marker_pos::Vector{Observable{Point2f}}

    # Energy panel (left axis)
    energy_t::Observable{Vector{Float64}}
    total_energy::Observable{Vector{Float64}}
    kinetic_energy::Observable{Vector{Float64}}
    potential_energy::Observable{Vector{Float64}}

    # Momentum panel (right axis on energy panel)
    momentum_mag::Observable{Vector{Float64}}

    # Angular momentum panel
    angular_t::Observable{Vector{Float64}}
    angular_momentum::Observable{Vector{Float64}}

    # Phase space panel
    phase_x::Observable{Vector{Float64}}
    phase_y::Observable{Vector{Float64}}
    phase_marker::Observable{Point2f}

    # Error labels
    error_texts::Vector{Observable{String}}

    # Info display
    time_text::Observable{String}
    step_text::Observable{String}
end

function _update_observables!(state::AnimationState, obs::PlotObservables)
    buf = state.buffer
    prob = state.prob
    n = prob.system.n_particles

    # Linearize time
    t_lin = _linearize(buf.t, buf)
    dw = state.display_window[]

    # Trajectory tails (last tail_length entries)
    tail_len = min(state.tail_length[], buf.count)
    for particle in 1:n
        x_full = _linearize_col(buf.positions[particle], 1, buf)
        y_full = _linearize_col(buf.positions[particle], 2, buf)
        start_idx = max(1, length(x_full) - tail_len + 1)
        obs.traj_x[particle][] = x_full[start_idx:end]
        obs.traj_y[particle][] = y_full[start_idx:end]
        if !isempty(x_full)
            obs.marker_pos[particle][] = Point2f(x_full[end], y_full[end])
        end
    end

    # Energy (trimmed to display window)
    obs.energy_t[] = _trim_to_window(t_lin, dw)
    obs.total_energy[] = _trim_to_window(_linearize(buf.total_energy, buf), dw)
    obs.kinetic_energy[] = _trim_to_window(_linearize(buf.kinetic_energy, buf), dw)
    obs.potential_energy[] = _trim_to_window(_linearize(buf.potential_energy, buf), dw)

    # Momentum (trimmed to display window)
    obs.momentum_mag[] = _trim_to_window(_linearize(buf.linear_momentum_mag, buf), dw)

    # Angular momentum (trimmed to display window)
    obs.angular_t[] = _trim_to_window(t_lin, dw)
    obs.angular_momentum[] = _trim_to_window(_linearize(buf.angular_momentum, buf), dw)

    # Phase space
    _update_phase_space!(state, obs, buf)

    # Errors
    errs = _compute_errors(state)
    for i in 1:6
        obs.error_texts[i][] = errs[i]
    end

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

    # Update marker to current point
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
            # One-way scaling for x: only expand, never shrink
            cur_xl, cur_xr = _get_xlims(ax)
            new_xl = min(new_xl, cur_xl)
            new_xr = max(new_xr, cur_xr)
        end

        # One-way scaling for y: only expand, never shrink
        cur_yb, cur_yt = _get_ylims(ax)
        new_yb = min(new_yb, cur_yb)
        new_yt = max(new_yt, cur_yt)

        xlims!(ax, new_xl, new_xr)
        ylims!(ax, new_yb, new_yt)
    end
    return nothing
end

function _autoscale_dual!(ax_left::Axis, ax_right::Axis, x_obs::Observable,
                          y_left_obs::Observable, y_right_obs::Observable;
                          padding::Float64 = 0.05, track_x::Bool = false)
    on(y_left_obs) do yd
        xd = x_obs[]
        if isempty(xd) || isempty(yd) || length(xd) < 2
            return
        end
        xmin, xmax = extrema(xd)
        dx = xmax - xmin
        if dx < eps(Float64)
            dx = 1.0
        end
        new_xl = xmin - padding * dx
        new_xr = xmax + padding * dx

        if !track_x
            # One-way x scaling (shared)
            cur_xl, cur_xr = _get_xlims(ax_left)
            new_xl = min(new_xl, cur_xl)
            new_xr = max(new_xr, cur_xr)
        end
        xlims!(ax_left, new_xl, new_xr)

        # Left y axis — one-way
        ymin, ymax = extrema(yd)
        dy = ymax - ymin
        if dy < eps(Float64)
            dy = max(abs(ymin) * 0.1, eps(Float64))
        end
        new_yb = ymin - padding * dy
        new_yt = ymax + padding * dy

        cur_yb, cur_yt = _get_ylims(ax_left)
        new_yb = min(new_yb, cur_yb)
        new_yt = max(new_yt, cur_yt)
        ylims!(ax_left, new_yb, new_yt)

        # Right axis — one-way
        yr = y_right_obs[]
        if !isempty(yr)
            yrmin, yrmax = extrema(yr)
            dyr = yrmax - yrmin
            if dyr < eps(Float64)
                dyr = max(abs(yrmin) * 0.1, eps(Float64))
            end
            new_rb = yrmin - padding * dyr
            new_rt = yrmax + padding * dyr

            cur_rb, cur_rt = _get_ylims(ax_right)
            new_rb = min(new_rb, cur_rb)
            new_rt = max(new_rt, cur_rt)
            ylims!(ax_right, new_rb, new_rt)
        end
    end
    return nothing
end

function _autoscale_trajectory!(ax::Axis, traj_xs::Vector{Observable{Vector{Float64}}},
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
        dx = xmax - xmin
        dy = ymax - ymin
        span = max(dx, dy)
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

        # One-way scaling: only expand, never shrink
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

# =============================================================================
# Figure Construction
# =============================================================================

const PARTICLE_COLORS = [:steelblue, :firebrick, :forestgreen, :darkorange, :purple,
                         :teal, :crimson, :darkgoldenrod]

function _build_figure(state::AnimationState, obs::PlotObservables;
                       figure_size::Tuple{Int,Int} = (1400, 900))
    fig = Figure(; size = figure_size, backgroundcolor = :white,
        fonts = Attributes(
            regular = "TeX Gyre Heros Makie",
            bold = "TeX Gyre Heros Makie Bold",
            italic = "TeX Gyre Heros Makie Italic",
            bold_italic = "TeX Gyre Heros Makie Bold Italic",
            mono = "TeX Gyre Heros Makie",
        ))

    prob = state.prob
    n = prob.system.n_particles

    # =========================================================================
    # Top-left: Trajectory panel (2D)
    # =========================================================================
    ax_traj = Axis(fig[1:2, 1];
        title = "Particle Trajectories",
        xlabel = "x", ylabel = "y",
    )

    for particle in 1:n
        col = PARTICLE_COLORS[mod1(particle, length(PARTICLE_COLORS))]
        lines!(ax_traj, obs.traj_x[particle], obs.traj_y[particle];
            color = col, linewidth = 1.5, label = "P$particle")
        scatter!(ax_traj, obs.marker_pos[particle];
            color = col, markersize = 12)
    end
    axislegend(ax_traj; position = :rt, framevisible = false, labelsize = 10)

    # =========================================================================
    # Top-right: Energy (left axis) + |P| (right axis)
    # =========================================================================
    ax_energy = Axis(fig[1:2, 2];
        title = "Energy & Momentum",
        xlabel = "t", ylabel = "Energy",
    )
    ax_momentum = Axis(fig[1:2, 2];
        ylabel = "|P|",
        yaxisposition = :right,
        yticklabelcolor = :purple,
        ylabelcolor = :purple,
    )
    hidespines!(ax_momentum)
    hidexdecorations!(ax_momentum)
    linkxaxes!(ax_energy, ax_momentum)

    lines!(ax_energy, obs.energy_t, obs.total_energy;
        color = :black, linewidth = 2, label = "Total E")
    lines!(ax_energy, obs.energy_t, obs.kinetic_energy;
        color = :steelblue, linewidth = 1.5, label = "Kinetic T")
    lines!(ax_energy, obs.energy_t, obs.potential_energy;
        color = :firebrick, linewidth = 1.5, label = "Potential U")

    lines!(ax_momentum, obs.energy_t, obs.momentum_mag;
        color = :purple, linewidth = 1.5, linestyle = :dash, label = "|P|")

    axislegend(ax_energy; position = :lt, framevisible = false, labelsize = 10)

    # =========================================================================
    # Bottom-left: Angular Momentum
    # =========================================================================
    ax_angular = Axis(fig[3:4, 1];
        title = "Angular Momentum",
        xlabel = "t", ylabel = "Lz",
    )

    lines!(ax_angular, obs.angular_t, obs.angular_momentum;
        color = :black, linewidth = 1.5)
    hlines!(ax_angular, [state.L0]; color = :gray, linestyle = :dash, linewidth = 0.5)

    # =========================================================================
    # Bottom-right: Phase Space
    # =========================================================================
    ax_phase = Axis(fig[3:4, 2];
        title = _phase_title(state),
        xlabel = _phase_xlabel(state),
        ylabel = _phase_ylabel(state),
    )

    lines!(ax_phase, obs.phase_x, obs.phase_y;
        color = :black, linewidth = 1.0)
    scatter!(ax_phase, obs.phase_marker;
        color = :firebrick, markersize = 8)

    # Phase space selector options (menu placed in bottom row below)
    pairs = [(i, j) for i in 1:n for j in (i+1):n]
    pair_labels = ["Pair ($i,$j)" for (i, j) in pairs]
    particle_labels = ["Particle $i" for i in 1:n]
    all_labels = vcat(pair_labels, particle_labels)

    # =========================================================================
    # Bottom row: Controls + Sliders | Phase menu + Error display
    # =========================================================================
    controls_grid = fig[5, 1] = GridLayout()
    info_grid = fig[5, 2] = GridLayout()
    rowsize!(fig.layout, 5, Fixed(120))

    # --- Left side: Playback controls (sub-row 1) + Sliders (sub-row 2) ---
    play_label = @lift($(state.is_playing) ? "|| Pause" : "> Play")
    play_btn = Button(controls_grid[1, 1]; label = play_label, width = 90)
    reset_btn = Button(controls_grid[1, 2]; label = "Reset", width = 90)
    Label(controls_grid[1, 3], obs.time_text; fontsize = 11)
    Label(controls_grid[1, 4], obs.step_text; fontsize = 11)

    on(play_btn.clicks) do _
        state.is_playing[] = !state.is_playing[]
    end

    on(reset_btn.clicks) do _
        _reset_animation!(state, obs)
    end

    # Sliders
    Label(controls_grid[2, 1], "Trail:"; fontsize = 11, halign = :right)
    trail_slider = Slider(controls_grid[2, 2]; range = 10:10:state.buffer_size,
                          startvalue = state.tail_length[])
    Label(controls_grid[2, 3], "Window:"; fontsize = 11, halign = :right)
    window_slider = Slider(controls_grid[2, 4]; range = 50:10:state.buffer_size,
                           startvalue = state.display_window[])
    Label(controls_grid[2, 5], "Speed:"; fontsize = 11, halign = :right)
    speed_range = _log_linear_range(1000)
    speed_slider = Slider(controls_grid[2, 6]; range = speed_range,
                          startvalue = state.compute_batch[])

    on(trail_slider.value) do val
        state.tail_length[] = val
    end

    on(window_slider.value) do val
        state.display_window[] = val
    end

    on(speed_slider.value) do val
        state.compute_batch[] = val
    end

    # --- Right side: Phase menu (sub-row 1) + Error display (2 rows) ---
    # Phase space selector
    phase_grid = info_grid[1, 1] = GridLayout()
    Label(phase_grid[1, 1], "Phase:"; fontsize = 12, halign = :right)
    menu = Menu(phase_grid[1, 2]; options = all_labels, default = state.phase_selection[],
                width = 200)

    on(menu.selection) do sel
        state.phase_selection[] = sel
        ax_phase.title[] = _phase_title(state)
        ax_phase.xlabel[] = _phase_xlabel(state)
        ax_phase.ylabel[] = _phase_ylabel(state)
        _update_phase_space!(state, obs, state.buffer)
    end

    # Error display — 2 rows x 3 pairs, fontsize 12
    error_grid = info_grid[1:2, 2] = GridLayout()
    error_labels = ["E local:", "E global:", "|P| local:", "|P| global:", "Lz local:", "Lz global:"]
    for idx in 1:6
        row = idx <= 3 ? 1 : 2
        col_in_row = idx <= 3 ? idx : idx - 3
        col_base = (col_in_row - 1) * 2 + 1
        Label(error_grid[row, col_base], error_labels[idx]; fontsize = 12, halign = :right)
        Label(error_grid[row, col_base + 1], obs.error_texts[idx];
              fontsize = 12, halign = :left, color = :firebrick)
    end

    # =========================================================================
    # Auto-scaling setup
    # =========================================================================
    _autoscale_trajectory!(ax_traj, obs.traj_x, obs.traj_y)
    _autoscale_dual!(ax_energy, ax_momentum, obs.energy_t, obs.total_energy, obs.momentum_mag;
                     track_x = true)
    _autoscale!(ax_angular, obs.angular_t, obs.angular_momentum; track_x = true)
    _autoscale!(ax_phase, obs.phase_x, obs.phase_y)

    # Reset y-axis one-way expansion when display window changes
    on(state.display_window) do _
        for ax in [ax_energy, ax_momentum, ax_angular]
            autolimits!(ax)
        end
    end

    # Store axes so _reset_animation! can clear one-way limits
    state.axes = [ax_traj, ax_energy, ax_momentum, ax_angular, ax_phase]

    return fig
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
# Reset
# =============================================================================

function _reset_animation!(state::AnimationState, obs::PlotObservables)
    state.is_playing[] = false
    reset!(state.buffer)

    # Clear one-way axis limits so they re-expand from initial state
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

# Reset the integrator so it can keep stepping, reusing its pre-allocated
# history arrays.  Only the bookkeeping fields are touched — the physical
# state (t, q, p) and regularization hysteresis are preserved.
function _recycle_integrator!(integ::WeberIntegrator)
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

function _make_extended_problem(prob::WeberProblem)
    max_time = prob.tspan[1] + 1000 * prob.dt
    reg = prob.regularization
    zol = prob.zollner

    return WeberProblem(
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
        regularization_enabled = reg.enabled,
        regularization_r_on = reg.r_on,
        regularization_r_off = reg.r_off,
        regularization_r_on_factor = reg.r_on_factor,
        regularization_r_off_factor = reg.r_off_factor,
        regularization_max_substeps = reg.max_substeps,
        regularization_constraint_tolerance = reg.constraint_tolerance,
        regularization_g_floor = reg.g_floor,
        regularization_chain_enabled = reg.chain_enabled,
        regularization_backend = reg.backend,
        regularization_warn_on_fallback = reg.warn_on_fallback,
        regularization_collision_bounce_radius = reg.collision_bounce_radius,
        zollner_enabled = zol.enabled,
        zollner_a = zol.a,
    )
end

# =============================================================================
# Public API: Streaming Mode
# =============================================================================

function WeberElectrodynamics.animate_weber(
    prob::WeberProblem;
    buffer_size::Int = 2000,
    tail_length::Int = 200,
    compute_batch::Int = 1,
    initial_pair::Tuple{Int,Int} = (1, min(2, prob.system.n_particles)),
    phase_mode::Symbol = :pair,
    initial_particle::Int = 1,
    initial_component::Int = 1,
    figure_size::Tuple{Int,Int} = (1400, 900),
    alg::SymmetricProjectionIntegrator = SymmetricProjectionIntegrator(),
)
    @assert prob.system.dims >= 2 "Animation viewer requires at least 2D (got $(prob.system.dims)D)"
    @assert buffer_size > 0 "buffer_size must be positive"
    @assert tail_length > 0 "tail_length must be positive"
    @assert compute_batch > 0 "compute_batch must be positive"
    @assert tail_length <= buffer_size "tail_length must be <= buffer_size"

    # Create extended-time problem for streaming
    extended_prob = _make_extended_problem(prob)
    integrator = init(extended_prob, alg)

    # Compute initial reference values
    q0, p0 = prob.q_initial, prob.p_initial
    _, _, E0 = _compute_step_energy(q0, p0, prob)
    P0_mag = _compute_step_linear_momentum(p0, prob.system.dims, prob.system.n_particles)
    L0 = prob.system.dims >= 2 ? _compute_step_angular_momentum_2d(q0, p0, prob.system.n_particles) : 0.0

    # Phase selection
    phase_sel = phase_mode == :pair ?
        "Pair ($(initial_pair[1]),$(initial_pair[2]))" :
        "Particle $initial_particle"

    # Create buffer and state
    buffer = RollingBuffer(buffer_size, prob)
    source = StreamingSource(integrator, 0)

    state = AnimationState(
        source, prob, buffer,
        Observable(false), nothing,
        Observable(tail_length), Observable(buffer_size), Observable(compute_batch),
        buffer_size,
        E0, P0_mag, L0,
        Observable(phase_sel), Observable(initial_component),
        alg, extended_prob,
        Axis[],
    )

    # Push initial state
    push_step!(buffer, integrator.t, integrator.q, integrator.p, prob)

    # Create observables and figure
    obs = _create_observables(prob)
    _update_observables!(state, obs)

    fig = _build_figure(state, obs; figure_size = figure_size)
    _start_animation!(state, obs)

    # Force a native window — disable inline so the backend opens a real screen
    Makie.inline!(false)
    screen = display(fig)

    return screen
end

# =============================================================================
# Public API: Replay Mode
# =============================================================================

function WeberElectrodynamics.animate_weber(
    sol::WeberSolution;
    buffer_size::Int = 2000,
    tail_length::Int = 200,
    compute_batch::Int = 1,
    initial_pair::Tuple{Int,Int} = (1, min(2, sol.prob.system.n_particles)),
    phase_mode::Symbol = :pair,
    initial_particle::Int = 1,
    initial_component::Int = 1,
    figure_size::Tuple{Int,Int} = (1400, 900),
    stride::Int = 1,
)
    prob = sol.prob
    @assert prob.system.dims >= 2 "Animation viewer requires at least 2D (got $(prob.system.dims)D)"
    @assert buffer_size > 0 "buffer_size must be positive"
    @assert tail_length > 0 "tail_length must be positive"
    @assert compute_batch > 0 "compute_batch must be positive"
    @assert tail_length <= buffer_size "tail_length must be <= buffer_size"
    @assert stride > 0 "stride must be positive"

    # Compute initial reference values
    q0, p0 = sol.q[1], sol.p[1]
    _, _, E0 = _compute_step_energy(q0, p0, prob)
    P0_mag = _compute_step_linear_momentum(p0, prob.system.dims, prob.system.n_particles)
    L0 = prob.system.dims >= 2 ? _compute_step_angular_momentum_2d(q0, p0, prob.system.n_particles) : 0.0

    # Phase selection
    phase_sel = phase_mode == :pair ?
        "Pair ($(initial_pair[1]),$(initial_pair[2]))" :
        "Particle $initial_particle"

    # Create buffer and state
    buffer = RollingBuffer(buffer_size, prob)
    source = ReplaySource(sol, stride, 1)

    state = AnimationState(
        source, prob, buffer,
        Observable(false), nothing,
        Observable(tail_length), Observable(buffer_size), Observable(compute_batch),
        buffer_size,
        E0, P0_mag, L0,
        Observable(phase_sel), Observable(initial_component),
        SymmetricProjectionIntegrator(), nothing,
        Axis[],
    )

    # Push initial state
    push_step!(buffer, sol.t[1], sol.q[1], sol.p[1], prob)

    # Create observables and figure
    obs = _create_observables(prob)
    _update_observables!(state, obs)

    fig = _build_figure(state, obs; figure_size = figure_size)
    _start_animation!(state, obs)

    # Force a native window — disable inline so the backend opens a real screen
    Makie.inline!(false)
    screen = display(fig)

    return screen
end

# =============================================================================
# Observable Factory
# =============================================================================

function _create_observables(prob::WeberProblem)
    n = prob.system.n_particles

    traj_x = [Observable(Float64[]) for _ in 1:n]
    traj_y = [Observable(Float64[]) for _ in 1:n]
    marker_pos = [Observable(Point2f(0, 0)) for _ in 1:n]

    energy_t = Observable(Float64[])
    total_energy = Observable(Float64[])
    kinetic_energy = Observable(Float64[])
    potential_energy = Observable(Float64[])

    momentum_mag = Observable(Float64[])

    angular_t = Observable(Float64[])
    angular_momentum = Observable(Float64[])

    phase_x = Observable(Float64[])
    phase_y = Observable(Float64[])
    phase_marker = Observable(Point2f(0, 0))

    error_texts = [Observable("--") for _ in 1:6]

    time_text = Observable("t = 0.0000")
    step_text = Observable("step = 0")

    PlotObservables(
        traj_x, traj_y, marker_pos,
        energy_t, total_energy, kinetic_energy, potential_energy,
        momentum_mag,
        angular_t, angular_momentum,
        phase_x, phase_y, phase_marker,
        error_texts,
        time_text, step_text,
    )
end

end # module
