using LinearAlgebra

abstract type WeberAlgorithm end

const REG_MODE_NONE = UInt8(0)
const REG_MODE_PAIR = UInt8(1)
const REG_MODE_CHAIN = UInt8(2)

const REG_BACKEND_ADAPTIVE = :adaptive_cartesian
const REG_BACKEND_LIFTED = :lifted_pair
const REG_BACKEND_DISABLED = :disabled
const REG_BACKEND_MIXED = :mixed

"""
    SymmetricProjectionIntegrator(; relaxation=0.25)

Symplectic integrator based on Strang-splitting with symmetric projection in
extended phase space.

Implements the method of Jayawardana & Ohsawa (2021) for Weber's
velocity-dependent Hamiltonian. The projection step solves a fixed-point
iteration whose convergence is controlled by `relaxation`.

# Keywords
- `relaxation=0.25`: Fixed-point relaxation factor ω ∈ (0, 1]. Smaller values
  slow convergence but improve stability near close encounters.

# Fields
- `relaxation::Float64`: Stored relaxation factor.
"""
struct SymmetricProjectionIntegrator <: WeberAlgorithm
    relaxation::Float64

    function SymmetricProjectionIntegrator(; relaxation::Real = 0.25)
        @assert 0 < relaxation <= 1 "Relaxation must be in (0, 1], got $relaxation"
        new(Float64(relaxation))
    end
end

# =============================================================================
# Zöllner Extension (research/experimental — see docs/src/zollner.md)
# =============================================================================

"""
    ZollnerOptions(; enabled=false, a=0.0)

Configuration for the Zöllner electrogravitational extension to Weber's force law.

When enabled, unlike-sign charge pairs receive a coupling factor κ = 1 + a,
producing an emergent attractive correction to the Weber potential.
Like-sign pairs are unaffected (κ = 1).

# Keywords
- `enabled=false`: Activate the Zöllner mismatch.
- `a=0.0`: Mismatch parameter; must be positive when `enabled=true`.

# Fields
- `enabled::Bool`: Whether the extension is active.
- `a::Float64`: Zöllner mismatch parameter a.
"""
struct ZollnerOptions
    enabled::Bool
    a::Float64

    function ZollnerOptions(; enabled::Bool = false, a::Real = 0.0)
        if enabled
            @assert a > 0 "Zöllner mismatch parameter a must be positive when enabled, got $a"
        end
        new(enabled, Float64(a))
    end
end

# Compute per-pair Zöllner coupling factors κ_ij.
# κ_ij = 1+a for unlike-sign charge pairs, 1.0 for like-sign pairs.
# Returns a vector of length n_particles*(n_particles-1)/2, ordered by (i<j).
# Note: sign(0.0) == 0.0, so a neutral particle (q=0) is treated as unlike
# any charged particle and receives κ = 1+a when Zöllner is enabled.
function _compute_zollner_kappas(
    charges::Vector{Float64},
    zollner::ZollnerOptions,
    n_particles::Int,
)::Vector{Float64}
    n_pairs = n_particles * (n_particles - 1) ÷ 2
    kappas = ones(Float64, n_pairs)
    if zollner.enabled
        idx = 1
        @inbounds for i = 1:n_particles
            for j = (i+1):n_particles
                if sign(charges[i]) != sign(charges[j])
                    kappas[idx] = 1.0 + zollner.a
                end
                idx += 1
            end
        end
    end
    return kappas
end

# =============================================================================
# Regularization (optional, advanced — see docs/src/regularization.md)
# =============================================================================

"""
    RegularizationOptions(; kwargs...)

Configuration for close-encounter regularization.

When two particles approach within `r_on`, the integrator switches from
Cartesian coordinates to a regularized representation (Levi-Civita/KS or
adaptive Cartesian substeps) and back once the separation exceeds `r_off`.

# Keywords
- `enabled=false`: Enable regularization globally.
- `r_on=nothing`: Absolute activation distance. If `nothing`, computed as
  `r_on_factor × min_initial_separation`.
- `r_off=nothing`: Absolute deactivation distance. If `nothing`, computed as
  `r_off_factor × min_initial_separation`. Must be greater than `r_on`.
- `r_on_factor=0.15`: Scale factor for automatic `r_on` (positive).
- `r_off_factor=0.25`: Scale factor for automatic `r_off` (positive).
- `max_substeps=512`: Maximum regularized substeps per macro-step.
- `constraint_tolerance=1e-12`: Convergence threshold for the projection constraint.
- `g_floor=1e-12`: Minimum regularization scale to prevent division by zero.
- `chain_enabled=true`: Allow chain regularization for multi-particle close encounters.
- `backend=:lifted_pair`: Regularization backend. `:lifted_pair` uses
  Levi-Civita/KS (2D only; auto-falls back to `:adaptive_cartesian` for 3D).
  `:adaptive_cartesian` works for all dimensions.
- `warn_on_fallback=true`: Emit a warning when the backend is automatically changed.
- `collision_bounce_radius=0.0`: Reflect pairs that come closer than this distance
  before each macro-step (0.0 = disabled). Intended for head-on (ℓ=0) collisions.

# Fields
See keyword documentation above; each keyword maps directly to a stored field.
"""
struct RegularizationOptions
    enabled::Bool
    r_on::Union{Nothing,Float64}
    r_off::Union{Nothing,Float64}
    r_on_factor::Float64
    r_off_factor::Float64
    max_substeps::Int
    constraint_tolerance::Float64
    g_floor::Float64
    chain_enabled::Bool
    backend::Symbol
    warn_on_fallback::Bool
    collision_bounce_radius::Float64

    function RegularizationOptions(
        ;
        enabled::Bool = false,
        r_on::Union{Nothing,Real} = nothing,
        r_off::Union{Nothing,Real} = nothing,
        r_on_factor::Real = 0.15,
        r_off_factor::Real = 0.25,
        max_substeps::Integer = 512,
        constraint_tolerance::Real = 1e-12,
        g_floor::Real = 1e-12,
        chain_enabled::Bool = true,
        backend::Symbol = REG_BACKEND_LIFTED,
        warn_on_fallback::Bool = true,
        collision_bounce_radius::Real = 0.0,
    )
        if !isnothing(r_on)
            @assert r_on > 0 "regularization_r_on must be positive"
        end
        if !isnothing(r_off)
            @assert r_off > 0 "regularization_r_off must be positive"
        end
        @assert r_on_factor > 0 "regularization_r_on_factor must be positive"
        @assert r_off_factor > 0 "regularization_r_off_factor must be positive"
        @assert max_substeps > 0 "regularization_max_substeps must be positive"
        @assert constraint_tolerance > 0 "regularization_constraint_tolerance must be positive"
        @assert g_floor > 0 "regularization_g_floor must be positive"
        if !isnothing(r_on) && !isnothing(r_off)
            @assert r_off > r_on "regularization_r_off must be greater than regularization_r_on"
        end
        @assert backend in (REG_BACKEND_ADAPTIVE, REG_BACKEND_LIFTED) "regularization_backend must be :adaptive_cartesian or :lifted_pair"
        @assert collision_bounce_radius >= 0 "collision_bounce_radius must be non-negative"

        new(
            enabled,
            isnothing(r_on) ? nothing : Float64(r_on),
            isnothing(r_off) ? nothing : Float64(r_off),
            Float64(r_on_factor),
            Float64(r_off_factor),
            Int(max_substeps),
            Float64(constraint_tolerance),
            Float64(g_floor),
            chain_enabled,
            backend,
            warn_on_fallback,
            Float64(collision_bounce_radius),
        )
    end
end

"""
    RegularizationDiagnostics

Diagnostics collected during a `solve!` call describing regularization usage.

Returned as `WeberSolution.regularization`. All step counts refer to macro-steps
of the outer `SymmetricProjectionIntegrator`.

# Fields
- `enabled::Bool`: Whether regularization was active for this solve.
- `requested_backend::Symbol`: Backend requested in `RegularizationOptions`.
- `used_backend::Symbol`: Backend actually used (may differ due to fallback).
- `activation_count::Int`: Number of times regularization was switched on.
- `deactivation_count::Int`: Number of times regularization was switched off.
- `active_steps::Int`: Steps taken while any regularization was active.
- `pair_steps::Int`: Steps handled by the pair regularization path.
- `adaptive_pair_steps::Int`: Steps handled by `:adaptive_cartesian`.
- `lifted_pair_steps::Int`: Steps handled by `:lifted_pair` (Levi-Civita/KS).
- `chain_steps::Int`: Steps handled by the chain (multi-pair) path.
- `unregularized_steps::Int`: Steps taken in plain Cartesian coordinates.
- `backend_fallback_steps::Int`: Steps where a fallback backend was used.
- `total_substeps::Int`: Total regularized substeps across all macro-steps.
- `max_substeps_used::Int`: Largest substep count used in any single macro-step.
- `max_constraint_violation::Float64`: Maximum projection constraint residual observed.
- `min_encounter_distance::Float64`: Minimum pairwise separation seen during the solve.
- `mode_history::Vector{UInt8}`: Per-step regularization mode (0=none, 1=pair, 2=chain).
"""
mutable struct RegularizationDiagnostics
    enabled::Bool
    requested_backend::Symbol
    used_backend::Symbol
    activation_count::Int
    deactivation_count::Int
    active_steps::Int
    pair_steps::Int
    adaptive_pair_steps::Int
    lifted_pair_steps::Int
    chain_steps::Int
    unregularized_steps::Int
    backend_fallback_steps::Int
    total_substeps::Int
    max_substeps_used::Int
    max_constraint_violation::Float64
    min_encounter_distance::Float64
    mode_history::Vector{UInt8}

    function RegularizationDiagnostics(
        enabled::Bool,
        n_steps::Int,
        requested_backend::Symbol,
        used_backend::Symbol,
    )
        new(
            enabled,
            requested_backend,
            used_backend,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0.0,
            Inf,
            fill(REG_MODE_NONE, n_steps),
        )
    end
end

mutable struct RegularizationBuffers
    n_particles::Int
    dims::Int
    dof::Int
    n_pairs::Int

    pair_i::Vector{Int}
    pair_j::Vector{Int}
    pair_distance::Vector{Float64}

    adjacency::Matrix{Bool}
    visited::Vector{Bool}
    queue::Vector{Int}

    component_nodes::Vector{Int}
    component_mask::Vector{Bool}
    active_nodes::Vector{Int}
    active_count::Int

    is_active::Bool
    active_mode::UInt8
    active_anchor_i::Int
    active_anchor_j::Int

    r_on::Float64
    r_off::Float64

    effective_backend::Symbol
    backend_fallback::Bool

    chain_order::Vector{Int}
    chain_used::Vector{Bool}

    rel_q::Vector{Float64}
    rel_p::Vector{Float64}
    rel_qdot::Vector{Float64}
    rel_pdot::Vector{Float64}

    pair_R::Vector{Float64}
    pair_P::Vector{Float64}
    pair_R_mid::Vector{Float64}
    pair_P_mid::Vector{Float64}
    pair_Rdot::Vector{Float64}
    pair_Pdot::Vector{Float64}

    lc_u::Vector{Float64}
    lc_U::Vector{Float64}
    lc_u_mid::Vector{Float64}
    lc_U_mid::Vector{Float64}
    lc_du_tau::Vector{Float64}
    lc_dU_tau::Vector{Float64}
    lc_du_tau_mid::Vector{Float64}
    lc_dU_tau_mid::Vector{Float64}

    temp_rel_q::Vector{Float64}
    temp_rel_p::Vector{Float64}

    ks_u::Vector{Float64}
    ks_U::Vector{Float64}
    ks_J::Matrix{Float64}
    ks_n::Vector{Float64}

    params_pair::Vector{Float64}
    dq_full::Vector{Float64}
    dp_full::Vector{Float64}
    dq_pair::Vector{Float64}
    dp_pair::Vector{Float64}
    dq_ext::Vector{Float64}
    dp_ext::Vector{Float64}
    q_mid::Vector{Float64}
    p_mid::Vector{Float64}

    function RegularizationBuffers(
        n_particles::Int,
        dims::Int,
        dof::Int,
        r_on::Float64,
        r_off::Float64,
        effective_backend::Symbol,
        backend_fallback::Bool,
    )
        n_pairs = n_particles * (n_particles - 1) ÷ 2
        pair_i = Vector{Int}(undef, n_pairs)
        pair_j = Vector{Int}(undef, n_pairs)
        idx = 1
        @inbounds for i = 1:n_particles
            for j = (i+1):n_particles
                pair_i[idx] = i
                pair_j[idx] = j
                idx += 1
            end
        end

        new(
            n_particles,
            dims,
            dof,
            n_pairs,
            pair_i,
            pair_j,
            Vector{Float64}(undef, n_pairs),
            Matrix{Bool}(undef, n_particles, n_particles),
            Vector{Bool}(undef, n_particles),
            Vector{Int}(undef, n_particles),
            Vector{Int}(undef, n_particles),
            Vector{Bool}(undef, n_particles),
            Vector{Int}(undef, n_particles),
            0,
            false,
            REG_MODE_NONE,
            0,
            0,
            r_on,
            r_off,
            effective_backend,
            backend_fallback,
            Vector{Int}(undef, n_particles),
            Vector{Bool}(undef, n_particles),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, 2),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, max(dims, 3)),
            Vector{Float64}(undef, 4),
            Vector{Float64}(undef, 4),
            Matrix{Float64}(undef, 3, 4),
            Vector{Float64}(undef, 4),
            Vector{Float64}(undef, 2n_particles + 1 + n_pairs),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
            Vector{Float64}(undef, dof),
        )
    end
end

"""
    WeberProblem(system, tspan, q_initial, p_initial; kwargs...)

Fully specified n-body Weber electrodynamics problem ready for integration.

Packages the compiled `WeberSystem`, initial conditions, physical parameters,
and solver/regularization options into a single immutable structure.

# Arguments
- `system::WeberSystem`: Pre-built symbolic + compiled Hamiltonian system.
- `tspan::Tuple{Real,Real}`: Integration interval `(t_start, t_end)`.
- `q_initial::AbstractVector`: Flattened initial positions, length = `n_particles × dims`.
- `p_initial::AbstractVector`: Flattened initial momenta, length = `n_particles × dims`.

# Keywords
- `masses`: Particle masses (all positive), length `n_particles`.
- `charges`: Particle charges, length `n_particles`.
- `c`: Speed of light (positive).
- `dt`: Fixed macro time step (positive).
- `convergence_tolerance=1e-13`: Fixed-point convergence threshold for projection.
- `maximum_iterations=100`: Maximum projection iterations per step.
- `regularization_enabled=false`, `regularization_r_on`, `regularization_r_off`,
  `regularization_r_on_factor=0.15`, `regularization_r_off_factor=0.25`,
  `regularization_max_substeps=512`, `regularization_constraint_tolerance=1e-12`,
  `regularization_g_floor=1e-12`, `regularization_chain_enabled=true`,
  `regularization_backend=:lifted_pair`, `regularization_warn_on_fallback=true`,
  `regularization_collision_bounce_radius=0.0`:
  All forwarded to `RegularizationOptions`; see its documentation.
- `zollner_enabled=false`, `zollner_a=0.0`: Forwarded to `ZollnerOptions`.

# Fields
- `system::WeberSystem`: Compiled Hamiltonian system.
- `tspan::Tuple{Float64,Float64}`: Integration interval.
- `q_initial`, `p_initial::Vector{Float64}`: Initial phase-space point.
- `masses`, `charges::Vector{Float64}`: Physical parameters.
- `c::Float64`: Speed of light.
- `kappas::Vector{Float64}`: Per-pair Zöllner coupling factors κ_ij.
- `params::Vector{Float64}`: Packed parameter vector `[masses; charges; c; kappas]`.
- `dt::Float64`: Fixed step size.
- `convergence_tolerance::Float64`: Projection convergence threshold.
- `maximum_iterations::Int`: Maximum projection iterations per step.
- `regularization::RegularizationOptions`: Regularization configuration.
- `zollner::ZollnerOptions`: Zöllner extension configuration.
"""
struct WeberProblem
    system::WeberSystem
    tspan::Tuple{Float64,Float64}
    q_initial::Vector{Float64}
    p_initial::Vector{Float64}
    masses::Vector{Float64}
    charges::Vector{Float64}
    c::Float64
    kappas::Vector{Float64}
    params::Vector{Float64}
    dt::Float64
    convergence_tolerance::Float64
    maximum_iterations::Int
    regularization::RegularizationOptions
    zollner::ZollnerOptions

    function WeberProblem(
        system::WeberSystem,
        tspan::Tuple{Real,Real},
        q_initial::AbstractVector,
        p_initial::AbstractVector;
        masses::AbstractVector{<:Real},
        charges::AbstractVector{<:Real},
        c::Real,
        dt::Real,
        convergence_tolerance::Real = 1e-13,
        maximum_iterations::Integer = 100,
        regularization_enabled::Bool = false,
        regularization_r_on::Union{Nothing,Real} = nothing,
        regularization_r_off::Union{Nothing,Real} = nothing,
        regularization_r_on_factor::Real = 0.15,
        regularization_r_off_factor::Real = 0.25,
        regularization_max_substeps::Integer = 512,
        regularization_constraint_tolerance::Real = 1e-12,
        regularization_g_floor::Real = 1e-12,
        regularization_chain_enabled::Bool = true,
        regularization_backend::Symbol = REG_BACKEND_LIFTED,
        regularization_warn_on_fallback::Bool = true,
        regularization_collision_bounce_radius::Real = 0.0,
        zollner_enabled::Bool = false,
        zollner_a::Real = 0.0,
    )
        n_particles = system.n_particles
        dof = system.degrees_of_freedom
        @assert length(q_initial) == dof "q_initial must have length $dof (got $(length(q_initial)))"
        @assert length(p_initial) == dof "p_initial must have length $dof (got $(length(p_initial)))"
        @assert length(masses) == n_particles "masses must have length $n_particles (got $(length(masses)))"
        @assert length(charges) == n_particles "charges must have length $n_particles (got $(length(charges)))"
        @assert all(m -> m > 0, masses) "All masses must be positive"
        @assert c > 0 "Speed of light must be positive"
        @assert tspan[2] > tspan[1] "End time must be greater than start time"
        @assert dt > 0 "Time step must be positive"
        @assert convergence_tolerance > 0 "Convergence tolerance must be positive"
        @assert maximum_iterations > 0 "Maximum iterations must be positive"

        masses_f64 = Vector{Float64}(masses)
        charges_f64 = Vector{Float64}(charges)
        c_f64 = Float64(c)

        zollner = ZollnerOptions(enabled = zollner_enabled, a = Float64(zollner_a))
        kappas = _compute_zollner_kappas(charges_f64, zollner, n_particles)
        params = vcat(masses_f64, charges_f64, [c_f64], kappas)

        regularization = RegularizationOptions(
            enabled = regularization_enabled,
            r_on = regularization_r_on,
            r_off = regularization_r_off,
            r_on_factor = regularization_r_on_factor,
            r_off_factor = regularization_r_off_factor,
            max_substeps = regularization_max_substeps,
            constraint_tolerance = regularization_constraint_tolerance,
            g_floor = regularization_g_floor,
            chain_enabled = regularization_chain_enabled,
            backend = regularization_backend,
            warn_on_fallback = regularization_warn_on_fallback,
            collision_bounce_radius = regularization_collision_bounce_radius,
        )

        new(
            system,
            (Float64(tspan[1]), Float64(tspan[2])),
            Vector{Float64}(q_initial),
            Vector{Float64}(p_initial),
            masses_f64,
            charges_f64,
            c_f64,
            kappas,
            params,
            Float64(dt),
            Float64(convergence_tolerance),
            Int(maximum_iterations),
            regularization,
            zollner,
        )
    end
end

"""
    WeberSolution

Result returned by `solve` or `solve!`.

Supports Julia iteration (`for (t, q, p) in sol`), integer indexing (`sol[i]`),
and `length(sol)`. Each index returns a `(t, q, p)` tuple.

# Fields
- `t::Vector{Float64}`: Time points.
- `q::Vector{Vector{Float64}}`: Flattened position snapshots, one per time point.
- `p::Vector{Vector{Float64}}`: Flattened momentum snapshots, one per time point.
- `prob::WeberProblem`: The originating problem definition.
- `retcode::Symbol`: `:Success` on normal completion, `:Failure` if the
  projection fixed-point failed to converge.
- `regularization::RegularizationDiagnostics`: Regularization usage statistics.
"""
struct WeberSolution
    t::Vector{Float64}
    q::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    prob::WeberProblem
    retcode::Symbol
    regularization::RegularizationDiagnostics
end

Base.length(sol::WeberSolution) = length(sol.t)
Base.getindex(sol::WeberSolution, i::Int) = (sol.t[i], sol.q[i], sol.p[i])
Base.firstindex(sol::WeberSolution) = 1
Base.lastindex(sol::WeberSolution) = length(sol)

function Base.iterate(sol::WeberSolution, state = 1)
    state > length(sol) && return nothing
    return (sol[state], state + 1)
end

function Base.show(io::IO, sol::WeberSolution)
    print(io, "WeberSolution with $(length(sol)) timesteps (retcode: $(sol.retcode))")
end

function Base.show(io::IO, ::MIME"text/plain", sol::WeberSolution)
    println(io, "WeberSolution")
    println(io, "  retcode: $(sol.retcode)")
    println(io, "  t: $(sol.t[1]) → $(sol.t[end]) ($(length(sol)) points)")
    println(io, "  DOF: $(length(sol.q[1]))")
    if sol.regularization.enabled
        println(io, "  Regularization backend: requested=$(sol.regularization.requested_backend), used=$(sol.regularization.used_backend)")
        println(io, "  Regularization steps: pair=$(sol.regularization.pair_steps), chain=$(sol.regularization.chain_steps), cartesian=$(sol.regularization.unregularized_steps)")
    end
end

@inline function _compute_min_pair_distance(
    q::Vector{Float64},
    n_particles::Int,
    dims::Int,
)::Float64
    min_r = Inf
    @inbounds for i = 1:n_particles
        qi_start = (i - 1) * dims
        for j = (i + 1):n_particles
            qj_start = (j - 1) * dims
            r2 = 0.0
            for d = 1:dims
                dq = q[qi_start+d] - q[qj_start+d]
                r2 += dq * dq
            end
            r = sqrt(r2)
            if r < min_r
                min_r = r
            end
        end
    end
    return min_r
end

@inline function _resolve_regularization_backend(
    dims::Int,
    options::RegularizationOptions,
)
    requested = options.backend
    if requested == REG_BACKEND_LIFTED && dims != 2
        return REG_BACKEND_ADAPTIVE, true
    end
    return requested, false
end

mutable struct SymmetricProjectionBuffers
    d::Int
    A::Matrix{Float64}
    A_transpose::Transpose{Float64,Matrix{Float64}}
    Z::Vector{Float64}
    Ẑ::Vector{Float64}
    Z_result::Vector{Float64}
    position_buffer::Vector{Float64}
    auxiliary_position_buffer::Vector{Float64}
    momentum_buffer::Vector{Float64}
    auxiliary_momentum_buffer::Vector{Float64}
    ATμ::Vector{Float64}
    μ::Vector{Float64}
    μ_prev::Vector{Float64}
    f_μ::Vector{Float64}
    diff_buffer::Vector{Float64}
    regularization_buffers::RegularizationBuffers

    function SymmetricProjectionBuffers(prob::WeberProblem)
        d = prob.system.degrees_of_freedom
        n_particles = prob.system.n_particles
        dims = prob.system.dims

        A = zeros(Float64, 2d, 4d)
        @inbounds for i = 1:d
            A[i, i] = 1.0
            A[i, d+i] = -1.0
            A[d+i, 2d+i] = 1.0
            A[d+i, 3d+i] = -1.0
        end
        A_transpose = transpose(A)

        min_pair_distance = _compute_min_pair_distance(prob.q_initial, n_particles, dims)
        if !isfinite(min_pair_distance) || min_pair_distance <= 0
            min_pair_distance = 1.0
        end

        reg_options = prob.regularization
        r_on = isnothing(reg_options.r_on) ? (reg_options.r_on_factor * min_pair_distance) : reg_options.r_on
        r_off = isnothing(reg_options.r_off) ? (reg_options.r_off_factor * min_pair_distance) : reg_options.r_off
        r_on_value = Float64(max(r_on, reg_options.g_floor))
        r_off_value = Float64(max(r_off, r_on_value * 1.01))
        effective_backend, backend_fallback = _resolve_regularization_backend(dims, reg_options)
        regularization_buffers = RegularizationBuffers(
            n_particles,
            dims,
            d,
            r_on_value,
            r_off_value,
            effective_backend,
            backend_fallback,
        )

        new(
            d,
            A,
            A_transpose,
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, 4d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, d),
            Vector{Float64}(undef, 4d),
            zeros(Float64, 2d),
            Vector{Float64}(undef, 2d),
            Vector{Float64}(undef, 2d),
            Vector{Float64}(undef, 2d),
            regularization_buffers,
        )
    end
end

"""
    WeberIntegrator

Mutable step-by-step integrator returned by `init`.

Use `step!(integrator)` to advance one macro-step, or `solve!(integrator)` to
run to completion. The current state is accessible via `integrator.q`,
`integrator.p`, and `integrator.t`.

# Fields
- `prob::WeberProblem`: Problem definition.
- `alg::SymmetricProjectionIntegrator`: Algorithm parameters.
- `t::Float64`: Current time.
- `t_end::Float64`: Final time (`prob.tspan[2]`).
- `q::Vector{Float64}`: Current flattened positions.
- `p::Vector{Float64}`: Current flattened momenta.
- `step_count::Int`: Number of macro-steps completed so far.
- `buffers::SymmetricProjectionBuffers`: Pre-allocated workspace (internal).
- `diagnostics::RegularizationDiagnostics`: Live regularization statistics.
- `t_history::Vector{Float64}`: Pre-allocated time history array.
- `q_history::Vector{Vector{Float64}}`: Pre-allocated position history.
- `p_history::Vector{Vector{Float64}}`: Pre-allocated momentum history.
"""
mutable struct WeberIntegrator
    prob::WeberProblem
    alg::SymmetricProjectionIntegrator
    t::Float64
    t_end::Float64
    q::Vector{Float64}
    p::Vector{Float64}
    step_count::Int
    buffers::SymmetricProjectionBuffers
    diagnostics::RegularizationDiagnostics
    t_history::Vector{Float64}
    q_history::Vector{Vector{Float64}}
    p_history::Vector{Vector{Float64}}
end

function Base.show(io::IO, int::WeberIntegrator)
    print(io, "WeberIntegrator at t=$(int.t) (step $(int.step_count))")
end
