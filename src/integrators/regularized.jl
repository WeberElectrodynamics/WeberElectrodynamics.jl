"""
    RegularizedIntegrator{A}(base_alg::A; kwargs...)

Algorithm wrapper that activates Levi-Civita / KS regularization on top of a
base Hamiltonian algorithm. `kwargs` mirror [`RegularizationOptions`](@ref)
and are forwarded verbatim; `enabled` is set to `true` implicitly.

This is the user-facing entry point for regularized integration:

```julia
alg = RegularizedIntegrator(SymmetricProjectionIntegrator();
                            r_on_factor = 0.15, backend = :lifted_pair)
sol = solve(prob, alg)
```

The base algorithm's settings (e.g. `relaxation`) are preserved.
"""
struct RegularizedIntegrator{A<:HamiltonianAlgorithm} <: HamiltonianAlgorithm
    base_alg::A
    options::RegularizationOptions
end

function RegularizedIntegrator(
    base_alg::HamiltonianAlgorithm;
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
    options = RegularizationOptions(
        enabled = true,
        r_on = r_on,
        r_off = r_off,
        r_on_factor = r_on_factor,
        r_off_factor = r_off_factor,
        max_substeps = max_substeps,
        constraint_tolerance = constraint_tolerance,
        g_floor = g_floor,
        chain_enabled = chain_enabled,
        backend = backend,
        warn_on_fallback = warn_on_fallback,
        collision_bounce_radius = collision_bounce_radius,
    )
    return RegularizedIntegrator(base_alg, options)
end

"""
    base_algorithm(alg::HamiltonianAlgorithm) -> HamiltonianAlgorithm

Unwrap an algorithm one level. For bare algorithms this is the identity; for
[`RegularizedIntegrator`](@ref) it returns the wrapped base algorithm.
"""
base_algorithm(alg::HamiltonianAlgorithm) = alg
base_algorithm(alg::RegularizedIntegrator) = alg.base_alg

@inline function _step_core!(
    integrator::HamiltonianIntegrator,
    ::RegularizedIntegrator,
    dt_step::Float64,
)
    _step_regularized_dispatch!(integrator, dt_step)
    return nothing
end

function _allocate_cache(prob::HamiltonianProblem, alg::RegularizedIntegrator)
    buffers = SymmetricProjectionBuffers(prob, alg.options)
    rb = buffers.regularization_buffers
    if rb.backend_fallback && alg.options.warn_on_fallback
        @warn "RegularizationOptions(backend=:lifted_pair) is currently supported only for 2D; falling back to :adaptive_cartesian for $(prob.system.dims)D"
    end
    if !rb.backend_fallback && alg.options.warn_on_fallback && has_term(prob.system, :weber)
        @warn "RegularizedIntegrator handles close encounters by substepping/lifting the Coulomb singularity; Weber's velocity-dependent force is not analytically regularized" maxlog = 1
    end
    return buffers
end

function _resolve_callbacks(::HamiltonianProblem, alg::RegularizedIntegrator, cbs)
    user = _normalise_callbacks(cbs)
    bounce_r = alg.options.collision_bounce_radius
    if bounce_r > 0 && !any(c -> c isa CollisionBounce, user)
        return (CollisionBounce(bounce_r), user...)
    end
    return user
end
