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

The wrapper's options override `prob.regularization`; the base algorithm's
settings (e.g. `relaxation`) are preserved.
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

# Rewrite a HamiltonianProblem with a new `regularization` field. This is the
# Phase 3c.1 shim that lets `RegularizedIntegrator` delegate to the existing
# regularization dispatch code without refactoring hot-path option reads.
# Phase 3d will drop `prob.regularization` entirely, at which point this
# helper disappears and the dispatch reads options from the algorithm instead.
function _with_regularization(prob::HamiltonianProblem, reg::RegularizationOptions)
    return HamiltonianProblem(
        prob.system,
        prob.tspan,
        prob.q_initial,
        prob.p_initial;
        masses = prob.masses,
        charges = prob.charges,
        c = prob.c,
        dt = prob.dt,
        convergence_tolerance = prob.convergence_tolerance,
        maximum_iterations = prob.maximum_iterations,
        regularization = reg,
        zollner = prob.zollner,
    )
end

function CommonSolve.init(
    prob::HamiltonianProblem,
    alg::RegularizedIntegrator;
    callbacks = (),
)
    effective_prob = _with_regularization(prob, alg.options)
    return CommonSolve.init(effective_prob, alg.base_alg; callbacks = callbacks)
end
