"""
    HamiltonianCallback

Abstract supertype for per-step callbacks applied by the integrator. Concrete
subtypes can override any of:

- [`apply_pre_step!`](@ref)`(cb, integrator, dt_step)` — invoked before
  the algorithm's `_step_core!` advances the state.
- [`apply_post_step!`](@ref)`(cb, integrator, dt_step)` — invoked after
  the state and history have been updated.

Both default to no-ops. Callbacks are type-stably iterated from a tuple stored
on the integrator; there is no per-step allocation for an empty callback set.
"""
abstract type HamiltonianCallback end

"""
    apply_pre_step!(cb::HamiltonianCallback, integrator, dt_step)

Invoked on the integrator's `q`/`p` state immediately before
`_step_core!(integrator, alg, dt_step)` for the upcoming macro step. Default
implementation is a no-op.
"""
apply_pre_step!(::HamiltonianCallback, _, _::Float64) = nothing

"""
    apply_post_step!(cb::HamiltonianCallback, integrator, dt_step)

Invoked after the algorithm's step has completed and history has been
recorded. Default implementation is a no-op.
"""
apply_post_step!(::HamiltonianCallback, _, _::Float64) = nothing

"""
    CollisionBounce(radius)

Pre-step callback that reflects the relative coordinate of any pair closer
than `radius` through the origin, preserving COM and momenta (C⁰-continuation
of the ℓ = 0 head-on collision; Frauenfelder & Weber 2024). Works best with
the unregularized symplectic path. Under `RegularizedIntegrator`, the callback
fires only at macro-step boundaries; close approaches inside regularized
substeps are not reflected until the next outer step.

See `docs/src/regularization.md` and `theory/Regularization.md` for the
regularization and collision-continuation context.
"""
struct CollisionBounce <: HamiltonianCallback
    radius::Float64
    function CollisionBounce(radius::Real)
        @assert radius >= 0 "CollisionBounce radius must be non-negative"
        new(Float64(radius))
    end
end

function apply_pre_step!(cb::CollisionBounce, integrator, ::Float64)
    cb.radius > 0 || return nothing
    prob = integrator.prob
    _apply_collision_bounces!(
        integrator.q,
        masses(prob),
        prob.system.dims,
        prob.system.n_particles,
        cb.radius,
    )
    return nothing
end

# Type-stable fanout over a tuple of callbacks. Manually recursing over the
# tuple keeps @inline happy and avoids allocation when the set is empty.
@inline _run_pre_step!(::Tuple{}, _, _::Float64) = nothing
@inline function _run_pre_step!(cbs::Tuple, integrator, dt_step::Float64)
    apply_pre_step!(first(cbs), integrator, dt_step)
    _run_pre_step!(Base.tail(cbs), integrator, dt_step)
    return nothing
end

@inline _run_post_step!(::Tuple{}, _, _::Float64) = nothing
@inline function _run_post_step!(cbs::Tuple, integrator, dt_step::Float64)
    apply_post_step!(first(cbs), integrator, dt_step)
    _run_post_step!(Base.tail(cbs), integrator, dt_step)
    return nothing
end

# Normalise user-supplied callback inputs into a concrete tuple. Accepts
# nothing, a single callback, or any iterable of callbacks.
_normalise_callbacks(::Nothing) = ()
_normalise_callbacks(cb::HamiltonianCallback) = (cb,)
_normalise_callbacks(cbs::Tuple{Vararg{HamiltonianCallback}}) = cbs
function _normalise_callbacks(cbs)
    callbacks = ()
    for cb in Iterators.reverse(collect(cbs))
        cb isa HamiltonianCallback ||
            throw(ArgumentError("callbacks must all subtype HamiltonianCallback"))
        callbacks = (cb, callbacks...)
    end
    return callbacks
end
