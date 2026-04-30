# Regularization

Regularization is an **optional, advanced** feature. The core integrator runs
unregularized by default; opt in by wrapping the base algorithm in a
[`RegularizedIntegrator`](@ref) and passing it to `solve`.

## When to use it

The Weber potential contains a Coulomb singularity at zero separation. For most
simulations with well-separated particles this is not an issue. Enable
regularization when:

- Two or more particles are expected to have close encounters (separation
  approaching zero).
- You observe energy drift or `NaN` values near close passages without it.

Regularization handles the Coulomb/Kepler singularity only. Weber's
velocity-dependent correction is **not analytically regularized**; it is still
evaluated through the package's equations of motion during regularized
substeps.

## Enabling regularization

```julia
prob = HamiltonianProblem(
    sys, tspan, q0, p0;
    masses  = [1.0, 1.0],
    charges = [1.0, -1.0],
    c       = 10.0,
    dt      = 0.01,
)

alg = RegularizedIntegrator(SymmetricProjectionIntegrator())   # opt in
sol = solve(prob, alg)
```

All keyword arguments below go on the `RegularizedIntegrator(...)` constructor
itself. Its kwargs mirror [`RegularizationOptions`](@ref) (with `enabled = true`
set implicitly).

## Backends

| Backend | Symbol | Dimensions | Method |
|---------|--------|-----------|--------|
| Lifted pair (default) | `:lifted_pair` | 1D | Square-root chart with fictitious time |
| Lifted pair (default) | `:lifted_pair` | 2D | Levi-Civita pair chart with fictitious time |
| Lifted pair (default) | `:lifted_pair` | 3D | KS pair chart with constraint projection diagnostics |
| Adaptive Cartesian | `:adaptive_cartesian` | 1D, 2D, 3D | Cartesian close-encounter substeps |
| Chain mode | automatic | 1D, 2D, 3D | Adaptive Cartesian substeps over multi-particle close clusters |

`RegularizedIntegrator(; kwargs...)` is shorthand for wrapping
`SymmetricProjectionIntegrator()`.

!!! note "3D KS constraint"
    The `:lifted_pair` backend uses a KS pair chart and projects the KS momentum
    constraint at the start, midpoint, and end of each lifted substep. The
    diagnostics report both `max_constraint_violation` and KS projection counts.
    The `:adaptive_cartesian` backend still uses KS lifts for diagnostics only
    before running Cartesian substeps.

```julia
# Explicit Cartesian fallback choice
alg = RegularizedIntegrator(
    SymmetricProjectionIntegrator();
    backend = :adaptive_cartesian,
)
sol = solve(prob, alg)
```

## Activation hysteresis

Regularization switches on when any pair separation drops below `r_on` and
switches off once all separations in the active component exceed `r_off`.
The thresholds are computed from the initial minimum separation by default:

```
r_on  = r_on_factor  × min_initial_separation   (default factor: 0.15)
r_off = r_off_factor × min_initial_separation   (default factor: 0.25)
```

You can override them directly:

```julia
alg = RegularizedIntegrator(
    SymmetricProjectionIntegrator();
    r_on  = 0.05,
    r_off = 0.10,
)
sol = solve(prob, alg)
```

## Chain mode

When three or more particles form a connected close-encounter cluster,
regularization enters chain mode (adaptive Cartesian substeps for the whole
cluster). Chain mode is enabled by default; disable with
`RegularizedIntegrator(...; chain_enabled = false)`.

Chain mode does not implement analytic chain-coordinate regularization; it is a
Cartesian fallback for multi-pair clusters.

## Collision bounce

For head-on (`L = 0`) collisions or pass-through events, a reflection boundary
can be applied before each macro-step. The bounce is geometric rather than
charge-sign-specific: it reflects any pair inside the radius, including a
two-body unlike-charge pair when explicitly enabled. Under
`RegularizedIntegrator`, the same radius is also checked after regularized
substeps, so close approaches inside a regularized macro-step do not have to
wait for the next outer step.

Collision bounce is a [`CollisionBounce`](@ref) callback; pass it through the
`callbacks` kwarg of `solve` (or `init`):

```julia
sol = solve(prob, SymmetricProjectionIntegrator();
            callbacks = CollisionBounce(0.02))   # reflect at r < 0.02
```

As a convenience, `RegularizedIntegrator` also accepts a `collision_bounce_radius`
kwarg and synthesises a matching callback automatically:

```julia
alg = RegularizedIntegrator(SymmetricProjectionIntegrator();
                            collision_bounce_radius = 0.02)
sol = solve(prob, alg)
```

Collision bounce is intended for C0-continuable head-on cases. It preserves the
isolated two-body pair energy exactly, but it is not generally energy-preserving
for `N > 2` because reflected particles move relative to third bodies. It does
not make generic nonzero-angular-momentum collisions regular.

## Diagnostics

Every `HamiltonianSolution` carries a `regularization::RegularizationDiagnostics`
field with step counts, backend used, minimum encounter distance, and more.

```julia
sol = solve(prob, SymmetricProjectionIntegrator())
d = sol.regularization
println("Backend used: ", d.used_backend)
println("Activation count: ", d.activation_count)
println("Min encounter distance: ", d.min_encounter_distance)
```

Key fields and what they mean:

| Field | Meaning | Worth investigating if… |
|-------|---------|------------------------|
| `activation_count` | How many times regularization switched on | Very high → r_off is too close to r_on |
| `min_encounter_distance` | Closest particle approach over the run | Near 0 → singularity risk; consider bounce radius |
| `max_constraint_violation` | Peak KS constraint residual (3D only) | > 1e-8 → reduce `dt` or `max_substeps` |
| `ks_constraint_projection_count` | Number of KS constraint projections | Useful for confirming the 3D lifted path ran |
| `ks_constraint_iteration_count` | Total KS projection iterations | Higher than projections if future iterative correction is enabled |
| `backend_fallback_steps` | Steps that used the fallback backend | > 0 only when a requested backend is unavailable |
| `total_substeps` | Total regularization micro-steps taken | Very large → consider wider r_on/r_off thresholds |

## API reference

```@docs
RegularizationOptions
RegularizationDiagnostics
```
