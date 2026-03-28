# Regularization

Regularization is an **optional, advanced** feature. The core integrator runs
unregularized by default; pass `regularization_enabled = true` to opt in.

## When to use it

The Weber potential contains a Coulomb singularity at zero separation. For most
simulations with well-separated particles this is not an issue. Enable
regularization when:

- Two or more particles are expected to have close encounters (separation
  approaching zero).
- You observe energy drift or `NaN` values near close passages without it.

Regularization handles the Coulomb/Kepler singularity only. Weber's
velocity-dependent correction is **not** regularized — the LC and KS transforms
apply to the conservative part of the force.

## Enabling regularization

```julia
prob = WeberProblem(
    sys, tspan, q0, p0;
    masses  = [1.0, 1.0],
    charges = [1.0, -1.0],
    c       = 10.0,
    dt      = 0.01,
    regularization_enabled = true,   # opt in
)
```

## Backends

| Backend | Symbol | Dimensions | Method |
|---------|--------|-----------|--------|
| Levi-Civita (default) | `:lifted_pair` | **2D only** | Lifts the pair to ℝ⁴ fictitious time |
| Adaptive Cartesian | `:adaptive_cartesian` | 2D and 3D | Cartesian sub-stepping with KS constraint |

For 3D problems, `:lifted_pair` automatically falls back to `:adaptive_cartesian`
(with an optional warning controlled by `regularization_warn_on_fallback`).

```julia
# Explicit 3D choice
prob = WeberProblem(sys, tspan, q0, p0; ...
    regularization_enabled  = true,
    regularization_backend  = :adaptive_cartesian,
)
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
prob = WeberProblem(sys, tspan, q0, p0; ...
    regularization_enabled = true,
    regularization_r_on    = 0.05,
    regularization_r_off   = 0.10,
)
```

## Chain mode

When three or more particles form a connected close-encounter cluster,
regularization falls back to chain mode (adaptive Cartesian substeps for the
whole cluster). Chain mode is enabled by default; disable with
`regularization_chain_enabled = false`.

## Collision bounce

For head-on (ℓ = 0) collisions between like-charge pairs, a reflection
boundary can be applied before each macro-step. This avoids the non-regularizable
ℓ ≠ 0 singularity (where particles reach r = 0 at infinite speed).

```julia
prob = WeberProblem(sys, tspan, q0, p0; ...
    regularization_collision_bounce_radius = 0.02,  # reflect at r < 0.02
)
```

Collision bounce works best **without** Levi-Civita regularization (the
unregularized symplectic integrator keeps energy bounded across the bounce).
See [CollisionBounceRegularization](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/research/exploratory/CollisionBounceRegularization.md) for details.

## Diagnostics

Every `WeberSolution` carries a `regularization::RegularizationDiagnostics`
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
| `backend_fallback_steps` | Steps that used the fallback backend | > 0 when `:lifted_pair` was requested on a 3D problem |
| `total_substeps` | Total regularization micro-steps taken | Very large → consider wider r_on/r_off thresholds |

## API reference

```@docs
RegularizationOptions
RegularizationDiagnostics
```
