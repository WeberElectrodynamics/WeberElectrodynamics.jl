# Problem

## Problem definition

```@docs
HamiltonianProblem
```

## Accessors

Read-only accessors on a `HamiltonianProblem`. `masses`, `charges`, and
`kappas` return views into the backing `params` vector (layout
`[m₁…mₙ, q₁…qₙ, c, κ…]`), so they are O(1) and allocation-free.

| Function | Returns |
|----------|---------|
| `masses(prob)` | `AbstractVector{Float64}` of particle masses |
| `charges(prob)` | `AbstractVector{Float64}` of particle charges |
| `speed_of_light(prob)` | `Float64` — the `c` passed at construction |
| `kappas(prob)` | `AbstractVector{Float64}` — per-pair Zöllner κ values (all `1.0` when Zöllner is disabled) |
| `params(prob)` | the full packed `Vector{Float64}` consumed by the compiled EOMs |

See also `n_particles(prob)` and `dims(prob)`, documented on the [System](system.md) page.

## Optional features

Configuration types for optional features are documented on their own pages:

- [Regularization](../regularization.md) — `RegularizationOptions`, `RegularizationDiagnostics`
- [Zöllner Extension](../zollner.md) — `ZollnerOptions`
