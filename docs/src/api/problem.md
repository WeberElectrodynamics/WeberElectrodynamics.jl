# Problem

## Problem definition

```@docs
HamiltonianProblem
```

## Accessors

Read-only accessors on a `HamiltonianProblem`. The backing `params` vector
has layout `[m₁…mₙ, q₁…qₙ, c]` (length `2N + 1`). `masses(prob)` and
`charges(prob)` are O(1) views into `params(prob)`; `speed_of_light(prob)`
reads the trailing `c` element. The per-pair Zöllner coupling `kappas`
lives on its own field of length `N(N−1)/2`, separate from `params`.

| Function | Returns |
|----------|---------|
| `masses(prob)` | `AbstractVector{Float64}` view of particle masses (length `N`) |
| `charges(prob)` | `AbstractVector{Float64}` view of particle charges (length `N`) |
| `speed_of_light(prob)` | `Float64` — the `c` passed at construction |
| `kappas(prob)` | `Vector{Float64}` of per-pair Zöllner κ values (all `1.0` when Zöllner is disabled; `Float64[]` for `N = 1`) |
| `kappa(prob, i, j)` | `Float64` — coupling for pair `(i, j)` with `i < j` |
| `params(prob)` | the full packed `Vector{Float64}` consumed by the compiled EOMs |

See also `n_particles(prob)` and `dims(prob)`, documented on the [System](system.md) page.

```@docs
kappa
```

## Optional features

Configuration types for optional features are documented on their own pages:

- [Regularization](../regularization.md) — `RegularizationOptions`, `RegularizationDiagnostics`
- [Zöllner Extension](../zollner.md) — `ZollnerOptions`
