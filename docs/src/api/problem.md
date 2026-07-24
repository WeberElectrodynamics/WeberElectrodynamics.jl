# Problem

## Problem definition

```@docs
HamiltonianProblem
```

## Accessors

Read-only accessors on a `HamiltonianProblem`. The backing `params` vector
has layout `[m₁…mₙ, q₁…qₙ, c]` (length `2N + 1`). `masses(prob)` and
`charges(prob)` are O(1) views into `params(prob)`; `speed_of_light(prob)`
reads the trailing `c` element.

| Function | Returns |
|----------|---------|
| `masses(prob)` | `AbstractVector{Float64}` view of particle masses (length `N`) |
| `charges(prob)` | `AbstractVector{Float64}` view of particle charges (length `N`) |
| `speed_of_light(prob)` | `Float64` — the `c` passed at construction |
| `params(prob)` | the full packed `Vector{Float64}` consumed by the compiled EOMs |

See also `n_particles(prob)` and `dims(prob)`, documented on the [System](system.md) page.

```@docs
masses
charges
speed_of_light
params
```

## Initial-condition helpers

These helpers return named tuples with flattened `q` and `p` vectors ready for
`HamiltonianProblem`.

```@docs
center_of_mass_frame
two_body_initial_conditions
polygon_initial_conditions
rigid_rotation_initial_conditions
```

## Optional features

Configuration types for optional features are documented on their own pages:

- [Regularization](../regularization.md) — `RegularizationOptions`, `RegularizationDiagnostics`
