# System

The `HamiltonianSystem` encapsulates the symbolic and compiled Weber Hamiltonian for
a given `(n_particles, dims)` configuration. Construct it once and reuse it
across many `HamiltonianProblem` instances.

```@docs
HamiltonianSystem
HamiltonianSystem(::Int, ::Int)
```

## Term builders

```@docs
weber_term
zollner_term
```

## Term introspection

```@docs
NamedTerm
term_names
has_term
get_term
```

## Metadata accessors

Shape accessors are defined for both `HamiltonianSystem` and `HamiltonianProblem`.

```@docs
n_particles
```

| Function | Returns |
|----------|---------|
| `dims(sys)` / `dims(prob)` | spatial dimension (`1`, `2`, or `3`) |
| `degrees_of_freedom(sys)` | total canonical DOF (`n_particles * dims`) |
