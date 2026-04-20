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

```@docs
n_particles
```
