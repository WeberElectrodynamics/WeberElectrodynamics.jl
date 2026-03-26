# System

The `WeberSystem` encapsulates the symbolic and compiled Weber Hamiltonian for
a given `(n_particles, dims)` configuration. Construct it once and reuse it
across many `WeberProblem` instances.

```@docs
WeberSystem
WeberSystem(::Int, ::Int)
```
