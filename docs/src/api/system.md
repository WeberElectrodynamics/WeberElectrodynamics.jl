# System

The `HamiltonianSystem` encapsulates the symbolic and compiled Weber Hamiltonian for
a given `(n_particles, dims)` configuration. Construct it once and reuse it
across many `HamiltonianProblem` instances.

```@docs
HamiltonianSystem
HamiltonianSystem(::Int, ::Int)
HamiltonianSystem(::Any, ::AbstractVector, ::AbstractVector)
```

The generic constructor takes a pre-built symbolic Hamiltonian `H` and the
phase-space variables, then derives Hamilton's equations and compiles them.
Use it to assemble a custom Hamiltonian from term builders, e.g.

```julia
H = kinetic_term(p_vars; masses = m_vars, n_particles, dims) +
    coulomb_term(q_vars; charges = q_charges, n_particles, dims)

sys = HamiltonianSystem(H, q_vars, p_vars;
    param_symbols = vcat(m_vars, q_charges, [c_var]),
    t = t_var, n_particles = n_particles, dims = dims,
)
```

## Term builders

Composable building blocks for symbolic Hamiltonians.

```@docs
weber_term
kinetic_term
coulomb_term
```

## Term introspection

```@docs
NamedTerm
term_names
has_term
get_term
```

A `NamedTerm` may attach an optional `pair_decomposition(i, j, q, p, params)`
closure that returns a `NamedTuple` summarising the term's contribution to a
single pair — used by the energy/force statistics to avoid re-deriving
per-term decompositions. The default `:weber` term supplies its own closure;
custom terms may attach any `NamedTuple` shape, but
the built-in statistics consume specific field names (`coulomb`, `velocity`,
`rdot`, `r`).

## Metadata accessors

Shape accessors are defined for both `HamiltonianSystem` and `HamiltonianProblem`.

```@docs
n_particles
n_pairs
pair_indices
```

| Function | Returns |
|----------|---------|
| `dims(sys)` / `dims(prob)` | spatial dimension (`1`, `2`, or `3`) |
| `degrees_of_freedom(sys)` | total canonical DOF (`n_particles * dims`) |
