# System

The `HamiltonianSystem` encapsulates the compiled Hamiltonian for a given
`(n_particles, dims)` configuration. Construct it once and reuse it across many
`HamiltonianProblem` instances.

There are two construction paths — an **analytic** one used by the built-in
Weber system, and a **symbolic** one for custom Hamiltonians. See
[Custom Hamiltonians](@ref) for the full comparison, and
[The Weber Hamiltonian](@ref) for why the Weber system is analytic.

```@docs
HamiltonianSystem
HamiltonianSystem(::Int, ::Int)
HamiltonianSystem(::Any, ::AbstractVector, ::AbstractVector)
has_symbolic_hamiltonian
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

## Physical velocities

```@docs
physical_velocities
WeberCriticalRadiusError
```

Canonical momentum in Weber electrodynamics is `p_i = ∂L/∂v_i`, which is not
`m_i v_i`. Always convert with `physical_velocities`, never `p ./ masses(prob)`.

## Term builders

Composable building blocks for symbolic Hamiltonians.

```@docs
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

A `NamedTerm` may attach optional `pair_decomposition(q, p, params)` and
`kinetic_energy(q, p, params)` closures — used by the energy/force statistics
to avoid re-deriving per-term decompositions. Both take the whole state rather
than one pair, because recovering physical velocities for a velocity-dependent
Hamiltonian is a single coupled solve over all pairs. `pair_decomposition`
returns per-pair vectors ordered as `pair_indices`. The default `:weber` term
supplies both; custom terms may attach any `NamedTuple` shape, but the built-in
statistics consume specific field names (`coulomb`, `velocity`, `rdot`, `r`).

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
