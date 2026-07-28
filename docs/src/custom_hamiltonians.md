# Custom Hamiltonians

The package ships one Hamiltonian baked in — the exact canonical Weber system
built by `HamiltonianSystem(n_particles, dims)` — and one extension point for
everything else.

## The two construction paths

| | Analytic | Symbolic |
|---|---|---|
| Constructor | `HamiltonianSystem(n, dims)` | `HamiltonianSystem(H, q_vars, p_vars; …)` |
| Equations of motion | hand-derived closed forms | `Symbolics.derivative` + `build_function` |
| `hamiltonian_symbolic` | `nothing` | the expression you passed |
| [`has_symbolic_hamiltonian`](@ref) | `false` | `true` |
| Used for | the built-in Weber system | your Hamiltonian |

Both fill the same compiled slots, so the solver, statistics, plotting, and
animation never need to know which path produced a system:

```
sys.dq_dt_compiled(out, q, p, t, params)
sys.dp_dt_compiled(out, q, p, t, params)
sys.hamiltonian_compiled(q, p, t, params)
```

The default Weber system takes the analytic path because recovering physical
velocities from canonical momenta requires a coupled linear solve over particle
pairs at every evaluation — see [The Weber Hamiltonian](@ref). That solve has no
practical closed symbolic form for general `n`, and symbolic elimination of the
canonical mass matrix fails outright.

## Building a symbolic Hamiltonian

Compose from the exported term builders, or write the expression directly:

```julia
using WeberElectrodynamics, Symbolics

@variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 q1 q2 c t

q_vars = [x1, y1, x2, y2]
p_vars = [px1, py1, px2, py2]

H = kinetic_term(p_vars; masses = [m1, m2], n_particles = 2, dims = 2) +
    coulomb_term(q_vars; charges = [q1, q2], n_particles = 2, dims = 2)

sys = HamiltonianSystem(H, q_vars, p_vars;
    param_symbols = [m1, m2, q1, q2, c],
    t = t, n_particles = 2, dims = 2,
)
```

`param_symbols` must follow the package layout `[m₁…mₙ, q₁…qₙ, c]` of length
`2n+1`, because [`HamiltonianProblem`](@ref) packs `masses`, `charges`, and `c`
into that order. Include `c` even when the Hamiltonian does not use it.

The resulting system drops straight into the normal pipeline:

```julia
prob = HamiltonianProblem(sys, (0.0, 10.0), q0, p0;
    masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 1e-3)
sol = solve(prob, SymmetricProjectionIntegrator())
```

### Term builders

[`kinetic_term`](@ref) and [`coulomb_term`](@ref) are documented on the
[System](@ref) API page.

## Naming components with `NamedTerm`

Pass a `terms` vector to keep the decomposition of `H` queryable:

```julia
kinetic = kinetic_term(p_vars; masses = [m1, m2], n_particles = 2, dims = 2)
coulomb = coulomb_term(q_vars; charges = [q1, q2], n_particles = 2, dims = 2)

sys = HamiltonianSystem(kinetic + coulomb, q_vars, p_vars;
    param_symbols = [m1, m2, q1, q2, c], t = t, n_particles = 2, dims = 2,
    terms = [NamedTerm(:kinetic, kinetic), NamedTerm(:coulomb, coulomb)],
)

term_names(sys)          # [:kinetic, :coulomb]
get_term(sys, :coulomb)  # NamedTerm
```

Without a `terms` vector the system stores a single `NamedTerm(:hamiltonian, H)`.

### Statistics hooks

A `NamedTerm` may carry two optional closures. Both take the whole state, not a
single pair, because velocity-dependent Hamiltonians need one coupled solve per
state rather than one per pair:

- `pair_decomposition(q, p, params)` returns a `NamedTuple` of per-pair vectors
  ordered as [`pair_indices`](@ref).
- `kinetic_energy(q, p, params)` returns the **physical** kinetic energy
  `Σ ½ mᵢ|vᵢ|²`.

`compute_energy_timeseries` uses both when the system carries a `:weber` term
whose `pair_decomposition` supplies the fields `coulomb`, `velocity`, `rdot`,
and `r`. Systems without them still get a valid total-energy timeseries from the
compiled Hamiltonian; `pair_energies` is then empty and `kinetic_energy` falls
back to the canonical split `Σ|pᵢ|²/(2mᵢ)`, which is correct only for
Hamiltonians whose kinetic term is genuinely `Σ|p|²/(2m)`.

```julia
custom_pd = (q, p, params) -> begin
    r = sqrt((q[1] - q[3])^2 + (q[2] - q[4])^2)
    return (coulomb = [params[3] * params[4] / r], velocity = [0.0],
            rdot = [0.0], r = [r])
end

NamedTerm(:my_term, H; pair_decomposition = custom_pd)
```

## Writing an analytic system

If your Hamiltonian also resists symbolic differentiation, construct a
`HamiltonianSystem` directly with `nothing` in the three symbolic fields and
your own closures in the compiled ones — this is exactly what the Weber
constructor does. The field order is

```
n_particles, dims, q_symbols, p_symbols, t_symbol, param_symbols,
hamiltonian_symbolic, dq_dt_symbolic, dp_dt_symbolic,
dq_dt_compiled, dp_dt_compiled, hamiltonian_compiled,
degrees_of_freedom, terms
```

Your `dq_dt_compiled` and `dp_dt_compiled` must be genuine partial derivatives
of your `hamiltonian_compiled`, or the symplectic integrator's conservation
guarantees do not apply. Verify with central differences before trusting
results — `test/test_canonical_weber.jl` shows the pattern used for the built-in
Weber system.

Give each compiled entry point its own workspace if you cache buffers in
closures; the built-in system does, since a shared workspace would not be
reentrant.

## Caveats

- Regularization backends (`:lifted_pair`, `:adaptive_cartesian`) assume a
  Coulomb-like pair singularity. They still run for arbitrary Hamiltonians, but
  the accuracy argument only holds for Kepler-like pair potentials.
- `save_solution`/`load_solution` archive default Weber systems only; a custom
  system cannot be reconstructed from the archive metadata.
- The compiled `t` argument is reserved for time-dependent terms and is
  currently unused by the built-in Weber system.
