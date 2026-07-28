"""
    NamedTerm(name, H_symbolic; pair_decomposition = nothing, kinetic_energy = nothing)

Named component of a composite Hamiltonian: a symbol that identifies the
physical role (e.g. `:weber`) and the expression contributing to the total `H`.
Stored on `HamiltonianSystem` as `terms::Vector{NamedTerm}` so statistics and
plotting can query individual contributions without re-deriving them from the
compiled aggregate.

`H_symbolic` is `nothing` for analytic systems such as the default Weber
system, whose Hamiltonian has no practical closed symbolic form.

The optional `pair_decomposition` closure is evaluated once per state as
`pair_decomposition(q, p, params)` and returns a `NamedTuple` of **per-pair
vectors** ordered as [`pair_indices`](@ref). Taking the whole state at once
matters for velocity-dependent Hamiltonians: recovering the physical velocities
of a Weber system is a single coupled solve over all pairs, so a per-pair
signature would either repeat that solve or invite a wrong per-pair shortcut.

The optional `kinetic_energy` closure is evaluated as
`kinetic_energy(q, p, params)` and returns the term's **physical** kinetic
energy `Σ ½ m_i |v_i|²`. It exists because canonical momentum is not `m_i v_i`
for a velocity-dependent Hamiltonian, so `Σ |p_i|²/(2m_i)` is not the kinetic
energy.

# Fields
- `name::Symbol`: Identifier used by `get_term`, `has_term`, `term_names`.
- `H_symbolic`: Symbolic expression contributing to the aggregate Hamiltonian,
  or `nothing` for analytic terms.
- `pair_decomposition::Union{Function,Nothing}`: Optional closure
  `(q, p, params) -> NamedTuple` of per-pair vectors.
- `kinetic_energy::Union{Function,Nothing}`: Optional closure
  `(q, p, params) -> Float64` giving physical kinetic energy.
"""
struct NamedTerm{H,F,K}
    name::Symbol
    H_symbolic::H
    pair_decomposition::F
    kinetic_energy::K
end

NamedTerm(
    name::Symbol,
    H_symbolic;
    pair_decomposition = nothing,
    kinetic_energy = nothing,
) = NamedTerm(name, H_symbolic, pair_decomposition, kinetic_energy)

"""
    term_names(sys::HamiltonianSystem) -> Vector{Symbol}

Return the names of every `NamedTerm` on `sys`, in insertion order.
"""
term_names(terms::AbstractVector{<:NamedTerm}) = Symbol[t.name for t in terms]

"""
    has_term(sys::HamiltonianSystem, name::Symbol) -> Bool

Return `true` if `sys` carries a `NamedTerm` called `name`.
"""
has_term(terms::AbstractVector{<:NamedTerm}, name::Symbol) =
    any(t -> t.name === name, terms)

"""
    get_term(sys::HamiltonianSystem, name::Symbol) -> NamedTerm

Return the `NamedTerm` called `name`. Throws `KeyError` if absent.
"""
function get_term(terms::AbstractVector{<:NamedTerm}, name::Symbol)
    idx = findfirst(t -> t.name === name, terms)
    idx === nothing && throw(KeyError(name))
    return terms[idx]
end
