"""
    NamedTerm(name, H_symbolic; pair_decomposition = nothing)

Named component of a composite Hamiltonian: a symbol that identifies the
physical role (e.g. `:weber`, `:zollner`) and the symbolic expression
contributing to the total `H`. Stored on `HamiltonianSystem` as
`terms::Vector{NamedTerm}` so statistics and plotting can query individual
contributions without re-deriving them from the compiled aggregate.

The optional `pair_decomposition` closure returns a per-pair decomposition of
the term's contribution, evaluated on concrete `(i, j, q, p, params)` inputs.
The shape of the returned tuple is term-specific; the Weber and Zöllner
builders will attach compatible closures in Phase 4 when statistics migrate
onto this interface.

# Fields
- `name::Symbol`: Identifier used by `get_term`, `has_term`, `term_names`.
- `H_symbolic`: Symbolic expression contributing to the aggregate Hamiltonian.
- `pair_decomposition::Union{Function,Nothing}`: Optional closure
  `(i, j, q, p, params) -> NamedTuple` for per-pair statistics.
"""
struct NamedTerm{H,F}
    name::Symbol
    H_symbolic::H
    pair_decomposition::F
end

NamedTerm(name::Symbol, H_symbolic; pair_decomposition = nothing) =
    NamedTerm(name, H_symbolic, pair_decomposition)

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
