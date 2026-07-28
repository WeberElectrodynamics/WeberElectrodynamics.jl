"""
    kinetic_term(p_vars, masses; n_particles, dims) -> Num

Build the symbolic non-relativistic kinetic energy `Σᵢ ‖pᵢ‖² / (2 mᵢ)` from
pre-constructed symbolic inputs.

Returns a `Symbolics.Num` expression usable with `Symbolics.derivative` and
`Symbolics.build_function`. Intended for users composing custom Hamiltonians
such as `H = kinetic_term(...) + coulomb_term(...)` without the Weber velocity
term.

# Arguments
- `p_vars`: Momentum symbolic variables, length `n_particles*dims`.

# Keywords
- `masses`: Per-particle mass symbolic variables (length `n_particles`).
- `n_particles::Int`, `dims::Int`: Problem shape.
"""
function kinetic_term(
    p_vars::AbstractVector;
    masses::AbstractVector,
    n_particles::Int,
    dims::Int,
)
    H = zero(eltype(p_vars))

    @inbounds for i = 1:n_particles
        p_start = (i - 1) * dims + 1
        p_end = i * dims
        p_squared = sum(p_vars[p_start:p_end] .^ 2)
        H = H + p_squared / (2 * masses[i])
    end

    return H
end

"""
    coulomb_term(q_vars; charges, n_particles, dims) -> Num

Build the symbolic pairwise Coulomb potential `Σᵢ<ⱼ qᵢ qⱼ / rᵢⱼ` from
pre-constructed symbolic inputs.

Returns a `Symbolics.Num` expression usable with `Symbolics.derivative` and
`Symbolics.build_function`. Unlike the built-in Weber system, this carries no
velocity-dependent correction (no `c`) — it is the limit of pure 1/r Coulomb.

# Arguments
- `q_vars`: Coordinate symbolic variables, length `n_particles*dims`.

# Keywords
- `charges`: Per-particle charge symbolic variables (length `n_particles`).
- `n_particles::Int`, `dims::Int`: Problem shape.
"""
function coulomb_term(
    q_vars::AbstractVector;
    charges::AbstractVector,
    n_particles::Int,
    dims::Int,
)
    H = zero(eltype(q_vars))

    @inbounds for i = 1:n_particles
        for j = (i+1):n_particles
            qi_start = (i - 1) * dims + 1
            qj_start = (j - 1) * dims + 1

            r_squared = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                r_squared = r_squared + dq^2
            end
            r = sqrt(r_squared)

            H = H + charges[i] * charges[j] / r
        end
    end

    return H
end
