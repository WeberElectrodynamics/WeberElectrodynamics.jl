"""
    weber_term(q_vars, p_vars;
               masses, charges, c, kappas,
               n_particles, dims) -> Num

Build the symbolic Weber Hamiltonian (kinetic term plus κ-weighted pairwise
velocity-dependent potential) from pre-constructed symbolic inputs.

Returns a `Symbolics.Num` expression usable with `Symbolics.derivative` and
`Symbolics.build_function`. This is the symbolic builder at the heart of the
default `HamiltonianSystem(n_particles, dims)` constructor; calling it
directly lets a user compose a custom Hamiltonian, e.g.
`H = weber_term(q,p; …) + zollner_term(…)`.

# Arguments
- `q_vars`, `p_vars`: Phase-space symbolic variables, each length `n_particles*dims`.

# Keywords
- `masses`: Per-particle mass symbolic variables (length `n_particles`).
- `charges`: Per-particle charge symbolic variables (length `n_particles`).
- `c`: Speed-of-light symbolic variable.
- `kappas`: Per-pair κ coupling symbolic variables (length `n_particles*(n_particles-1)/2`,
  ordered by `i<j` and indexed via `_pair_index`).
- `n_particles::Int`, `dims::Int`: Problem shape.
"""
function weber_term(
    q_vars::AbstractVector,
    p_vars::AbstractVector;
    masses::AbstractVector,
    charges::AbstractVector,
    c,
    kappas::AbstractVector,
    n_particles::Int,
    dims::Int,
)
    H = zero(eltype(q_vars))

    @inbounds for i = 1:n_particles
        p_start = (i - 1) * dims + 1
        p_end = i * dims
        p_squared = sum(p_vars[p_start:p_end] .^ 2)
        H = H + p_squared / (2 * masses[i])
    end

    c_squared = c^2
    @inbounds for i = 1:n_particles
        for j = (i+1):n_particles
            qi_start = (i - 1) * dims + 1
            qj_start = (j - 1) * dims + 1

            pi_start = (i - 1) * dims + 1
            pj_start = (j - 1) * dims + 1

            r_squared = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                r_squared = r_squared + dq^2
            end
            r = sqrt(r_squared)

            r_dot_v = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                dv = p_vars[pi_start+d-1] / masses[i] - p_vars[pj_start+d-1] / masses[j]
                r_dot_v = r_dot_v + dq * dv
            end
            r_dot = r_dot_v / r

            pair_idx = _pair_index(i, j, n_particles)
            kappa = kappas[pair_idx]
            k = kappa * charges[i] * charges[j]
            U_ij = k / r * (1 - r_dot^2 / (2 * c_squared))
            H = H + U_ij
        end
    end

    return H
end

# Concrete per-pair numeric decomposition of the Weber pair potential, matching
# the symbolic `weber_term` body. The closure is attached to the `:weber`
# NamedTerm by `HamiltonianSystem(n, dims)` and queried by statistics.
#
# Returns a NamedTuple `(coulomb, velocity, rdot, r)` where
#   coulomb  = κ·q_i·q_j / r                                    (Coulomb part)
#   velocity = −coulomb · ṙ² / (2c²)                            (Weber velocity part)
#   rdot     = ṙ = (q_i − q_j)·(v_i − v_j) / r                  (radial velocity)
#   r        = |q_i − q_j|                                      (pair distance)
# `coulomb + velocity` is the full Weber-pair contribution to H (κ-scaled).
function _weber_pair_decomposition(
    i::Int,
    j::Int,
    q::AbstractVector{Float64},
    p::AbstractVector{Float64},
    params::AbstractVector{Float64},
    kappas::AbstractVector{Float64},
    n_particles::Int,
    dims::Int,
)
    @inbounds mi = params[i]
    @inbounds mj = params[j]
    @inbounds qi = params[n_particles+i]
    @inbounds qj = params[n_particles+j]
    @inbounds c = params[2*n_particles+1]
    @inbounds κ = kappas[_pair_index(i, j, n_particles)]

    qi_start = (i - 1) * dims + 1
    qj_start = (j - 1) * dims + 1

    r_squared = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        r_squared += dq * dq
    end
    r = sqrt(r_squared)

    r_dot_v = 0.0
    @inbounds for d = 0:(dims-1)
        dq = q[qi_start+d] - q[qj_start+d]
        dv = p[qi_start+d] / mi - p[qj_start+d] / mj
        r_dot_v += dq * dv
    end
    rdot = r_dot_v / r

    coulomb = κ * qi * qj / r
    velocity = -coulomb * rdot^2 / (2 * c^2)
    return (coulomb = coulomb, velocity = velocity, rdot = rdot, r = r)
end
