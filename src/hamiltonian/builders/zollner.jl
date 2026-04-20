"""
    zollner_term(q_vars, p_vars;
                 masses, charges, c, kappas,
                 n_particles, dims) -> Num

Build the symbolic Zöllner correction to the Weber Hamiltonian.

Returns `Σ_{i<j} (kappas[pair_idx] - 1) · U_weber(i,j)`, where `U_weber(i,j)`
is the standard Weber pair potential `q_i q_j / r · (1 − ṙ² / (2c²))`. This
is the extra potential relative to the κ=1 (pure Weber) case, so

    weber_term(…; kappas = κ)
    == weber_term(…; kappas = ones) + zollner_term(…; kappas = κ)

up to Symbolics.jl expression rewriting. When all κ = 1 the correction is
identically zero.

The physical Zöllner mismatch (κ_ij = 1+a for unlike-sign pairs, 1 otherwise)
is expressed by the concrete κ values injected via `params` at solve time;
see `_compute_zollner_kappas`.

# Arguments
- `q_vars`, `p_vars`: Phase-space symbolic variables.

# Keywords
- `masses`: Per-particle mass symbolic variables (enter through `ṙ = r̂ · v`,
  with `v_i = p_i / m_i`).
- `charges`: Per-particle charge symbolic variables.
- `c`: Speed-of-light symbolic variable.
- `kappas`: Per-pair κ coupling symbolic variables
  (length `n_particles*(n_particles-1)/2`, ordered by `i<j` and indexed via
  `_pair_index`).
- `n_particles::Int`, `dims::Int`: Problem shape.
"""
function zollner_term(
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
            extra = kappas[pair_idx] - 1
            U_ij = extra * charges[i] * charges[j] / r * (1 - r_dot^2 / (2 * c_squared))
            H = H + U_ij
        end
    end

    return H
end
