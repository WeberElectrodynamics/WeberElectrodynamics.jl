# =============================================================================
# Exact canonical Weber system
# =============================================================================
#
# The Weber Lagrangian
#
#   L = Σ_i ½ m_i |v_i|² − Σ_{i<j} q_i q_j / r_ij · (1 + ṙ_ij²/(2c²))
#
# is velocity dependent, so its canonical momentum p_i = ∂L/∂v_i is NOT m_i v_i:
#
#   p_i = m_i v_i − Σ_{j≠i} (q_i q_j / c²) · ṙ_ij / r_ij² · (r_i − r_j)
#
# Writing k_ij = q_i q_j / (c² r_ij), this is the linear system p = M(q) v with
# the symmetric canonical mass matrix
#
#   M_ii = m_i I − Σ_{j≠i} k_ij r̂_ij r̂_ijᵀ,      M_ij = + k_ij r̂_ij r̂_ijᵀ.
#
# Inverting M directly is impractical symbolically, but the correction is a sum
# of n(n−1)/2 rank-one pair terms, so the inverse reduces to a linear system in
# the *pair radial velocities* alone:
#
#   ṙ_b = s_b + Σ_a G_ba k_a ṙ_a,   i.e.   (I − G K) ṙ = s
#
# where, for pairs a = (i,j) and b = (k,l),
#
#   s_b  = r̂_b · (p_k/m_k − p_l/m_l)          ("naive" radial rate)
#   G_ba = (r̂_b · r̂_a) · (δ_ki/m_k − δ_kj/m_k − δ_li/m_l + δ_lj/m_l)
#
# Once ṙ is known the physical velocities follow explicitly,
#
#   v_i = p_i/m_i + (1/m_i) Σ_{j≠i} k_ij ṙ_ij r̂_ij,
#
# and the exact canonical Hamiltonian (the Legendre transform of L) is
#
#   H = Σ_i |p_i|²/(2m_i) + ½ Σ_{i<j} k_ij ṙ_ij s_ij + Σ_{i<j} q_i q_j / r_ij
#
# which equals the velocity-space energy Σ_i ½ m_i |v_i|² + Σ_{i<j} U_ij.
# The canonical equations are then
#
#   q̇_i = v_i
#   ṗ_i = Σ_{j≠i} q_i q_j / r_ij² · [ r̂_ij (1 + 3ṙ_ij²/(2c²)) − ṙ_ij (v_i − v_j)/c² ]
#
# Every formula above is verified independently in
# `papers/Computational-Weber-Electrodynamics/verify_formulas.py`.
# =============================================================================

"""
    WeberCriticalRadiusError(pair, r, rho)

Thrown when the canonical mass matrix of the Weber system is singular, i.e.
when the linear system relating canonical momenta to physical velocities cannot
be solved.

For a like-charge two-particle system this happens exactly at Weber's critical
radius `rho = q₁q₂/(μc²)`, where `μ` is the reduced mass: there the effective
radial inertia `μ − q₁q₂/(rc²)` vanishes and a canonical momentum no longer
determines a finite physical velocity. Below `rho` the effective inertia is
negative and the dynamics remain well defined; only the crossing itself is
singular.

# Fields
- `pair::Tuple{Int,Int}`: Particle pair whose elimination step failed, or
  `(0, 0)` when the failure cannot be attributed to a single pair.
- `r::Float64`: Separation of that pair at the failure point.
- `rho::Float64`: Critical radius of that pair, or `NaN` when undefined.
"""
struct WeberCriticalRadiusError <: Exception
    pair::Tuple{Int,Int}
    r::Float64
    rho::Float64
end

function Base.showerror(io::IO, e::WeberCriticalRadiusError)
    print(io, "WeberCriticalRadiusError: the canonical Weber mass matrix is singular")
    if e.pair != (0, 0)
        print(io, " at pair $(e.pair) with r = $(e.r)")
        isfinite(e.rho) && print(io, " (critical radius rho = $(e.rho))")
    end
    print(
        io,
        ". Canonical momenta do not determine finite physical velocities there; ",
        "reduce the step size, enable regularization, or start away from the ",
        "critical radius.",
    )
end

"""
    WeberWorkspace(n_particles, dims)

Pre-allocated workspace for one exact canonical Weber evaluation.

Holds the pair table, the per-pair geometry, the `n_pairs × n_pairs` radial
system, and the recovered physical velocities. A workspace is **not**
thread-safe: each compiled Weber function owns its own instance.

# Fields
- `n::Int`, `dims::Int`, `npairs::Int`: Problem shape (`npairs = n(n−1)/2`).
- `pair_i`, `pair_j::Vector{Int}`: Pair table in `pair_indices` order.
- `r::Vector{Float64}`: Pair separations `r_ij`.
- `rhat::Matrix{Float64}`: `dims × npairs` unit separation vectors.
- `kco::Vector{Float64}`: Pair coefficients `q_i q_j / (c² r_ij)`.
- `s::Vector{Float64}`: Naive radial rates `r̂_ij · (p_i/m_i − p_j/m_j)`.
- `rdot::Vector{Float64}`: Physical pair radial velocities `ṙ_ij`.
- `v::Vector{Float64}`: Physical particle velocities, length `n × dims`.
- `A::Matrix{Float64}`: Radial system matrix `I − G K`.
"""
struct WeberWorkspace
    n::Int
    dims::Int
    npairs::Int
    pair_i::Vector{Int}
    pair_j::Vector{Int}
    r::Vector{Float64}
    rhat::Matrix{Float64}
    kco::Vector{Float64}
    s::Vector{Float64}
    rdot::Vector{Float64}
    v::Vector{Float64}
    A::Matrix{Float64}
end

function WeberWorkspace(n::Int, dims::Int)
    npairs = n * (n - 1) ÷ 2
    pair_i = Vector{Int}(undef, npairs)
    pair_j = Vector{Int}(undef, npairs)
    idx = 1
    @inbounds for i = 1:n
        for j = (i+1):n
            pair_i[idx] = i
            pair_j[idx] = j
            idx += 1
        end
    end
    return WeberWorkspace(
        n,
        dims,
        npairs,
        pair_i,
        pair_j,
        Vector{Float64}(undef, npairs),
        Matrix{Float64}(undef, dims, npairs),
        Vector{Float64}(undef, npairs),
        Vector{Float64}(undef, npairs),
        Vector{Float64}(undef, npairs),
        Vector{Float64}(undef, n * dims),
        Matrix{Float64}(undef, npairs, npairs),
    )
end

# Dense LU with partial pivoting, specialised for the small symmetric-structure
# radial system. Solves `A x = b` in place, leaving the solution in `b`.
# Throws `WeberCriticalRadiusError` when a pivot collapses.
function _weber_linear_solve!(ws::WeberWorkspace, params::AbstractVector{<:Real})
    P = ws.npairs
    A = ws.A
    b = ws.rdot

    scale = 1.0
    @inbounds for j = 1:P, i = 1:P
        a = abs(A[i, j])
        a > scale && (scale = a)
    end
    tol = 1e-12 * scale

    @inbounds for k = 1:P
        piv = k
        amax = abs(A[k, k])
        for i = (k+1):P
            a = abs(A[i, k])
            if a > amax
                amax = a
                piv = i
            end
        end
        if amax <= tol
            throw(_weber_singular_error(ws, params, k))
        end
        if piv != k
            for j = k:P
                A[k, j], A[piv, j] = A[piv, j], A[k, j]
            end
            b[k], b[piv] = b[piv], b[k]
        end
        akk = A[k, k]
        for i = (k+1):P
            f = A[i, k] / akk
            f == 0.0 && continue
            A[i, k] = 0.0
            for j = (k+1):P
                A[i, j] -= f * A[k, j]
            end
            b[i] -= f * b[k]
        end
    end

    @inbounds for k = P:-1:1
        acc = b[k]
        for j = (k+1):P
            acc -= A[k, j] * b[j]
        end
        b[k] = acc / A[k, k]
    end

    return nothing
end

function _weber_singular_error(ws::WeberWorkspace, params::AbstractVector{<:Real}, k::Int)
    n = ws.n
    # Attribute the failure to the pair with the smallest |r − rho| when that is
    # meaningful; otherwise report the pair whose elimination step collapsed.
    a = clamp(k, 1, max(ws.npairs, 1))
    ws.npairs == 0 && return WeberCriticalRadiusError((0, 0), NaN, NaN)
    i = ws.pair_i[a]
    j = ws.pair_j[a]
    mi = Float64(params[i])
    mj = Float64(params[j])
    qi = Float64(params[n+i])
    qj = Float64(params[n+j])
    c = Float64(params[2n+1])
    mu = mi * mj / (mi + mj)
    rho = qi * qj / (mu * c^2)
    return WeberCriticalRadiusError((i, j), ws.r[a], rho)
end

"""
    _weber_state!(ws, q, p, params) -> ws

Populate `ws` with the exact canonical Weber state at `(q, p)`: pair geometry,
naive radial rates, physical pair radial velocities, and physical particle
velocities.

This is the single routine through which every Weber path — equations of
motion, Hamiltonian, energy statistics, force statistics, and the Makie
dashboard — obtains physical velocities. Nothing in the package may reconstruct
a velocity as `p_i/m_i`.
"""
function _weber_state!(
    ws::WeberWorkspace,
    q::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    n = ws.n
    d = ws.dims
    P = ws.npairs
    @inbounds c2 = Float64(params[2n+1])^2

    # Physical velocity starts from the kinetic term and picks up pair corrections.
    @inbounds for i = 1:n
        mi = Float64(params[i])
        base = (i - 1) * d
        for t = 1:d
            ws.v[base+t] = p[base+t] / mi
        end
    end

    P == 0 && return ws

    # Pair geometry, naive radial rates, pair coefficients.
    @inbounds for a = 1:P
        i = ws.pair_i[a]
        j = ws.pair_j[a]
        mi = Float64(params[i])
        mj = Float64(params[j])
        bi = (i - 1) * d
        bj = (j - 1) * d

        r2 = 0.0
        for t = 1:d
            dq = q[bi+t] - q[bj+t]
            ws.rhat[t, a] = dq
            r2 += dq * dq
        end
        r = sqrt(r2)
        ws.r[a] = r
        for t = 1:d
            ws.rhat[t, a] /= r
        end

        acc = 0.0
        for t = 1:d
            acc += ws.rhat[t, a] * (p[bi+t] / mi - p[bj+t] / mj)
        end
        ws.s[a] = acc
        ws.rdot[a] = acc

        qi = Float64(params[n+i])
        qj = Float64(params[n+j])
        ws.kco[a] = qi * qj / (c2 * r)
    end

    # Radial system A = I − G K.
    @inbounds for b = 1:P
        ib = ws.pair_i[b]
        jb = ws.pair_j[b]
        inv_mib = 1.0 / Float64(params[ib])
        inv_mjb = 1.0 / Float64(params[jb])
        for a = 1:P
            ia = ws.pair_i[a]
            ja = ws.pair_j[a]
            coup =
                (ib == ia ? inv_mib : 0.0) - (ib == ja ? inv_mib : 0.0) -
                (jb == ia ? inv_mjb : 0.0) + (jb == ja ? inv_mjb : 0.0)
            if coup == 0.0
                ws.A[b, a] = (a == b) ? 1.0 : 0.0
            else
                dot_ba = 0.0
                for t = 1:d
                    dot_ba += ws.rhat[t, b] * ws.rhat[t, a]
                end
                ws.A[b, a] = (a == b ? 1.0 : 0.0) - dot_ba * coup * ws.kco[a]
            end
        end
    end

    _weber_linear_solve!(ws, params)

    # v_i = p_i/m_i + (1/m_i) Σ_a σ_ia k_a ṙ_a r̂_a
    @inbounds for a = 1:P
        i = ws.pair_i[a]
        j = ws.pair_j[a]
        inv_mi = 1.0 / Float64(params[i])
        inv_mj = 1.0 / Float64(params[j])
        w = ws.kco[a] * ws.rdot[a]
        bi = (i - 1) * d
        bj = (j - 1) * d
        for t = 1:d
            contrib = w * ws.rhat[t, a]
            ws.v[bi+t] += contrib * inv_mi
            ws.v[bj+t] -= contrib * inv_mj
        end
    end

    return ws
end

"""
    weber_dq_dt!(out, ws, q, p, params)

Write the first canonical equation `q̇_i = v_i(q, p)` into `out`.

The coordinate rate is the *physical* velocity recovered from the canonical
momenta, not `p_i/m_i`.
"""
function weber_dq_dt!(
    out::AbstractVector{<:Real},
    ws::WeberWorkspace,
    q::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    _weber_state!(ws, q, p, params)
    @inbounds copyto!(out, ws.v)
    return out
end

"""
    weber_dp_dt!(out, ws, q, p, params)

Write the second canonical equation into `out`:

`ṗ_i = Σ_{j≠i} q_i q_j / r_ij² · [ r̂_ij (1 + 3ṙ_ij²/(2c²)) − ṙ_ij (v_i − v_j)/c² ]`

Both `1/c²` terms use the physical radial velocity and physical relative
velocity obtained from the canonical momenta.
"""
function weber_dp_dt!(
    out::AbstractVector{<:Real},
    ws::WeberWorkspace,
    q::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    _weber_state!(ws, q, p, params)

    n = ws.n
    d = ws.dims
    @inbounds c2 = Float64(params[2n+1])^2

    fill!(out, 0.0)

    @inbounds for a = 1:ws.npairs
        i = ws.pair_i[a]
        j = ws.pair_j[a]
        qi = Float64(params[n+i])
        qj = Float64(params[n+j])
        k = qi * qj
        k == 0.0 && continue

        r = ws.r[a]
        rdot = ws.rdot[a]
        pref = k / (r * r)
        radial = 1.0 + 3.0 * rdot * rdot / (2.0 * c2)

        bi = (i - 1) * d
        bj = (j - 1) * d
        for t = 1:d
            dv = ws.v[bi+t] - ws.v[bj+t]
            term = pref * (ws.rhat[t, a] * radial - rdot * dv / c2)
            out[bi+t] += term
            out[bj+t] -= term
        end
    end

    return out
end

"""
    weber_hamiltonian(ws, q, p, params) -> Float64

Evaluate the exact canonical Weber Hamiltonian

`H = Σ_i |p_i|²/(2m_i) + ½ Σ_{i<j} k_ij ṙ_ij s_ij + Σ_{i<j} q_i q_j / r_ij`

with `k_ij = q_i q_j/(c² r_ij)`, `s_ij` the naive radial rate, and `ṙ_ij` the
physical radial velocity. This is the Legendre transform of the Weber
Lagrangian and equals the velocity-space energy `Σ ½ m_i |v_i|² + Σ U_ij`.
"""
function weber_hamiltonian(
    ws::WeberWorkspace,
    q::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    params::AbstractVector{<:Real},
)
    _weber_state!(ws, q, p, params)

    n = ws.n
    d = ws.dims

    H = 0.0
    @inbounds for i = 1:n
        mi = Float64(params[i])
        base = (i - 1) * d
        acc = 0.0
        for t = 1:d
            acc += p[base+t]^2
        end
        H += acc / (2 * mi)
    end

    @inbounds for a = 1:ws.npairs
        i = ws.pair_i[a]
        j = ws.pair_j[a]
        qi = Float64(params[n+i])
        qj = Float64(params[n+j])
        H += 0.5 * ws.kco[a] * ws.rdot[a] * ws.s[a] + qi * qj / ws.r[a]
    end

    return H
end

# Per-state numeric decomposition of the Weber pair energies, attached to the
# `:weber` NamedTerm by `HamiltonianSystem(n, dims)` and consumed by the energy
# statistics. Returns per-pair vectors in `pair_indices` order:
#
#   coulomb  = q_i q_j / r                          (Coulomb part)
#   velocity = −coulomb · ṙ²/(2c²)                  (Weber velocity part)
#   rdot     = ṙ_ij, the PHYSICAL pair radial velocity
#   r        = |q_i − q_j|
#
# `Σ (coulomb + velocity) + Σ ½ m_i |v_i|²` reproduces the compiled Hamiltonian
# exactly, which is what `hamiltonian_validation_error` cross-checks.
function _weber_pair_decomposition(
    ws::WeberWorkspace,
    q::AbstractVector{Float64},
    p::AbstractVector{Float64},
    params::AbstractVector{Float64},
)
    _weber_state!(ws, q, p, params)

    n = ws.n
    P = ws.npairs
    @inbounds c2 = Float64(params[2n+1])^2

    coulomb = Vector{Float64}(undef, P)
    velocity = Vector{Float64}(undef, P)
    rdot = Vector{Float64}(undef, P)
    rvec = Vector{Float64}(undef, P)

    @inbounds for a = 1:P
        i = ws.pair_i[a]
        j = ws.pair_j[a]
        cb = Float64(params[n+i]) * Float64(params[n+j]) / ws.r[a]
        coulomb[a] = cb
        velocity[a] = -cb * ws.rdot[a]^2 / (2 * c2)
        rdot[a] = ws.rdot[a]
        rvec[a] = ws.r[a]
    end

    return (coulomb = coulomb, velocity = velocity, rdot = rdot, r = rvec)
end

"""
    physical_velocities(q, p, params; n_particles, dims) -> Vector{Float64}

Recover the physical particle velocities `v` from a canonical Weber state
`(q, p)`.

Because the Weber Lagrangian is velocity dependent, `p_i` is **not** `m_i v_i`
whenever any pair has nonzero radial velocity. Use this function — never
`p ./ masses` — to obtain velocities, kinetic energies, or radial rates from a
solution.

`params` uses the standard layout `[m₁…mₙ, q₁…qₙ, c]`. The returned vector is
flattened in the same `[v₁ₓ, v₁ᵧ, …, vₙₓ, vₙᵧ, …]` order as `q` and `p`.

Allocates a fresh workspace per call; for tight loops the package uses the
in-place internal path instead.
"""
function physical_velocities(
    q::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    params::AbstractVector{<:Real};
    n_particles::Int,
    dims::Int,
)
    dof = n_particles * dims
    length(q) == dof || throw(ArgumentError("q must have length $dof, got $(length(q))"))
    length(p) == dof || throw(ArgumentError("p must have length $dof, got $(length(p))"))
    length(params) == 2 * n_particles + 1 || throw(
        ArgumentError(
            "params must have length $(2 * n_particles + 1), got $(length(params))",
        ),
    )
    ws = WeberWorkspace(n_particles, dims)
    _weber_state!(ws, q, p, params)
    return copy(ws.v)
end

# Physical kinetic energy Σ ½ m_i |v_i|² at the current state.
function _weber_kinetic_energy(
    ws::WeberWorkspace,
    q::AbstractVector{Float64},
    p::AbstractVector{Float64},
    params::AbstractVector{Float64},
)
    _weber_state!(ws, q, p, params)
    n = ws.n
    d = ws.dims
    KE = 0.0
    @inbounds for i = 1:n
        mi = Float64(params[i])
        base = (i - 1) * d
        acc = 0.0
        for t = 1:d
            acc += ws.v[base+t]^2
        end
        KE += 0.5 * mi * acc
    end
    return KE
end
