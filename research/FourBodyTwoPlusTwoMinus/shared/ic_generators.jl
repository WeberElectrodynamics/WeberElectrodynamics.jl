"""
Initial condition generators for the 4-body Weber Hamiltonian with 2 positive
and 2 negative charges. All generators return `(q0, p0, masses, charges)` with
flat length-8 (2D) or length-12 (3D) q/p vectors.

Conventions:
- Particles labeled 1,2 are positive (+q), 3,4 are negative (-q).
- All generators are center-of-mass and zero-total-momentum by construction.
- `energy_fraction` η ∈ (0,1) sets T₀ = η·|U₀| at the initial instant
  (bound iff η < 1 for Coulomb-like configurations).
"""

using LinearAlgebra

# --------------------------------------------------------------------------
# Utility: given position matrix X (n×d), charges, compute instantaneous U.
# --------------------------------------------------------------------------
function _coulomb_U(X::AbstractMatrix, charges::AbstractVector)
    n = size(X, 1)
    U = 0.0
    @inbounds for i = 1:n, j = (i+1):n
        r = norm(view(X, i, :) .- view(X, j, :))
        U += charges[i] * charges[j] / r
    end
    return U
end

# --------------------------------------------------------------------------
# Remove COM and enforce zero total momentum.
# --------------------------------------------------------------------------
function _zero_com!(X::AbstractMatrix, masses::AbstractVector)
    com = (masses' * X) ./ sum(masses)
    X .-= com
    return X
end

function _zero_total_p!(P::AbstractMatrix)
    P .-= sum(P, dims = 1) ./ size(P, 1)
    return P
end

function _flatten(X::AbstractMatrix)
    n, d = size(X)
    out = Vector{Float64}(undef, n * d)
    @inbounds for i = 1:n, k = 1:d
        out[(i-1)*d+k] = X[i, k]
    end
    return out
end

# --------------------------------------------------------------------------
# Alternating square: +,−,+,− at vertices of a square with side s.
# Labeling: 1 at (+s/2, +s/2), 3 at (-s/2, +s/2), 2 at (-s/2, -s/2), 4 at (+s/2, -s/2)
# (so 1,2 are +; 3,4 are −; 1–3 and 2–4 are edges; 1–2 and 3–4 are diagonals).
# velocity_mode ∈ (:breathing, :rotating, :static)
# --------------------------------------------------------------------------
function alternating_square(;
    side::Float64 = 1.0,
    q::Float64 = 1.0,
    m::Float64 = 1.0,
    energy_fraction::Float64 = 0.25,
    velocity_mode::Symbol = :rotating,
)
    s = side / 2
    X = [
        s  s;
        -s -s;
        -s  s;
        s -s;
    ]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    U = _coulomb_U(X, charges)
    T_target = energy_fraction * abs(U)
    # per-particle kinetic
    speed = sqrt(2 * T_target / (4 * m))
    P = zeros(4, 2)
    if velocity_mode == :rotating
        # tangential velocity around COM
        for i = 1:4
            r = X[i, :]
            t = [-r[2], r[1]] / max(norm(r), eps())
            P[i, :] = m * speed * t
        end
    elseif velocity_mode == :breathing
        for i = 1:4
            r = X[i, :]
            P[i, :] = -m * speed * r / max(norm(r), eps())
        end
    end
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# --------------------------------------------------------------------------
# Two dimers: (+−)(+−) dyads separated by R along x, each dyad of length a.
# Each dyad orbits internally (bound 2-body) with tangential momenta.
# --------------------------------------------------------------------------
function two_dimers(;
    dyad_length::Float64 = 0.2,
    separation::Float64 = 2.0,
    q::Float64 = 1.0,
    m::Float64 = 1.0,
    intra_fraction::Float64 = 0.4,
    inter_velocity::Float64 = 0.0,
)
    a = dyad_length
    R = separation
    X = [
        R/2  a/2;    # +1
        -R/2  a/2;   # +2
        R/2 -a/2;    # -1
        -R/2 -a/2;   # -2
    ]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    # Intra-dyad tangential speed for each (+,−) pair; using Coulomb circular estimate.
    # Pairs: (1,3) right dyad, (2,4) left dyad. r = a.
    mu = m / 2
    v_circ = sqrt(intra_fraction * (q * q) / (mu * a))
    P = zeros(4, 2)
    # Right dyad: + at top (+y), − at bottom (-y); rotate about dyad center (R/2,0)
    P[1, :] = m * v_circ * [ 1.0, 0.0]
    P[3, :] = m * v_circ * [-1.0, 0.0]
    # Left dyad
    P[2, :] = m * v_circ * [-1.0, 0.0]
    P[4, :] = m * v_circ * [ 1.0, 0.0]
    # Optional relative motion between dyads
    P[1, 1] += m * inter_velocity
    P[3, 1] += m * inter_velocity
    P[2, 1] -= m * inter_velocity
    P[4, 1] -= m * inter_velocity
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# --------------------------------------------------------------------------
# Linear ABAB chain: + − + − collinear along x, equispaced.
# --------------------------------------------------------------------------
function linear_abab(;
    spacing::Float64 = 1.0,
    q::Float64 = 1.0,
    m::Float64 = 1.0,
    transverse_kick::Float64 = 0.0,
)
    d = spacing
    X = [
        -1.5d  0.0;   # +
        0.5d  0.0;   # +
        -0.5d  0.0;   # −
        1.5d  0.0;   # −
    ]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    P = zeros(4, 2)
    # Transverse staggered kick (optional)
    if transverse_kick != 0.0
        P[:, 2] = m * transverse_kick * [1.0, -1.0, 1.0, -1.0]
    end
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# --------------------------------------------------------------------------
# Linear AABB chain: + + − − collinear.
# --------------------------------------------------------------------------
function linear_aabb(;
    spacing::Float64 = 1.0,
    q::Float64 = 1.0,
    m::Float64 = 1.0,
    transverse_kick::Float64 = 0.0,
)
    d = spacing
    X = [
        -1.5d  0.0;   # +
        -0.5d  0.0;   # +
        0.5d  0.0;   # −
        1.5d  0.0;   # −
    ]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    P = zeros(4, 2)
    if transverse_kick != 0.0
        P[:, 2] = m * transverse_kick * [1.0, -1.0, 1.0, -1.0]
    end
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# --------------------------------------------------------------------------
# Rhombus: vertices at (±a, 0), (0, ±b) with charges + − + − cyclically.
# --------------------------------------------------------------------------
function rhombus(;
    a::Float64 = 1.0,
    b::Float64 = 0.6,
    q::Float64 = 1.0,
    m::Float64 = 1.0,
    energy_fraction::Float64 = 0.25,
    velocity_mode::Symbol = :rotating,
)
    X = [
        a  0.0;
        -a  0.0;
        0.0  b;
        0.0 -b;
    ]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    U = _coulomb_U(X, charges)
    T_target = energy_fraction * abs(U)
    speed = sqrt(2 * T_target / (4 * m))
    P = zeros(4, 2)
    if velocity_mode == :rotating
        for i = 1:4
            r = X[i, :]
            t = [-r[2], r[1]] / max(norm(r), eps())
            P[i, :] = m * speed * t
        end
    end
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# --------------------------------------------------------------------------
# Planetary-atom-with-moon: ++ nucleus pair at small separation, with a −
# orbiting close and another − far away.
# --------------------------------------------------------------------------
function planetary_atom_with_moon(;
    nucleus_sep::Float64 = 0.05,
    inner_orbit::Float64 = 0.3,
    outer_orbit::Float64 = 2.0,
    q::Float64 = 1.0,
    m::Float64 = 1.0,
)
    X = [
        -nucleus_sep/2  0.0;
        nucleus_sep/2  0.0;
        inner_orbit  0.0;
        outer_orbit  0.0;
    ]
    masses = fill(m, 4)
    charges = [q, q, -q, -q]
    # Inner electron circular around combined nucleus at (0,0) treated as 2q.
    v_in = sqrt(2 * q * q / (m * inner_orbit))
    v_out = sqrt(2 * q * q / (m * outer_orbit))
    P = zeros(4, 2)
    P[3, :] = m * v_in * [0.0, 1.0]
    P[4, :] = m * v_out * [0.0, 1.0]
    _zero_com!(X, masses)
    _zero_total_p!(P)
    return (_flatten(X), _flatten(P), masses, charges)
end

# Export
const IC_GENERATORS = Dict(
    :alternating_square => alternating_square,
    :two_dimers => two_dimers,
    :linear_abab => linear_abab,
    :linear_aabb => linear_aabb,
    :rhombus => rhombus,
    :planetary_atom_with_moon => planetary_atom_with_moon,
)
