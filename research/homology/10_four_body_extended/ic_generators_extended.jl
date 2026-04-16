"""
Extended IC generators for 4-body Weber beyond 2+/2-.
Covers: 3+/1- configurations, unequal masses, unequal charges,
and c-variation of known candidates.
"""

using LinearAlgebra

# --------------------------------------------------------------------------
# Utility (duplicated from shared for standalone use)
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

function _zero_com!(X::AbstractMatrix, masses::AbstractVector)
    com = (masses' * X) ./ sum(masses)
    X .-= com
    return X
end

function _zero_total_p!(P::AbstractMatrix, masses::AbstractVector)
    # Zero total momentum: subtract uniform velocity
    total_p = vec(sum(P, dims = 1))
    M = sum(masses)
    v_com = total_p ./ M
    for i in 1:size(P, 1)
        P[i, :] .-= masses[i] .* v_com
    end
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

# ==========================================================================
# 1. 3+/1- configurations
# ==========================================================================

"""
    triangular_trap_3p1m(; R, q_pos, q_neg, m, energy_fraction, velocity_mode, dims)

3 positive charges at vertices of equilateral triangle (radius R from center),
1 negative charge at center. 2D or 3D.
Particles: 1,2,3 are positive, 4 is negative.
"""
function triangular_trap_3p1m(;
    R::Float64 = 2.0,
    q_pos::Float64 = 1.0,
    q_neg::Float64 = -1.0,
    masses::Vector{Float64} = ones(4),
    energy_fraction::Float64 = 0.25,
    velocity_mode::Symbol = :rotating,
    dims::Int = 2,
)
    angles = [0.0, 2pi/3, 4pi/3]
    if dims == 2
        X = zeros(4, 2)
        for i in 1:3
            X[i, 1] = R * cos(angles[i])
            X[i, 2] = R * sin(angles[i])
        end
        # particle 4 at origin
    elseif dims == 3
        X = zeros(4, 3)
        for i in 1:3
            X[i, 1] = R * cos(angles[i])
            X[i, 2] = R * sin(angles[i])
        end
    end
    charges = [q_pos, q_pos, q_pos, q_neg]
    _zero_com!(X, masses)

    U = _coulomb_U(X, charges)
    T_target = energy_fraction * abs(U)
    P = zeros(4, dims)

    if velocity_mode == :rotating && abs(T_target) > 0
        # Tangential velocities for outer particles
        total_ke = 0.0
        dirs = zeros(3, dims)
        for i in 1:3
            r = X[i, :]
            if dims == 2
                dirs[i, :] = [-r[2], r[1]] / max(norm(r), eps())
            else
                dirs[i, :] = [-r[2], r[1], 0.0] / max(norm(r[1:2]), eps())
            end
        end
        # Equal speed for outer particles
        speed = sqrt(2 * T_target / (3 * masses[1]))  # approx
        for i in 1:3
            P[i, :] = masses[i] * speed * dirs[i, :]
        end
    elseif velocity_mode == :breathing && abs(T_target) > 0
        speed = sqrt(2 * T_target / (3 * masses[1]))
        for i in 1:3
            r = X[i, :]
            P[i, :] = -masses[i] * speed * r / max(norm(r), eps())
        end
    end

    _zero_total_p!(P, masses)
    return (_flatten(X), _flatten(P), masses, charges)
end

"""
    linear_3p1m(; spacing, q_pos, q_neg, masses, transverse_kick, dims)

Linear chain: + + - + along x-axis. The negative charge is between two positives.
"""
function linear_3p1m(;
    spacing::Float64 = 1.5,
    q_pos::Float64 = 1.0,
    q_neg::Float64 = -1.0,
    masses::Vector{Float64} = ones(4),
    transverse_kick::Float64 = 0.1,
    dims::Int = 2,
)
    d = spacing
    if dims == 2
        X = [
            -1.5d  0.0;   # +
             0.5d  0.0;   # +
            -0.5d  0.0;   # -  (between first two +)
             1.5d  0.0;   # +
        ]
        # Wait -- we want + + - + order.
        # Positions: -1.5d, -0.5d, 0.5d, 1.5d
        # Charges:      +,    +,    -,    +
        X = [
            -1.5d  0.0;
            -0.5d  0.0;
             0.5d  0.0;
             1.5d  0.0;
        ]
    else
        X = [
            -1.5d  0.0  0.0;
            -0.5d  0.0  0.0;
             0.5d  0.0  0.0;
             1.5d  0.0  0.0;
        ]
    end
    charges = [q_pos, q_pos, q_neg, q_pos]
    _zero_com!(X, masses)

    P = zeros(4, dims)
    if transverse_kick != 0.0
        P[:, 2] = [masses[i] * transverse_kick * (-1)^(i+1) for i in 1:4]
    end
    _zero_total_p!(P, masses)
    return (_flatten(X), _flatten(P), masses, charges)
end

"""
    tetrahedral_3p1m(; R, q_pos, q_neg, masses, energy_fraction, velocity_mode)

3D tetrahedron: 3 positive charges at base triangle, 1 negative at apex.
"""
function tetrahedral_3p1m(;
    R::Float64 = 2.0,
    q_pos::Float64 = 1.0,
    q_neg::Float64 = -1.0,
    masses::Vector{Float64} = ones(4),
    energy_fraction::Float64 = 0.25,
    velocity_mode::Symbol = :rotating,
)
    # Regular tetrahedron with edge length a = R*sqrt(3)
    # Base triangle in z=0 plane, apex above
    h = R * sqrt(2/3) * sqrt(3)  # height
    X = zeros(4, 3)
    angles = [0.0, 2pi/3, 4pi/3]
    for i in 1:3
        X[i, 1] = R * cos(angles[i])
        X[i, 2] = R * sin(angles[i])
        X[i, 3] = -h/4  # shift so COM is at 0
    end
    X[4, 3] = 3h/4
    charges = [q_pos, q_pos, q_pos, q_neg]
    _zero_com!(X, masses)

    U = _coulomb_U(X, charges)
    T_target = energy_fraction * abs(U)
    P = zeros(4, 3)

    if velocity_mode == :rotating && abs(T_target) > 0
        # Tangential velocities for base particles around z-axis
        speed = sqrt(2 * T_target / (3 * masses[1]))
        for i in 1:3
            r2d = X[i, 1:2]
            tang = [-r2d[2], r2d[1], 0.0] / max(norm(r2d), eps())
            P[i, :] = masses[i] * speed * tang
        end
    end

    _zero_total_p!(P, masses)
    return (_flatten(X), _flatten(P), masses, charges)
end

# ==========================================================================
# 2. Symmetric double-orbiter (from prior study, parameterized for c-variation)
# ==========================================================================

"""
    symmetric_double_orbiter(; r_pp, R, orb, z_kick, masses, charges, dims)

(+,+) on x-axis at +/-r_pp/2, (-,-) on y-axis at (0,+/-R).
Negative pair gets tangential momenta +/-orb*sqrt(2/R).
Optional z_kick for 3D.
"""
function symmetric_double_orbiter(;
    r_pp::Float64 = 4.0,
    R::Float64 = 4.0,
    orb::Float64 = 1.3,
    z_kick::Float64 = 0.0,
    masses::Vector{Float64} = [1.0, 1.0, 1.0, 1.0],
    charges::Vector{Float64} = [1.0, 1.0, -1.0, -1.0],
    dims::Int = 2,
)
    if dims == 2
        X = [
             r_pp/2  0.0;
            -r_pp/2  0.0;
             0.0     R;
             0.0    -R;
        ]
        P = zeros(4, 2)
        v_orb = orb * sqrt(2.0 / R)
        P[3, :] = [masses[3] * v_orb, 0.0]
        P[4, :] = [-masses[4] * v_orb, 0.0]
    else
        X = [
             r_pp/2  0.0  0.0;
            -r_pp/2  0.0  0.0;
             0.0     R    0.0;
             0.0    -R    0.0;
        ]
        P = zeros(4, 3)
        v_orb = orb * sqrt(2.0 / R)
        P[3, :] = [masses[3] * v_orb, 0.0, 0.0]
        P[4, :] = [-masses[4] * v_orb, 0.0, 0.0]
        if z_kick != 0.0
            P[1, 3] = masses[1] * z_kick
            P[2, 3] = -masses[2] * z_kick
        end
    end
    _zero_com!(X, masses)
    _zero_total_p!(P, masses)
    return (_flatten(X), _flatten(P), masses, charges)
end

# ==========================================================================
# 3. Rhombus (from prior study, parameterized for c-variation)
# ==========================================================================

"""
    rhombus_extended(; a, b, q, masses, energy_fraction, velocity_mode)
"""
function rhombus_extended(;
    a::Float64 = 1.5,
    b::Float64 = 1.45,
    masses::Vector{Float64} = [1.0, 1.0, 1.0, 1.0],
    charges::Vector{Float64} = [1.0, 1.0, -1.0, -1.0],
    energy_fraction::Float64 = 0.75,
    velocity_mode::Symbol = :rotating,
)
    X = [
        a  0.0;
       -a  0.0;
        0.0  b;
        0.0 -b;
    ]
    _zero_com!(X, masses)
    U = _coulomb_U(X, charges)
    T_target = energy_fraction * abs(U)
    speed = sqrt(2 * T_target / sum(masses))
    P = zeros(4, 2)
    if velocity_mode == :rotating
        for i = 1:4
            r = X[i, :]
            t = [-r[2], r[1]] / max(norm(r), eps())
            P[i, :] = masses[i] * speed * t
        end
    end
    _zero_total_p!(P, masses)
    return (_flatten(X), _flatten(P), masses, charges)
end

# ==========================================================================
# 4. Unequal mass variants of alternating square
# ==========================================================================

"""
    alternating_square_unequal(; side, masses, charges, energy_fraction, velocity_mode)
"""
function alternating_square_unequal(;
    side::Float64 = 1.5,
    masses::Vector{Float64} = [1.0, 1.0, 1.0, 1.0],
    charges::Vector{Float64} = [1.0, 1.0, -1.0, -1.0],
    energy_fraction::Float64 = 0.5,
    velocity_mode::Symbol = :rotating,
)
    s = side / 2
    # +,-,+,- at vertices: 1(+) at (s,s), 2(+) at (-s,-s), 3(-) at (-s,s), 4(-) at (s,-s)
    X = [
        s  s;
       -s -s;
       -s  s;
        s -s;
    ]
    _zero_com!(X, masses)
    U = _coulomb_U(X, charges)
    T_target = energy_fraction * abs(U)

    P = zeros(4, 2)
    if velocity_mode == :rotating && abs(T_target) > 0
        # Tangential, but mass-weighted
        total_mr2 = sum(masses[i] * dot(X[i,:], X[i,:]) for i in 1:4)
        omega = sqrt(2 * T_target / total_mr2)
        for i in 1:4
            r = X[i, :]
            t = [-r[2], r[1]]
            P[i, :] = masses[i] * omega * t
        end
    elseif velocity_mode == :breathing && abs(T_target) > 0
        total_mr2 = sum(masses[i] * dot(X[i,:], X[i,:]) for i in 1:4)
        omega = sqrt(2 * T_target / total_mr2)
        for i in 1:4
            r = X[i, :]
            P[i, :] = -masses[i] * omega * r
        end
    end
    _zero_total_p!(P, masses)
    return (_flatten(X), _flatten(P), masses, charges)
end
