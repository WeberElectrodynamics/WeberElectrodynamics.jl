"""
Initial condition generators for the 3-body Weber Hamiltonian.
All generators return `(q0, p0, masses, charges)` with flat vectors.

Conventions:
- COM at origin, zero total momentum by construction.
- `energy_fraction` η ∈ (0,1) sets T₀ = η·|U₀| (bound iff η < 1 for attractive configs).
"""

using LinearAlgebra

# --------------------------------------------------------------------------
# Utility
# --------------------------------------------------------------------------
function _coulomb_U_3(X::AbstractMatrix, charges::AbstractVector)
    n = size(X, 1)
    U = 0.0
    @inbounds for i = 1:n, j = (i+1):n
        r = norm(view(X, i, :) .- view(X, j, :))
        U += charges[i] * charges[j] / r
    end
    return U
end

function _zero_com_3!(X::AbstractMatrix, masses::AbstractVector)
    com = (masses' * X) ./ sum(masses)
    X .-= com
    return X
end

function _zero_total_p_3!(P::AbstractMatrix, masses::AbstractVector)
    # Remove total momentum while preserving COM frame
    p_tot = zeros(1, size(P, 2))
    for i = 1:size(P, 1)
        p_tot .+= P[i:i, :]
    end
    # Distribute removal proportional to mass (preserves relative velocities better)
    M = sum(masses)
    for i = 1:size(P, 1)
        P[i:i, :] .-= (masses[i] / M) .* p_tot
    end
    return P
end

function _flatten_3(X::AbstractMatrix)
    n, d = size(X)
    out = Vector{Float64}(undef, n * d)
    @inbounds for i = 1:n, k = 1:d
        out[(i-1)*d+k] = X[i, k]
    end
    return out
end

# --------------------------------------------------------------------------
# Equilateral triangle — rotating or breathing
# charges[1:3] and masses[1:3] are user-specified
# --------------------------------------------------------------------------
function equilateral_triangle(;
    side::Float64 = 2.0,
    charges::Vector{Float64} = [1.0, 1.0, -1.0],
    masses::Vector{Float64} = [1.0, 1.0, 1.0],
    energy_fraction::Float64 = 0.5,
    velocity_mode::Symbol = :rotating,
)
    s = side
    h = s * sqrt(3) / 2
    X = [
        0.0       h*2/3;
        -s/2     -h/3;
         s/2     -h/3;
    ]
    _zero_com_3!(X, masses)
    U = _coulomb_U_3(X, charges)
    T_target = energy_fraction * abs(U)
    # Per-particle KE = T_target/3 = m_i v_i^2 / 2
    P = zeros(3, 2)
    if velocity_mode == :rotating
        for i = 1:3
            r = X[i, :]
            rn = norm(r)
            rn < 1e-14 && continue
            t = [-r[2], r[1]] / rn
            # v_i = speed_i, T = sum m_i v_i^2 /2
            # For equal mass: speed = sqrt(2T/(3m))
            # For unequal: distribute KE equally among particles
            speed_i = sqrt(2 * T_target / (3 * masses[i]))
            P[i, :] = masses[i] * speed_i * t
        end
    elseif velocity_mode == :breathing
        for i = 1:3
            r = X[i, :]
            rn = norm(r)
            rn < 1e-14 && continue
            speed_i = sqrt(2 * T_target / (3 * masses[i]))
            P[i, :] = -masses[i] * speed_i * r / rn
        end
    end
    _zero_total_p_3!(P, masses)
    return (_flatten_3(X), _flatten_3(P), masses, charges)
end

# --------------------------------------------------------------------------
# Collinear configuration along x-axis
# particles at -d, 0, +d (symmetric)
# --------------------------------------------------------------------------
function collinear_symmetric(;
    spacing::Float64 = 2.0,
    charges::Vector{Float64} = [1.0, -1.0, 1.0],
    masses::Vector{Float64} = [1.0, 1.0, 1.0],
    transverse_kick::Float64 = 0.0,
    energy_fraction::Float64 = 0.0,
    velocity_mode::Symbol = :static,
)
    d = spacing
    X = [
        -d  0.0;
         0.0  0.0;
         d  0.0;
    ]
    _zero_com_3!(X, masses)
    U = _coulomb_U_3(X, charges)
    P = zeros(3, 2)
    if velocity_mode == :transverse_stagger
        # Staggered transverse kicks
        T_target = energy_fraction * abs(U)
        speed = sqrt(2 * T_target / (3 * masses[1]))
        P[:, 2] = [speed * masses[1], -speed * masses[2], speed * masses[3]]
    elseif transverse_kick != 0.0
        P[:, 2] = [transverse_kick * masses[1], -transverse_kick * masses[2], transverse_kick * masses[3]]
    end
    _zero_total_p_3!(P, masses)
    return (_flatten_3(X), _flatten_3(P), masses, charges)
end

# --------------------------------------------------------------------------
# Planetary atom: nucleus pair at small separation, orbiter at distance R
# Particles 1,2 form nucleus (like charges), particle 3 is orbiter (opposite)
# --------------------------------------------------------------------------
function planetary_atom_3(;
    nucleus_sep::Float64 = 0.05,
    orbit_radius::Float64 = 1.0,
    charges::Vector{Float64} = [1.0, 1.0, -1.0],
    masses::Vector{Float64} = [1.0, 1.0, 1.0],
    eta::Float64 = 0.7,
    c::Float64 = 4.0,
)
    X = [
        -nucleus_sep/2  0.0;
         nucleus_sep/2  0.0;
         0.0  orbit_radius;
    ]
    _zero_com_3!(X, masses)
    # Compute circular speed for orbiter around effective nucleus (total charge of 1+2)
    Q_eff = abs(charges[1] + charges[2])
    m_nuc = masses[1] + masses[2]
    mu_orb = m_nuc * masses[3] / (m_nuc + masses[3])
    v_circ = sqrt(Q_eff / (mu_orb * orbit_radius))
    v_orb = eta * v_circ
    P = zeros(3, 2)
    # Tangential velocity for orbiter
    P[3, :] = masses[3] * v_orb * [-1.0, 0.0]
    _zero_total_p_3!(P, masses)
    return (_flatten_3(X), _flatten_3(P), masses, charges)
end

# --------------------------------------------------------------------------
# Isosceles configuration: two particles symmetric, one on axis
# --------------------------------------------------------------------------
function isosceles(;
    base::Float64 = 2.0,
    height::Float64 = 1.5,
    charges::Vector{Float64} = [1.0, 1.0, -1.0],
    masses::Vector{Float64} = [1.0, 1.0, 1.0],
    energy_fraction::Float64 = 0.5,
    velocity_mode::Symbol = :rotating,
)
    X = [
        -base/2  0.0;
         base/2  0.0;
         0.0     height;
    ]
    _zero_com_3!(X, masses)
    U = _coulomb_U_3(X, charges)
    T_target = energy_fraction * abs(U)
    P = zeros(3, 2)
    if velocity_mode == :rotating
        for i = 1:3
            r = X[i, :]
            rn = norm(r)
            rn < 1e-14 && continue
            t = [-r[2], r[1]] / rn
            speed_i = sqrt(2 * T_target / (3 * masses[i]))
            P[i, :] = masses[i] * speed_i * t
        end
    end
    _zero_total_p_3!(P, masses)
    return (_flatten_3(X), _flatten_3(P), masses, charges)
end

# --------------------------------------------------------------------------
# Asymmetric: helium-like configurations
# --------------------------------------------------------------------------
function helium_like(;
    orbit_radius::Float64 = 1.5,
    angle_sep::Float64 = Float64(π),   # angle between the two electrons
    charges::Vector{Float64} = [2.0, -1.0, -1.0],
    masses::Vector{Float64} = [10.0, 1.0, 1.0],
    energy_fraction::Float64 = 0.5,
)
    # Nucleus at origin, two orbiters at angle_sep apart
    theta1 = π/2
    theta2 = theta1 + angle_sep
    X = [
        0.0                    0.0;
        orbit_radius*cos(theta1)  orbit_radius*sin(theta1);
        orbit_radius*cos(theta2)  orbit_radius*sin(theta2);
    ]
    _zero_com_3!(X, masses)
    U = _coulomb_U_3(X, charges)
    T_target = energy_fraction * abs(U)
    P = zeros(3, 2)
    # Tangential velocities for orbiters
    for i = 2:3
        r = X[i, :] .- X[1, :]  # relative to nucleus position
        rn = norm(r)
        rn < 1e-14 && continue
        t = [-r[2], r[1]] / rn
        speed_i = sqrt(2 * T_target / (3 * masses[i]))
        P[i, :] = masses[i] * speed_i * t
    end
    _zero_total_p_3!(P, masses)
    return (_flatten_3(X), _flatten_3(P), masses, charges)
end

# --------------------------------------------------------------------------
# Same-sign triangle (all positive or all negative)
# Sub-critical: side < rho = q^2 / (mu * c^2)
# --------------------------------------------------------------------------
function same_sign_triangle(;
    side::Float64 = 0.1,
    charge::Float64 = 1.0,
    mass::Float64 = 1.0,
    c::Float64 = 1.0,
    energy_fraction::Float64 = 0.25,
    velocity_mode::Symbol = :rotating,
)
    charges = fill(charge, 3)
    masses = fill(mass, 3)
    return equilateral_triangle(;
        side=side, charges=charges, masses=masses,
        energy_fraction=energy_fraction, velocity_mode=velocity_mode,
    )
end

# --------------------------------------------------------------------------
# Random perturbation of any IC
# --------------------------------------------------------------------------
function perturb_ic(q0, p0, masses, charges; delta_q=0.01, delta_p=0.01, seed=42)
    rng = Random.MersenneTwister(seed)
    q_pert = q0 .+ delta_q .* randn(rng, length(q0))
    p_pert = p0 .+ delta_p .* randn(rng, length(p0))
    # Re-zero COM and total momentum
    n = length(masses)
    d = length(q0) ÷ n
    X = reshape(copy(q_pert), d, n)'
    P = reshape(copy(p_pert), d, n)'
    _zero_com_3!(X, masses)
    _zero_total_p_3!(P, masses)
    return (_flatten_3(X), _flatten_3(P), masses, charges)
end
