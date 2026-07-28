# =============================================================================
# Test Problem Builders
# =============================================================================

"""Build a 2D two-body Weber problem."""
function make_weber_problem(;
    dt = 0.001,
    tspan = (0.0, 1.0),
    c = 4.0,
    m1 = 1.0,
    m2 = 0.1,
    k = -0.1,
)
    # Set charges so q1*q2 = k
    q1 = sqrt(abs(k))
    q2 = -sign(k) * sqrt(abs(k))  # q1*q2 = -|k|*sign(k) = k for k < 0

    system = HamiltonianSystem(2, 2)

    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))
    q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p0 = [0.0, m1 * (-m2 / M * v_circ * 0.9), 0.0, m2 * (m1 / M * v_circ * 0.9)]
    HamiltonianProblem(
        system,
        tspan,
        q0,
        p0;
        masses = [m1, m2],
        charges = [q1, q2],
        c = c,
        dt = dt,
    )
end

"""
Build a 2D two-body Coulomb-like problem (Weber with very large c).
The velocity-dependent term becomes negligible when c is large.
"""
function make_coulomb_like_problem(;
    dt = 0.01,
    tspan = (0.0, 1.0),
    m1 = 1.0,
    m2 = 0.5,
    k = 1.0,
)
    # Use very large c to make Weber ≈ Coulomb
    # For attractive Coulomb: U = -k/r, so q1*q2 = -k
    q1 = sqrt(k)
    q2 = -sqrt(k)  # q1*q2 = -k (attractive)

    system = HamiltonianSystem(2, 2)

    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(k * M / (m1 * m2 * r0))
    q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p0 = [0.0, m1 * (-m2 / M * v_circ), 0.0, m2 * (m1 / M * v_circ)]
    HamiltonianProblem(
        system,
        tspan,
        q0,
        p0;
        masses = [m1, m2],
        charges = [q1, q2],
        c = 1e10,
        dt = dt,
    )
end

# =============================================================================
# Test Energy Functions
# =============================================================================

"""
    weber_two_body_velocities_2d(q, p, masses, charges, c) -> (v1, v2, rdot, r)

Physical velocities of a 2-body 2D Weber system from the canonical state.

Independent reference implementation: it uses the **closed-form two-particle
scalar inverse** rather than the package's general `n`-pair linear solve, so
tests built on it check the implementation instead of restating it.

Canonical momentum is `p_i = ∂L/∂v_i`, and for two particles the radial
relation collapses to `p_r = (μ − q₁q₂/(r c²)) ṙ` with `p_r = μ s`,
`s = r̂·(p₁/m₁ − p₂/m₂)`. Hence

```
rdot = μ s / (μ − q₁q₂/(r c²)),   α = q₁q₂ rdot r̂ /(r c²),
v₁   = (p₁ + α)/m₁,               v₂ = (p₂ − α)/m₂.
```
"""
function weber_two_body_velocities_2d(q, p, masses, charges, c)
    m1, m2 = masses
    qc1, qc2 = charges
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    dx, dy = x1 - x2, y1 - y2
    r = sqrt(dx^2 + dy^2)
    rhx, rhy = dx / r, dy / r

    mu = m1 * m2 / (m1 + m2)
    k = qc1 * qc2 / (c^2 * r)
    s = rhx * (px1 / m1 - px2 / m2) + rhy * (py1 / m1 - py2 / m2)
    rdot = mu * s / (mu - k)

    ax = k * rdot * rhx
    ay = k * rdot * rhy

    v1 = ((px1 + ax) / m1, (py1 + ay) / m1)
    v2 = ((px2 - ax) / m2, (py2 - ay) / m2)
    return v1, v2, rdot, r
end

"""
Compute Weber energy for a 2-body 2D system.

The conserved energy is the velocity-space form `E = T_phys + U_Weber`, built
from the physical velocities recovered above. It is *not*
`Σ|p|²/(2m) + U(p/m)`: canonical momentum is not kinetic momentum whenever the
pair has nonzero radial velocity.
"""
function weber_energy_2body_2d(q, p, masses, charges, c)
    m1, m2 = masses
    qc1, qc2 = charges

    v1, v2, rdot, r = weber_two_body_velocities_2d(q, p, masses, charges, c)

    KE = 0.5 * m1 * (v1[1]^2 + v1[2]^2) + 0.5 * m2 * (v2[1]^2 + v2[2]^2)
    PE = qc1 * qc2 / r * (1 - rdot^2 / (2 * c^2))

    return KE + PE
end

"""Compute Coulomb-like energy (Weber with large c, velocity term ignored)."""
function coulomb_like_energy_2body_2d(q, p, masses, charges)
    m1, m2 = masses
    q1, q2 = charges
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    KE = (px1^2 + py1^2) / (2m1) + (px2^2 + py2^2) / (2m2)
    r = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    PE = q1 * q2 / r

    return KE + PE
end
