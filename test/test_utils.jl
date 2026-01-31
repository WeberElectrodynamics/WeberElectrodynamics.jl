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

    system = WeberSystem(2, 2)

    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))
    q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p0 = [0.0, m1 * (-m2 / M * v_circ * 0.9), 0.0, m2 * (m1 / M * v_circ * 0.9)]
    WeberProblem(system, tspan, q0, p0; masses = [m1, m2], charges = [q1, q2], c = c, dt = dt)
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

    system = WeberSystem(2, 2)

    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(k * M / (m1 * m2 * r0))
    q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p0 = [0.0, m1 * (-m2 / M * v_circ), 0.0, m2 * (m1 / M * v_circ)]
    WeberProblem(system, tspan, q0, p0; masses = [m1, m2], charges = [q1, q2], c = 1e10, dt = dt)
end

# =============================================================================
# Test Energy Functions
# =============================================================================

"""
Compute Weber energy for a 2-body 2D system.
This is the actual physics energy function for testing.
"""
function weber_energy_2body_2d(q, p, masses, charges, c)
    m1, m2 = masses
    q1, q2 = charges
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    # Kinetic energy
    KE = (px1^2 + py1^2) / (2m1) + (px2^2 + py2^2) / (2m2)

    # Weber potential
    dx, dy = x1 - x2, y1 - y2
    r = sqrt(dx^2 + dy^2)
    vx1, vy1 = px1 / m1, py1 / m1
    vx2, vy2 = px2 / m2, py2 / m2
    rdot = (dx * (vx1 - vx2) + dy * (vy1 - vy2)) / r
    PE = q1 * q2 / r * (1 - rdot^2 / (2 * c^2))

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
