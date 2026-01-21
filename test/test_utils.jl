# =============================================================================
# Test Hamiltonians
# =============================================================================

"""Simple harmonic oscillator: H = p²/2m + kq²/2"""
function harmonic_oscillator_H(q, p, params)
    m, k = params
    sum(p .^ 2) / (2m) + k * sum(q .^ 2) / 2
end

"""Two-body Coulomb (no velocity dependence): H = KE - k/r"""
function coulomb_H(q, p, params)
    m1, m2, k = params
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p
    KE = (px1^2 + py1^2) / (2m1) + (px2^2 + py2^2) / (2m2)
    r = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    KE - k / r
end

"""Two-body Weber (velocity-dependent): H = KE + U_weber"""
function weber_H(q, p, params)
    m1, m2, k, c = params
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p
    KE = (px1^2 + py1^2) / (2m1) + (px2^2 + py2^2) / (2m2)
    dx, dy = x1 - x2, y1 - y2
    r = sqrt(dx^2 + dy^2)
    vx1, vy1 = px1 / m1, py1 / m1
    vx2, vy2 = px2 / m2, py2 / m2
    rdot = (dx * (vx1 - vx2) + dy * (vy1 - vy2)) / r
    PE = k / r * (1 - rdot^2 / (2 * c^2))
    KE + PE
end

# =============================================================================
# Problem Builders
# =============================================================================

"""Build a simple 1D harmonic oscillator problem for fast tests."""
function make_harmonic_problem(; dt=0.01, tspan=(0.0, 1.0), m=1.0, k=1.0, q0=1.0, p0=0.0)
    H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])
    WeberProblem(H, tspan, [q0], [p0]; params=[m, k], dt=dt)
end

"""Build a 2D two-body Coulomb problem."""
function make_coulomb_problem(; dt=0.01, tspan=(0.0, 1.0), m1=1.0, m2=0.5, k=1.0)
    H = compile_hamiltonian(coulomb_H, 2, 2; parameter_names=[:m1, :m2, :k])
    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(k * M / (m1 * m2 * r0))
    q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p0 = [0.0, m1 * (-m2 / M * v_circ), 0.0, m2 * (m1 / M * v_circ)]
    WeberProblem(H, tspan, q0, p0; params=[m1, m2, k], dt=dt)
end

"""Build a 2D two-body Weber problem."""
function make_weber_problem(; dt=0.001, tspan=(0.0, 1.0), c=4.0, m1=1.0, m2=0.1, k=-0.1)
    H = compile_hamiltonian(weber_H, 2, 2; parameter_names=[:m1, :m2, :k, :c])
    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))
    q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p0 = [0.0, m1 * (-m2 / M * v_circ * 0.9), 0.0, m2 * (m1 / M * v_circ * 0.9)]
    WeberProblem(H, tspan, q0, p0; params=[m1, m2, k, c], dt=dt)
end
