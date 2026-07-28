#!/usr/bin/env python3
"""
Symbolic verification of mathematical formulas from:
"Computational Weber electrodynamics: symplectic n-body integration"

Verifies all key derived results using SymPy.
Each check prints PASS/FAIL with the verification method used.

The canonical momentum of the Weber Lagrangian is p_i = dL/dv_i, which is NOT
m_i*v_i: the velocity-dependent pair term contributes a radial correction. Every
group below therefore distinguishes

    s_ij  = rhat_ij . (p_i/m_i - p_j/m_j)     "naive" radial rate (canonical)
    rdot_ij                                    physical radial velocity

which coincide only in the Coulomb limit or at zero radial velocity.

Requires SymPy. Provision it with:

    uv venv .venv && uv pip install --python .venv/bin/python sympy

Run: .venv/bin/python verify_formulas.py
"""

import sys
from random import uniform, seed

import mpmath as mp
from sympy import (
    symbols, sqrt, diff, simplify, expand, trigsimp, cancel,
    Matrix, Rational, Integer, lambdify,
)

seed(42)
mp.mp.dps = 40

# ============================================================
# SECTION 1: Symbol Definitions (shared across all groups)
# ============================================================

# Parameters
c      = symbols('c', positive=True)
q1, q2 = symbols('q1 q2', real=True)
m1, m2 = symbols('m1 m2', positive=True)

# Absolute positions and velocities (2 particles, 2D — main body of paper)
x1, y1, x2, y2     = symbols('x1 y1 x2 y2', real=True)
xd1, yd1, xd2, yd2 = symbols('xd1 yd1 xd2 yd2', real=True)  # physical velocities
px1, py1, px2, py2 = symbols('px1 py1 px2 py2', real=True)  # canonical momenta

# 3D relative variables (for Appendix A.1 radial acceleration identities)
x, y, z           = symbols('x y z', real=True)
xd, yd, zd        = symbols('xd yd zd', real=True)
xdd, ydd, zdd     = symbols('xdd ydd zdd', real=True)

# ---------- Derived 2D expressions ----------
x_rel  = x1 - x2
y_rel  = y1 - y2
xd_rel = xd1 - xd2
yd_rel = yd1 - yd2
r12    = sqrt(x_rel**2 + y_rel**2)

# Physical radial velocity in velocity space (Lagrangian variables)
rdot12_v = (x_rel*xd_rel + y_rel*yd_rel) / r12

# ---------- Potentials (velocity space) ----------
U = q1*q2/r12 * (1 - rdot12_v**2 / (2*c**2))   # Weber potential (Eq. potential)
S = q1*q2/r12 * (1 + rdot12_v**2 / (2*c**2))   # Auxiliary potential

# ---------- Lagrangian (velocity space) ----------
T = m1*(xd1**2 + yd1**2)/2 + m2*(xd2**2 + yd2**2)/2
L = T - S

# ---------- Canonical two-particle machinery, exact ----------
# Reduced mass, pair coupling, naive radial rate, physical radial velocity.
mu    = m1*m2/(m1 + m2)
k12   = q1*q2/(c**2*r12)                       # pair coefficient q1q2/(c^2 r)
s12_p = (x_rel*(px1/m1 - px2/m2) + y_rel*(py1/m1 - py2/m2)) / r12
rdot12_p = mu*s12_p/(mu - k12)                 # exact scalar inverse (Sec. inverse)

# Weber momentum correction alpha, evaluated with the PHYSICAL rdot.
alpha_x = k12*rdot12_p*x_rel/r12
alpha_y = k12*rdot12_p*y_rel/r12

# Physical velocities recovered from the canonical momenta.
vx1_p = (px1 + alpha_x)/m1
vy1_p = (py1 + alpha_y)/m1
vx2_p = (px2 - alpha_x)/m2
vy2_p = (py2 - alpha_y)/m2

# ---------- Exact canonical Hamiltonian in (q,p) space ----------
T_p  = px1**2/(2*m1) + py1**2/(2*m1) + px2**2/(2*m2) + py2**2/(2*m2)
H_qp = T_p + k12*rdot12_p*s12_p/2 + q1*q2/r12

# Forward map p(q, v) from the Lagrangian, used to cross-check H_qp against E.
alpha_x_v = k12*rdot12_v*x_rel/r12
alpha_y_v = k12*rdot12_v*y_rel/r12
FORWARD_MAP = [
    (px1, m1*xd1 - alpha_x_v), (py1, m1*yd1 - alpha_y_v),
    (px2, m2*xd2 + alpha_x_v), (py2, m2*yd2 + alpha_y_v),
]


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def _num_vals_2d():
    """Random 2D (q,p) values with particles separated and away from rho."""
    return {
        x1: uniform(1.0, 3.0),  x2: uniform(-3.0, -1.0),  # |x1-x2| >= 2
        y1: uniform(0.5, 2.0),  y2: uniform(0.5, 2.0),
        px1: uniform(-0.5, 0.5), py1: uniform(-0.5, 0.5),
        px2: uniform(-0.5, 0.5), py2: uniform(-0.5, 0.5),
        m1: uniform(0.5, 2.0),   m2: uniform(0.5, 2.0),
        q1: uniform(0.2, 1.5),   q2: uniform(0.2, 1.5),
        c:  uniform(3.0, 8.0),
    }

def _num_vals_vel():
    """Like _num_vals_2d but with velocity symbols instead of momenta."""
    return {
        x1: uniform(1.0, 3.0),  x2: uniform(-3.0, -1.0),
        y1: uniform(0.5, 2.0),  y2: uniform(0.5, 2.0),
        xd1: uniform(-0.3, 0.3), yd1: uniform(-0.3, 0.3),
        xd2: uniform(-0.3, 0.3), yd2: uniform(-0.3, 0.3),
        m1: uniform(0.5, 2.0),   m2: uniform(0.5, 2.0),
        q1: uniform(0.2, 1.5),   q2: uniform(0.2, 1.5),
        c:  uniform(3.0, 8.0),
    }

def _num_vals_3d():
    """Random 3D relative-coordinate values for appendix checks."""
    return {
        x: uniform(0.5, 2.0), y: uniform(0.5, 2.0), z: uniform(0.5, 2.0),
        xd: uniform(-0.3, 0.3), yd: uniform(-0.3, 0.3), zd: uniform(-0.3, 0.3),
        xdd: uniform(-0.1, 0.1), ydd: uniform(-0.1, 0.1), zdd: uniform(-0.1, 0.1),
    }

def _residuals(expr, val_func, n_trials):
    """Evaluate expr at random points in 40-digit arithmetic; return residuals."""
    probe = val_func()
    order = list(probe.keys())
    f = lambdify(order, expr, modules="mpmath")
    out = []
    for trial in range(n_trials):
        vals = probe if trial == 0 else val_func()
        try:
            out.append(abs(mp.mpf(f(*[mp.mpf(vals[s]) for s in order]))))
        except Exception:
            continue  # skip degenerate sample (e.g. sqrt of negative)
    return out

def numerical_check(expr, val_func, n_trials=25, tol=mp.mpf('1e-25')):
    """Verify expr ~ 0 at random numeric points in high-precision arithmetic."""
    try:
        res = _residuals(expr, val_func, n_trials)
    except Exception:
        return False
    return bool(res) and max(res) <= tol

def check_zero_fast(expr, label, val_func=None, tol=mp.mpf('1e-25')):
    """High-precision numerical check.

    Used where symbolic simplification of the exact canonical Hamiltonian is
    intractable. Residuals below 1e-25 across 25 random points in 40-digit
    arithmetic establish the identity far beyond float round-off.
    """
    vf = val_func or _num_vals_2d
    try:
        res = _residuals(expr, vf, 25)
    except Exception as exc:
        return (False, label, f"FAILED — could not evaluate: {exc}")
    if not res:
        return (False, label, "FAILED — no valid sample points")
    worst = max(res)
    if worst <= tol:
        return (True, label, f"numerically verified (max residual {mp.nstr(worst, 3)})")
    return (False, label, f"FAILED — max residual {mp.nstr(worst, 3)}")

def check_nonzero(expr, label, val_func=None, floor=mp.mpf('1e-6')):
    """Assert expr is genuinely NOT identically zero (used for the two
    'the old formula was wrong' checks)."""
    vf = val_func or _num_vals_2d
    try:
        res = _residuals(expr, vf, 25)
    except Exception as exc:
        return (False, label, f"FAILED — could not evaluate: {exc}")
    if not res:
        return (False, label, "FAILED — no valid sample points")
    best = max(res)
    if best >= floor:
        return (True, label, f"confirmed non-zero (max |residual| {mp.nstr(best, 3)})")
    return (False, label, f"FAILED — residual vanished ({mp.nstr(best, 3)})")

def check_zero(expr, label, val_func=None):
    """Attempt symbolic simplification; fall back to a numerical check."""
    vf = val_func or _num_vals_2d
    res = expr
    try:
        res = cancel(simplify(expand(expr)))
        if res == 0:
            return (True, label, "symbolic zero")
        res2 = trigsimp(res)
        if res2 == 0:
            return (True, label, "symbolic zero (trigsimp)")
    except Exception:
        pass  # proceed to numerical check

    if numerical_check(expr, vf):
        return (True, label, "numerically verified (symbolic simplification inconclusive)")
    return (False, label, f"FAILED — non-zero residual: {res}")

def check_bool(passed, label, detail_pass, detail_fail):
    return (bool(passed), label, detail_pass if passed else detail_fail)


# ============================================================
# GROUP A: Potential Identities (Eqs. potential, S)
# ============================================================

def verify_group_A():
    results = []

    # A.1: Coulomb limit — U -> q1q2/r12 when rdot = 0
    U_static = U.subs([(xd1, 0), (yd1, 0), (xd2, 0), (yd2, 0)])
    results.append(check_zero(
        simplify(U_static) - q1*q2/r12,
        "A.1  Coulomb limit: U(rdot=0) = q1*q2/r",
        _num_vals_vel,
    ))

    # A.2: U + S = 2*q1*q2/r  (the velocity terms cancel by construction)
    results.append(check_zero(
        U + S - 2*q1*q2/r12,
        "A.2  U + S = 2*q1*q2/r",
        _num_vals_vel,
    ))

    # A.3: S - U = q1*q2*rdot^2/(c^2*r)  (sign flip of the velocity term)
    results.append(check_zero(
        S - U - q1*q2*rdot12_v**2/(c**2*r12),
        "A.3  S - U = q1*q2*rdot^2/(c^2*r)",
        _num_vals_vel,
    ))

    return results


# ============================================================
# GROUP B: Euler-Lagrange -> Weber Force (Eq. euler_lagrange)
#
# Uses relative coordinates + chain-rule total time derivative.
# ============================================================

def verify_group_B():
    results = []

    xr, yr     = symbols('xr yr', real=True)
    xdr, ydr   = symbols('xdr ydr', real=True)
    xddr, yddr = symbols('xddr yddr', real=True)

    rr    = sqrt(xr**2 + yr**2)
    rdotr = (xr*xdr + yr*ydr) / rr
    rddotr = ((xr*xddr + yr*yddr)/rr
              + (xdr**2 + ydr**2)/rr
              - (xr*xdr + yr*ydr)**2/rr**3)

    S_rel = q1*q2/rr * (1 + rdotr**2 / (2*c**2))

    dS_dxdr = diff(S_rel, xdr)
    dt_dS_dxdr = (diff(dS_dxdr, xr)*xdr
                  + diff(dS_dxdr, yr)*ydr
                  + diff(dS_dxdr, xdr)*xddr
                  + diff(dS_dxdr, ydr)*yddr)
    dS_dxr = diff(S_rel, xr)

    EL = dt_dS_dxdr - dS_dxr
    F_expected = q1*q2 * xr/rr**3 * (1 - rdotr**2/(2*c**2) + rr*rddotr/c**2)

    def val_func_rel():
        return {
            xr: uniform(0.5, 2.0), yr: uniform(0.5, 2.0),
            xdr: uniform(-0.3, 0.3), ydr: uniform(-0.3, 0.3),
            xddr: uniform(-0.1, 0.1), yddr: uniform(-0.1, 0.1),
            q1: uniform(0.2, 1.5), q2: uniform(0.2, 1.5),
            c: uniform(3.0, 8.0),
        }

    results.append(check_zero(
        EL - F_expected,
        "B.1  Euler-Lagrange: d/dt(dS/dxdot) - dS/dx = Weber force (Eq. euler_lagrange)",
        val_func_rel,
    ))

    return results


# ============================================================
# GROUP C: Canonical momentum and the exact Legendre transform
#          (Eqs. canonical_momentum, H_def, H, h_box, hamiltonian_canonical)
# ============================================================

def verify_group_C():
    results = []

    # C.1: p_x1 = dL/dxd1 = m1*xd1 - alpha_x   (canonical, NOT m1*xd1)
    results.append(check_zero(
        diff(L, xd1) - (m1*xd1 - alpha_x_v),
        "C.1  Canonical momentum: p_x1 = dL/dxd1 = m1*xd1 - alpha_x  [Eq. canonical_momentum]",
        _num_vals_vel,
    ))

    # C.2: p_x2 = dL/dxd2 = m2*xd2 + alpha_x   (opposite correction sign)
    results.append(check_zero(
        diff(L, xd2) - (m2*xd2 + alpha_x_v),
        "C.2  Canonical momentum: p_x2 = dL/dxd2 = m2*xd2 + alpha_x  [Eq. canonical_momentum]",
        _num_vals_vel,
    ))

    # C.3: Canonical momentum is NOT kinetic momentum whenever rdot != 0.
    results.append(check_nonzero(
        diff(L, xd1) - m1*xd1,
        "C.3  p_x1 != m1*xd1 in general (the corrected finding)",
        _num_vals_vel,
    ))

    # C.4: Legendre transform (exact, velocity space): sum(v*dL/dv) - L = T + U
    H_legendre = (xd1*diff(L, xd1) + yd1*diff(L, yd1)
                  + xd2*diff(L, xd2) + yd2*diff(L, yd2)) - L
    results.append(check_zero(
        H_legendre - (T + U),
        "C.4  Legendre transform: E = T + U  (exact, velocity space)  [Eq. h_box]",
        _num_vals_vel,
    ))

    # C.5: Exact canonical Hamiltonian equals the velocity-space energy when the
    #      canonical momenta are supplied by the forward map p(q,v).
    #      This is the central replacement for the old v -> p/m substitution.
    results.append(check_zero_fast(
        H_qp.subs(FORWARD_MAP) - (T + U),
        "C.5  H(q, p(q,v)) = E(q,v)  (exact canonical Hamiltonian)  [Eq. hamiltonian_canonical]",
        _num_vals_vel,
    ))

    # C.6: The naive substitution v -> p/m does NOT give the canonical Hamiltonian.
    H_naive = T_p + q1*q2/r12 * (1 - s12_p**2/(2*c**2))
    results.append(check_nonzero(
        H_qp - H_naive,
        "C.6  H(q,p) != naive T(p) + U(s)  (the old formula is not the Legendre transform)",
        _num_vals_2d,
        floor=mp.mpf('1e-12'),
    ))

    # C.7: Two-particle scalar inverse  p_r = (mu - q1q2/(r c^2)) * rdot.
    #      Built from the velocity side so it is not true by construction.
    results.append(check_zero_fast(
        (mu*s12_p).subs(FORWARD_MAP) - (mu - k12)*rdot12_v,
        "C.7  p_r = (mu - q1q2/(r c^2)) * rdot  (two-particle scalar inverse)",
        _num_vals_vel,
    ))

    # C.8: Round trip — rdot recovered from (q,p) equals the physical rdot.
    results.append(check_zero_fast(
        rdot12_p.subs(FORWARD_MAP) - rdot12_v,
        "C.8  rdot(q,p) round-trips to the physical rdot",
        _num_vals_vel,
    ))

    # C.9: Coulomb limit — as c grows the exact H approaches T(p) + q1q2/r.
    def _num_vals_bigc():
        vals = _num_vals_2d()
        vals[c] = 1.0e9
        return vals
    results.append(check_zero_fast(
        H_qp - (T_p + q1*q2/r12),
        "C.9  Coulomb limit: H -> T(p) + q1*q2/r as c -> infinity",
        _num_vals_bigc,
        tol=mp.mpf('1e-15'),
    ))

    # C.10: Zero relative velocity => zero radial velocity => H = T(p) + q1q2/r.
    zero_rel = [(xd2, xd1), (yd2, yd1)]
    results.append(check_zero_fast(
        ((H_qp - (T_p + q1*q2/r12)).subs(FORWARD_MAP)).subs(zero_rel),
        "C.10 Zero relative velocity: H = T(p) + q1*q2/r",
        _num_vals_vel,
    ))

    return results


# ============================================================
# GROUP D: 2-Particle Hamilton's Equations
#          (Eqs. xdot_two, xdot_two_p2, pdot_two, pdot_two_p2)
# ============================================================

def verify_group_D():
    results = []

    # D.1: dH/dpx1 = (px1 + alpha_x)/m1   [Eq. xdot_two, corrected sign]
    results.append(check_zero_fast(
        diff(H_qp, px1) - vx1_p,
        "D.1  xdot1 = dH/dpx1 = (px1 + alpha_x)/m1  [Eq. xdot_two]",
    ))

    # D.2: dH/dpx2 = (px2 - alpha_x)/m2   [Eq. xdot_two_p2, corrected sign]
    results.append(check_zero_fast(
        diff(H_qp, px2) - vx2_p,
        "D.2  xdot2 = dH/dpx2 = (px2 - alpha_x)/m2  [Eq. xdot_two_p2]",
    ))

    # D.3: -dH/dx1 = q1q2/r^2 * [(x1-x2)/r*(1 + 3 rdot^2/(2c^2)) - rdot*(vx1-vx2)/c^2]
    #      [Eq. pdot_two, both 1/c^2 signs corrected]
    D3_expected = (q1*q2/r12**2
                   * (x_rel/r12 * (1 + Rational(3, 2)*rdot12_p**2/c**2)
                      - rdot12_p*(vx1_p - vx2_p)/c**2))
    results.append(check_zero_fast(
        -diff(H_qp, x1) - D3_expected,
        "D.3  pxdot1 = -dH/dx1 matches the corrected formula  [Eq. pdot_two]",
    ))

    # D.4: pxdot1 + pxdot2 = 0  [Eq. pdot_two_p2]
    results.append(check_zero_fast(
        -diff(H_qp, x1) - diff(H_qp, x2),
        "D.4  pxdot1 + pxdot2 = 0  (Newton's 3rd law, Eq. pdot_two_p2)",
    ))

    # D.5: the y-direction analogues hold with the same structure.
    D5_expected = (q1*q2/r12**2
                   * (y_rel/r12 * (1 + Rational(3, 2)*rdot12_p**2/c**2)
                      - rdot12_p*(vy1_p - vy2_p)/c**2))
    results.append(check_zero_fast(
        -diff(H_qp, y1) - D5_expected,
        "D.5  pydot1 = -dH/dy1 matches the corrected formula (y analogue)",
    ))

    return results


# ============================================================
# GROUP E: n-Particle EOM Consistency (n=2 explicit summation)
#          Verifies Eqs. xdot_expanded and pdot_expanded reduce to the
#          2-particle formulas when written with explicit sums.
# ============================================================

def verify_group_E():
    results = []

    # E.1a: n=2, i=1 — xdot_expanded carries a PLUS on the summed correction.
    xdot1_npart = (px1 + q1*q2/c**2 * rdot12_p*x_rel/r12**2) / m1
    results.append(check_zero_fast(
        xdot1_npart - diff(H_qp, px1),
        "E.1a xdot_expanded (n=2, i=1) equals dH/dpx1  [Eq. xdot_expanded]",
    ))

    # E.1b: n=2, i=2 — (x2-x1) = -(x1-x2) flips the correction back.
    xdot2_npart = (px2 + q2*q1/c**2 * rdot12_p*(x2 - x1)/r12**2) / m2
    results.append(check_zero_fast(
        xdot2_npart - diff(H_qp, px2),
        "E.1b xdot_expanded (n=2, i=2) equals dH/dpx2  [Eq. xdot_expanded]",
    ))

    # E.2: n=2, i=1 — pdot_expanded with both corrected 1/c^2 signs.
    pdot1_npart = (q1*q2/r12**2
                   * (x_rel/r12 * (1 + Rational(3, 2)*rdot12_p**2/c**2)
                      - rdot12_p*(vx1_p - vx2_p)/c**2))
    results.append(check_zero_fast(
        pdot1_npart - (-diff(H_qp, x1)),
        "E.2  pdot_expanded (n=2, i=1) equals -dH/dx1  [Eq. pdot_expanded]",
    ))

    return results


# ============================================================
# GROUP F: Radial Acceleration Identities (Appendix A.1)
# ============================================================

def verify_group_F():
    results = []

    r3d    = sqrt(x**2 + y**2 + z**2)
    rdot3d = (x*xd + y*yd + z*zd) / r3d

    rddot_F1 = ((x*xdd + y*ydd + z*zdd)/r3d
                + (xd**2 + yd**2 + zd**2)/r3d
                - (x*xd + y*yd + z*zd)**2/r3d**3)

    rddot_F3 = (1/r3d)*((x*xdd + y*ydd + z*zdd)
                        + (xd**2 + yd**2 + zd**2)
                        - rdot3d**2)

    results.append(check_zero(
        rddot_F1 - rddot_F3,
        "F.1  rddot Form1 = Form3  (three-form equivalence, Appendix A.1)",
        _num_vals_3d,
    ))

    def total_dt(f, pos, vel, acc):
        return (sum(diff(f, p)*v for p, v in zip(pos, vel))
                + sum(diff(f, v)*a for v, a in zip(vel, acc)))

    pos_syms = [x, y, z]
    vel_syms = [xd, yd, zd]
    acc_syms = [xdd, ydd, zdd]

    rdot_expr = total_dt(r3d, pos_syms, vel_syms, acc_syms)
    results.append(check_zero(
        rdot_expr - rdot3d,
        "F.2a rdot = d/dt(|r|) derived equals (r.v)/|r|",
        _num_vals_3d,
    ))

    rddot_derived = total_dt(rdot_expr, pos_syms, vel_syms, acc_syms)
    results.append(check_zero(
        rddot_derived - rddot_F1,
        "F.2b rddot Form1 = d2/dt2(|r|) derived from first principles",
        _num_vals_3d,
    ))

    return results


# ============================================================
# GROUP G: Conservation Laws (Section 2.2), corrected Hamiltonian
# ============================================================

def verify_group_G():
    results = []
    t_sym = symbols('t')

    results.append(check_bool(
        t_sym not in H_qp.free_symbols,
        "G.1  H(q,p) has no explicit time dependence",
        "PASS: t not in free_symbols",
        "FAIL: t found in H",
    ))

    pairs = [(x1, px1), (y1, py1), (x2, px2), (y2, py2)]
    dH_dt = sum(
        diff(H_qp, qi)*diff(H_qp, pi) - diff(H_qp, pi)*diff(H_qp, qi)
        for qi, pi in pairs
    )
    results.append(check_bool(
        dH_dt == 0,
        "G.2  dH/dt = 0 via Hamilton's equations (Poisson bracket cancellation)",
        "trivially zero (a*b - b*a = 0)",
        f"non-zero: {dH_dt}",
    ))

    results.append(check_zero_fast(
        -diff(H_qp, x1) - diff(H_qp, x2),
        "G.3  Total x-momentum rate = -dH/dx1 - dH/dx2 = 0",
    ))

    results.append(check_zero_fast(
        -diff(H_qp, y1) - diff(H_qp, y2),
        "G.4  Total y-momentum rate = -dH/dy1 - dH/dy2 = 0",
    ))

    return results


# ============================================================
# GROUP H: n=3 simultaneous velocity solve and canonical equations
#
# The exact canonical Hamiltonian requires solving a coupled linear system for
# the physical pair radial velocities, so for n=3 it is evaluated numerically at
# high precision (mpmath, 40 digits) and its gradients are checked by central
# differences. This is an independent check: the derivative formulas under test
# never enter the evaluation of H.
# ============================================================

PAIRS3 = [(0, 1), (0, 2), (1, 2)]


def _weber_pair_solve(qs, ps, ms, chs, cval, n, d):
    """Physical pair radial velocities rdot_a and particle velocities v_i.

    Solves (I - G K) rdot = s over the n(n-1)/2 pairs, where
      s_a   = rhat_a . (p_i/m_i - p_j/m_j)
      k_a   = q_i q_j / (c^2 r_a)
      G_ab  = (rhat_a . rhat_b) * (di_k/m_i - di_l/m_i - dj_k/m_j + dj_l/m_j)
    """
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    P = len(pairs)

    rhat, rmag, kco, s = [], [], [], []
    for (i, j) in pairs:
        dq = [qs[i*d + a] - qs[j*d + a] for a in range(d)]
        r = mp.sqrt(sum(t*t for t in dq))
        rh = [t/r for t in dq]
        rhat.append(rh)
        rmag.append(r)
        kco.append(chs[i]*chs[j]/(cval**2 * r))
        s.append(sum(rh[a]*(ps[i*d + a]/ms[i] - ps[j*d + a]/ms[j]) for a in range(d)))

    A = mp.zeros(P, P)
    for b in range(P):
        for a in range(P):
            ib, jb = pairs[b]
            ia, ja = pairs[a]
            dot = sum(rhat[b][t]*rhat[a][t] for t in range(d))
            coup = ((1/ms[ib] if ib == ia else 0) - (1/ms[ib] if ib == ja else 0)
                    - (1/ms[jb] if jb == ia else 0) + (1/ms[jb] if jb == ja else 0))
            A[b, a] = (1 if a == b else 0) - dot*coup*kco[a]

    rdot = mp.lu_solve(A, mp.matrix(s))

    v = [ps[i*d + a]/ms[i] for i in range(n) for a in range(d)]
    for a, (i, j) in enumerate(pairs):
        for t in range(d):
            v[i*d + t] += kco[a]*rdot[a]*rhat[a][t]/ms[i]
            v[j*d + t] -= kco[a]*rdot[a]*rhat[a][t]/ms[j]
    return pairs, rmag, rhat, kco, s, [rdot[a] for a in range(P)], v


def _weber_H(qs, ps, ms, chs, cval, n, d):
    """Exact canonical Weber Hamiltonian, evaluated numerically."""
    pairs, rmag, _, kco, s, rdot, _ = _weber_pair_solve(qs, ps, ms, chs, cval, n, d)
    H = sum(sum(ps[i*d + a]**2 for a in range(d))/(2*ms[i]) for i in range(n))
    for a, (i, j) in enumerate(pairs):
        H += kco[a]*rdot[a]*s[a]/2 + chs[i]*chs[j]/rmag[a]
    return H


def _weber_pdot(qs, ps, ms, chs, cval, n, d):
    """Corrected canonical momentum rates from the closed-form expression."""
    pairs, rmag, rhat, kco, _, rdot, v = _weber_pair_solve(qs, ps, ms, chs, cval, n, d)
    out = [mp.mpf(0)]*(n*d)
    for a, (i, j) in enumerate(pairs):
        pref = chs[i]*chs[j]/rmag[a]**2
        for t in range(d):
            dv = v[i*d + t] - v[j*d + t]
            term = pref*(rhat[a][t]*(1 + 3*rdot[a]**2/(2*cval**2))
                         - rdot[a]*dv/cval**2)
            out[i*d + t] += term
            out[j*d + t] -= term
    return out


def _weber_forward_p(qs, vs, ms, chs, cval, n, d):
    """Forward map p(q,v) from the Lagrangian: p_i = m_i v_i - sum_j alpha_ij."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    p = [ms[i]*vs[i*d + a] for i in range(n) for a in range(d)]
    for (i, j) in pairs:
        dq = [qs[i*d + a] - qs[j*d + a] for a in range(d)]
        r = mp.sqrt(sum(t*t for t in dq))
        rh = [t/r for t in dq]
        rd = sum(rh[a]*(vs[i*d + a] - vs[j*d + a]) for a in range(d))
        k = chs[i]*chs[j]/(cval**2 * r)
        for t in range(d):
            p[i*d + t] -= k*rd*rh[t]
            p[j*d + t] += k*rd*rh[t]
    return p


def _weber_energy_velocity(qs, vs, ms, chs, cval, n, d):
    """Velocity-space energy E = sum m v^2/2 + sum q_i q_j/r (1 - rdot^2/(2c^2))."""
    E = sum(ms[i]*sum(vs[i*d + a]**2 for a in range(d))/2 for i in range(n))
    for i in range(n):
        for j in range(i + 1, n):
            dq = [qs[i*d + a] - qs[j*d + a] for a in range(d)]
            r = mp.sqrt(sum(t*t for t in dq))
            rd = sum(dq[a]*(vs[i*d + a] - vs[j*d + a]) for a in range(d))/r
            E += chs[i]*chs[j]/r*(1 - rd**2/(2*cval**2))
    return E


def _sample_3body():
    qs = [mp.mpf(v) for v in (4.5, 2.4, 1.0, -1.0, -2.0, 2.6)]
    ps = [mp.mpf(v) for v in (0.31, -0.22, -0.17, 0.28, -0.14, -0.06)]
    ms = [mp.mpf(v) for v in (1.3, 0.7, 1.9)]
    chs = [mp.mpf(v) for v in (0.9, -1.1, 0.6)]
    return qs, ps, ms, chs, mp.mpf('5.0'), 3, 2


def _central_diff(f, vec, idx, h=mp.mpf('1e-12')):
    plus = list(vec); plus[idx] = plus[idx] + h
    minus = list(vec); minus[idx] = minus[idx] - h
    return (f(plus) - f(minus))/(2*h)


def verify_group_H():
    results = []
    tol = mp.mpf('1e-16')
    qs, ps, ms, chs, cval, n, d = _sample_3body()

    # H.1: the pair solve round-trips against the Lagrangian forward map.
    _, _, _, _, _, _, v = _weber_pair_solve(qs, ps, ms, chs, cval, n, d)
    p_round = _weber_forward_p(qs, v, ms, chs, cval, n, d)
    err = max(abs(p_round[i] - ps[i]) for i in range(n*d))
    results.append(check_bool(
        err < mp.mpf('1e-30'),
        "H.1  n=3 simultaneous velocity solve round-trips through p(q,v)",
        f"max |p_roundtrip - p| = {mp.nstr(err, 3)}",
        f"FAILED — max error {mp.nstr(err, 3)}",
    ))

    # H.2: dH/dp_i equals the recovered physical velocity (first canonical equation).
    fH_p = lambda pp: _weber_H(qs, pp, ms, chs, cval, n, d)
    err = max(abs(_central_diff(fH_p, ps, i) - v[i]) for i in range(n*d))
    results.append(check_bool(
        err < tol,
        "H.2  n=3: dH/dp_i = v_i  (coordinate rate is the physical velocity)",
        f"max |dH/dp - v| = {mp.nstr(err, 3)}",
        f"FAILED — max error {mp.nstr(err, 3)}",
    ))

    # H.3: -dH/dq_i equals the corrected canonical momentum rate.
    fH_q = lambda qq: _weber_H(qq, ps, ms, chs, cval, n, d)
    pdot = _weber_pdot(qs, ps, ms, chs, cval, n, d)
    err = max(abs(-_central_diff(fH_q, qs, i) - pdot[i]) for i in range(n*d))
    results.append(check_bool(
        err < tol,
        "H.3  n=3: -dH/dq_i equals the corrected pdot formula  [Eq. pdot_expanded]",
        f"max |-dH/dq - pdot| = {mp.nstr(err, 3)}",
        f"FAILED — max error {mp.nstr(err, 3)}",
    ))

    # H.4: total momentum rate vanishes for every component.
    err = max(abs(sum(pdot[i*d + t] for i in range(n))) for t in range(d))
    results.append(check_bool(
        err < mp.mpf('1e-30'),
        "H.4  n=3: sum_i pdot_i = 0  (linear momentum conservation)",
        f"max |sum pdot| = {mp.nstr(err, 3)}",
        f"FAILED — max error {mp.nstr(err, 3)}",
    ))

    # H.5: exact Legendre transform for n=3.
    vs = [mp.mpf(t) for t in ('0.21', '-0.13', '-0.09', '0.17', '0.05', '-0.11')]
    p_fw = _weber_forward_p(qs, vs, ms, chs, cval, n, d)
    err = abs(_weber_H(qs, p_fw, ms, chs, cval, n, d)
              - _weber_energy_velocity(qs, vs, ms, chs, cval, n, d))
    results.append(check_bool(
        err < mp.mpf('1e-30'),
        "H.5  n=3: H(q, p(q,v)) = E(q,v)  (exact Legendre transform)",
        f"|H - E| = {mp.nstr(err, 3)}",
        f"FAILED — residual {mp.nstr(err, 3)}",
    ))

    # H.6: Coulomb limit — H -> sum p^2/2m + sum q_i q_j / r as c -> oo.
    H_big_c = _weber_H(qs, ps, ms, chs, mp.mpf('1e12'), n, d)
    H_coulomb = sum(sum(ps[i*d + a]**2 for a in range(d))/(2*ms[i]) for i in range(n))
    for i in range(n):
        for j in range(i + 1, n):
            r = mp.sqrt(sum((qs[i*d + a] - qs[j*d + a])**2 for a in range(d)))
            H_coulomb += chs[i]*chs[j]/r
    results.append(check_bool(
        abs(H_big_c - H_coulomb) < mp.mpf('1e-20'),
        "H.6  n=3 Coulomb limit: H -> sum p^2/(2m) + sum q_i q_j / r",
        f"|H(c=1e12) - H_Coulomb| = {mp.nstr(abs(H_big_c - H_coulomb), 3)}",
        "FAILED — Coulomb limit not recovered",
    ))

    # H.7: zero radial velocity — H agrees with the Coulomb form at that instant.
    v_rigid = [mp.mpf('0.2'), mp.mpf('0.1')]*n   # identical velocities => all rdot = 0
    p_rigid = _weber_forward_p(qs, v_rigid, ms, chs, cval, n, d)
    H_rigid = _weber_H(qs, p_rigid, ms, chs, cval, n, d)
    H_ref = sum(sum(p_rigid[i*d + a]**2 for a in range(d))/(2*ms[i]) for i in range(n))
    for i in range(n):
        for j in range(i + 1, n):
            r = mp.sqrt(sum((qs[i*d + a] - qs[j*d + a])**2 for a in range(d)))
            H_ref += chs[i]*chs[j]/r
    results.append(check_bool(
        abs(H_rigid - H_ref) < mp.mpf('1e-30'),
        "H.7  n=3 zero-radial-velocity limit: H = sum p^2/(2m) + sum q_i q_j / r",
        f"|H - H_Coulomb| = {mp.nstr(abs(H_rigid - H_ref), 3)}",
        "FAILED — zero-radial limit not recovered",
    ))

    return results


# ============================================================
# GROUP I: Angular momentum conservation (corrected Hamiltonian)
# ============================================================

def verify_group_I():
    results = []

    def poisson_bracket_Lz(H_expr, coord_tuples):
        bracket = Integer(0)
        for xi, yi, pxi, pyi in coord_tuples:
            bracket += (diff(H_expr, xi)*(-yi) - diff(H_expr, pxi)*pyi
                        + diff(H_expr, yi)*xi - diff(H_expr, pyi)*(-pxi))
        return bracket

    results.append(check_zero_fast(
        poisson_bracket_Lz(H_qp, [(x1, y1, px1, py1), (x2, y2, px2, py2)]),
        "I.1  {H(q,p), Lz} = 0  (angular momentum conservation, n=2)",
        _num_vals_2d,
    ))

    # I.2: n=3 rotational invariance of the exact canonical Hamiltonian,
    #      checked numerically by rotating the whole phase-space point.
    qs, ps, ms, chs, cval, n, d = _sample_3body()
    H0 = _weber_H(qs, ps, ms, chs, cval, n, d)
    th = mp.mpf('0.7')
    ct, st = mp.cos(th), mp.sin(th)
    qr, pr = [], []
    for i in range(n):
        qr += [ct*qs[i*d] - st*qs[i*d+1], st*qs[i*d] + ct*qs[i*d+1]]
        pr += [ct*ps[i*d] - st*ps[i*d+1], st*ps[i*d] + ct*ps[i*d+1]]
    err = abs(_weber_H(qr, pr, ms, chs, cval, n, d) - H0)
    results.append(check_bool(
        err < mp.mpf('1e-30'),
        "I.2  n=3: H is invariant under rigid rotation of (q,p)",
        f"|H(rotated) - H| = {mp.nstr(err, 3)}",
        f"FAILED — residual {mp.nstr(err, 3)}",
    ))

    # I.3: n=3 translation invariance.
    shift = [mp.mpf('1.7'), mp.mpf('-0.9')]
    qt = [qs[i*d + t] + shift[t] for i in range(n) for t in range(d)]
    err = abs(_weber_H(qt, ps, ms, chs, cval, n, d) - H0)
    results.append(check_bool(
        err < mp.mpf('1e-30'),
        "I.3  n=3: H is invariant under rigid translation of q",
        f"|H(shifted) - H| = {mp.nstr(err, 3)}",
        f"FAILED — residual {mp.nstr(err, 3)}",
    ))

    return results


# ============================================================
# GROUP J: Canonical equations reproduce the mechanical Weber force
#
# m_i a_i = pdot_i + d(alpha_i)/dt. Substituting the corrected canonical
# momentum rate must return Weber's force law exactly.
# ============================================================

def verify_group_J():
    results = []

    xr, yr     = symbols('xr yr', real=True)
    xdr, ydr   = symbols('xdr ydr', real=True)
    xddr, yddr = symbols('xddr yddr', real=True)

    rr     = sqrt(xr**2 + yr**2)
    rdotr  = (xr*xdr + yr*ydr)/rr
    rddotr = ((xr*xddr + yr*yddr)/rr + (xdr**2 + ydr**2)/rr
              - (xr*xdr + yr*ydr)**2/rr**3)

    # alpha_x for the pair, in relative coordinates, with the physical rdot.
    alpha = q1*q2/c**2 * rdotr*xr/rr**2

    # Total time derivative of alpha along the trajectory.
    alpha_dot = (diff(alpha, xr)*xdr + diff(alpha, yr)*ydr
                 + diff(alpha, xdr)*xddr + diff(alpha, ydr)*yddr)

    # Corrected canonical momentum rate for particle 1, x component.
    pdot = q1*q2/rr**2*(xr/rr*(1 + Rational(3, 2)*rdotr**2/c**2)
                        - rdotr*xdr/c**2)

    # Weber force, x component.
    F_weber = q1*q2*xr/rr**3*(1 - rdotr**2/(2*c**2) + rr*rddotr/c**2)

    def val_func_rel():
        return {
            xr: uniform(0.5, 2.0), yr: uniform(0.5, 2.0),
            xdr: uniform(-0.3, 0.3), ydr: uniform(-0.3, 0.3),
            xddr: uniform(-0.1, 0.1), yddr: uniform(-0.1, 0.1),
            q1: uniform(0.2, 1.5), q2: uniform(0.2, 1.5),
            c: uniform(3.0, 8.0),
        }

    results.append(check_zero(
        (pdot + alpha_dot) - F_weber,
        "J.1  m*a = pdot + d(alpha)/dt reproduces Weber's force law",
        val_func_rel,
    ))

    return results


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 68)
    print("  Weber Electrodynamics — Symbolic Formula Verification (SymPy)")
    print("=" * 68)

    all_results = []

    print("\n--- Group A: Potential Identities ---", flush=True)
    all_results += verify_group_A()

    print("--- Group B: Euler-Lagrange Equation ---", flush=True)
    all_results += verify_group_B()

    print("--- Group C: Canonical Momentum and Legendre Transform ---", flush=True)
    all_results += verify_group_C()

    print("--- Group D: 2-Particle Hamilton's Equations ---", flush=True)
    all_results += verify_group_D()

    print("--- Group E: n-Particle EOM Consistency (n=2) ---", flush=True)
    all_results += verify_group_E()

    print("--- Group F: Radial Acceleration Identities ---", flush=True)
    all_results += verify_group_F()

    print("--- Group G: Conservation Laws ---", flush=True)
    all_results += verify_group_G()

    print("--- Group H: n=3 Simultaneous Velocity Solve ---", flush=True)
    all_results += verify_group_H()

    print("--- Group I: Rotational and Translational Invariance ---", flush=True)
    all_results += verify_group_I()

    print("--- Group J: Recovery of the Weber Force ---", flush=True)
    all_results += verify_group_J()

    print("\n" + "=" * 68)
    print("  VERIFICATION SUMMARY")
    print("=" * 68)
    n_pass = n_fail = 0
    for passed, label, detail in all_results:
        if passed:
            n_pass += 1
            mark = "PASS ✓"
        else:
            n_fail += 1
            mark = "FAIL ✗"
        print(f"  [{mark}]  {label}")
        if not passed or "numerically" in detail:
            print(f"            → {detail}")

    print()
    print(f"  Result: {n_pass}/{n_pass + n_fail} checks passed.", end="")
    if n_fail == 0:
        print("  All formulas verified. ✓")
    else:
        print(f"  {n_fail} check(s) FAILED.")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
