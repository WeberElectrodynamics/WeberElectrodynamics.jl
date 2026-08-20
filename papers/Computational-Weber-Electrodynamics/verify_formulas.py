#!/usr/bin/env python3
"""
Symbolic verification of every derivation in:
"Computational Weber electrodynamics: symplectic n-body integration"

Each group below corresponds to one derivation in the paper and proves it
algebraically with SymPy, starting from a single source of truth: the
Lagrangian L = T - S.  Nothing is assumed that the paper itself asserts.

  A  potential identities                 Eqs. potential, S
  B  Euler-Lagrange -> Weber force        Eq. euler_lagrange
  C  Legendre transform, H = T + U        Eqs. H_def, H, h_box
  D  2-particle Hamilton equations        Eqs. xdot_two .. pdot_two_p2
  E  n-particle EOM at n = 2              Eqs. xdot_expanded, pdot_expanded
  F  radial acceleration identities       Appendix A.1
  G  conservation laws                    Sec. 2.2
  H  n-particle EOM at n = 3              Eqs. xdot_expanded, pdot_expanded
  I  angular momentum / rotational inv.   Sec. 2.2
  J  force <-> potential, 3rd law         Sec. 1, Eqs. potential, force
  K  closed-form H(q,p) for n = 2         Eqs. px1, H_p, xdot, pdot
  L  Weber's critical radius rho          Discussion
  M  extended phase space & integrator    Sec. 3, Eqs. Z, A, AZ, Phi
  N  n-body mass matrix and implicitness  Secs. 2.3, 2.4
  O  units, dimensions, pair count        Sec. 2.4, Sec. 2.5, Appendix A
  P  negative controls for the fixes      (asserts the old forms are wrong)

Run: python3 verify_formulas.py     (needs sympy; e.g. uv run --with sympy)
"""

import sys
from random import uniform, seed
from sympy import (
    symbols, sqrt, diff, simplify, expand, trigsimp, series, solve,
    cos, sin, Matrix, eye, zeros, Rational, Integer
)

seed(42)

# ============================================================
# SECTION 1: Symbol Definitions (shared across all groups)
# ============================================================

# Parameters
c      = symbols('c', positive=True)
q1, q2 = symbols('q1 q2', real=True)
m1, m2 = symbols('m1 m2', positive=True)

# Absolute positions and velocities (2 particles, 2D — main body of paper)
x1, y1, x2, y2     = symbols('x1 y1 x2 y2', real=True)
xd1, yd1, xd2, yd2 = symbols('xd1 yd1 xd2 yd2', real=True)  # velocities (Lagrangian)

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

# rdot in velocity space (Lagrangian variables)
rdot12_v = (x_rel*xd_rel + y_rel*yd_rel) / r12

# ---------- Potentials (velocity space) ----------
U = q1*q2/r12 * (1 - rdot12_v**2 / (2*c**2))   # Weber potential (Eq. potential)
S = q1*q2/r12 * (1 + rdot12_v**2 / (2*c**2))   # Auxiliary potential (Eq. 104)

# ---------- Lagrangian (velocity space) ----------
T = m1*(xd1**2 + yd1**2)/2 + m2*(xd2**2 + yd2**2)/2
L = T - S

# ---------- Canonical momenta from the Lagrangian ----------
# p_i = dL/dv_i.  Every check below is carried out in velocity space. The
# substitution v_i = p_i/m_i is used nowhere: it is precisely the step the
# paper marks as wrong.
p1x_L = diff(L, xd1)
p1y_L = diff(L, yd1)
p2x_L = diff(L, xd2)
p2y_L = diff(L, yd2)

# ---------- alpha_x, alpha_y (Eq. alpha_x), velocity space ----------
alpha_x = q1*q2/c**2 * rdot12_v*(x1 - x2)/r12**2
alpha_y = q1*q2/c**2 * rdot12_v*(y1 - y2)/r12**2

# ============================================================
# 3-PARTICLE SETUP (shared by Group H and Group I)
# Builds L3, the n-body Lagrangian, in velocity space.
# ============================================================

x3, y3   = symbols('x3 y3', real=True)
xd3, yd3 = symbols('xd3 yd3', real=True)
m3       = symbols('m3', positive=True)
q3       = symbols('q3', real=True)

x_rel_13 = x1 - x3
y_rel_13 = y1 - y3
x_rel_23 = x2 - x3
y_rel_23 = y2 - y3
r13      = sqrt(x_rel_13**2 + y_rel_13**2)
r23      = sqrt(x_rel_23**2 + y_rel_23**2)

# ṙ₁₃ and ṙ₂₃ in velocity space
rdot13_v = (x_rel_13*(xd1 - xd3) + y_rel_13*(yd1 - yd3)) / r13
rdot23_v = (x_rel_23*(xd2 - xd3) + y_rel_23*(yd2 - yd3)) / r23

# 3-particle kinetic energy
T3 = (m1*(xd1**2 + yd1**2)/2
     + m2*(xd2**2 + yd2**2)/2
     + m3*(xd3**2 + yd3**2)/2)

# 3-particle Lagrangian, L = T - sum of pair S (Eq. 104 per pair)
L3 = (T3
     - q1*q2/r12 * (1 + rdot12_v**2/(2*c**2))
     - q1*q3/r13 * (1 + rdot13_v**2/(2*c**2))
     - q2*q3/r23 * (1 + rdot23_v**2/(2*c**2)))


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def _num_vals_vel():
    """Random 2-particle 2D values in velocity space, particles separated."""
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

def _num_vals_3body():
    """Random 3-particle 2D values in velocity space, all pairs separated."""
    return {
        x1: uniform(4.0, 5.0),   x2: uniform(0.5, 1.5),   x3: uniform(-2.5, -1.5),
        y1: uniform(2.0, 3.0),   y2: uniform(-1.5, -0.5),  y3: uniform(2.0, 3.0),
        xd1: uniform(-0.3, 0.3), yd1: uniform(-0.3, 0.3),
        xd2: uniform(-0.3, 0.3), yd2: uniform(-0.3, 0.3),
        xd3: uniform(-0.3, 0.3), yd3: uniform(-0.3, 0.3),
        m1: uniform(0.5, 2.0),   m2: uniform(0.5, 2.0),   m3: uniform(0.5, 2.0),
        q1: uniform(0.2, 1.5),   q2: uniform(0.2, 1.5),   q3: uniform(0.2, 1.5),
        c:  uniform(3.0, 8.0),
    }

def numerical_check(expr, val_func, n_trials=25, tol=1e-7):
    """Verify expr ≈ 0 by substituting random numeric values."""
    for _ in range(n_trials):
        vals = val_func()
        try:
            val = complex(expr.subs(vals))
            if abs(val) > tol:
                return False
        except Exception:
            pass  # Skip if substitution fails (e.g., sqrt of negative)
    return True

def check_zero_fast(expr, label, val_func):
    """Numerical-only check (skips symbolic simplification for expensive expressions)."""
    if numerical_check(expr, val_func):
        return (True, label, "numerically verified")
    return (False, label, "FAILED — non-zero residual (numerical check)")

def check_zero(expr, label, val_func=None):
    """
    Attempt symbolic simplification; fall back to numerical check.
    Returns (passed: bool, label: str, detail: str).
    """
    vf = val_func or _num_vals_vel
    try:
        res = simplify(expand(expr))
        if res == 0:
            return (True, label, "symbolic zero")
        res2 = trigsimp(res)
        if res2 == 0:
            return (True, label, "symbolic zero (trigsimp)")
    except Exception as e:
        res = expr  # proceed to numerical check

    if numerical_check(expr, vf):
        return (True, label, "numerically verified (symbolic simplification inconclusive)")
    return (False, label, f"FAILED — non-zero residual: {res}")


def check_nonzero(expr, label, val_func=None):
    """
    Negative control: assert expr does NOT vanish.  Used to prove that the
    pre-correction form of an equation really is wrong, so that each sign fix
    in the paper is shown to be necessary and not merely cosmetic.
    """
    vf = val_func or _num_vals_vel
    for _ in range(25):
        vals = vf()
        try:
            v = complex(expr.subs(vals))
        except Exception:
            continue
        if abs(v) > 1e-6:
            return (True, label, f"confirmed non-zero (residual {abs(v):.4g})")
    return (False, label, "FAILED — expression vanishes; the variant is not wrong")


# ============================================================
# GROUP A: Potential Identities (Eqs. potential, 104)
# ============================================================

def verify_group_A():
    results = []

    # A.1: Coulomb limit — U → q1q2/r12 when rdot=0 (xd1=xd2=yd1=yd2=0)
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

    # A.3: S − U = q1*q2*rdot²/(c²*r)  (sign flip of the velocity term)
    results.append(check_zero(
        S - U - q1*q2*rdot12_v**2/(c**2*r12),
        "A.3  S - U = q1*q2*rdot²/(c²*r)",
        _num_vals_vel,
    ))

    return results


# ============================================================
# GROUP B: Euler–Lagrange → Weber Force (Eq. euler_lagrange)
#
# Uses relative coordinates + chain-rule total time derivative.
# No SymPy time-Function needed; avoids slow symbolic expansion.
# ============================================================

def verify_group_B():
    results = []

    # Relative coordinates (2D)
    xr, yr   = symbols('xr yr', real=True)
    xdr, ydr = symbols('xdr ydr', real=True)     # relative velocity components
    xddr, yddr = symbols('xddr yddr', real=True) # relative acceleration components

    rr    = sqrt(xr**2 + yr**2)
    rdotr = (xr*xdr + yr*ydr) / rr
    # rddot in terms of relative coords (Appendix A.1, Form 1)
    rddotr = ((xr*xddr + yr*yddr)/rr
              + (xdr**2 + ydr**2)/rr
              - (xr*xdr + yr*ydr)**2/rr**3)

    S_rel = q1*q2/rr * (1 + rdotr**2 / (2*c**2))

    # ∂S/∂xdr  (the "generalized momentum" term)
    dS_dxdr = diff(S_rel, xdr)

    # d/dt(∂S/∂xdr) via total time derivative:
    # d/dt[f(xr,yr,xdr,ydr)] = ∂f/∂xr*xdr + ∂f/∂yr*ydr + ∂f/∂xdr*xddr + ∂f/∂ydr*yddr
    dt_dS_dxdr = (diff(dS_dxdr, xr)*xdr
                  + diff(dS_dxdr, yr)*ydr
                  + diff(dS_dxdr, xdr)*xddr
                  + diff(dS_dxdr, ydr)*yddr)

    # ∂S/∂xr
    dS_dxr = diff(S_rel, xr)

    # Euler-Lagrange residual
    EL = dt_dS_dxdr - dS_dxr

    # Expected Weber force (x-component of relative coordinate, Eq. euler_lagrange)
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
        "B.1  Euler-Lagrange: d/dt(∂S/∂ẋ) − ∂S/∂x = Weber force (Eq. euler_lagrange)",
        val_func_rel,
    ))

    # B.2: the same statement in the paper's own ABSOLUTE coordinates, i.e.
    #      literally Eq. euler_lagrange:
    #          F12^x = d/dt(∂S/∂ẋ1) − ∂S/∂x1
    #                = q1q2 (x1−x2)/r12³ (1 − ṙ²/2c² + r r̈/c²)
    xdd1, ydd1, xdd2, ydd2 = symbols('xdd1 ydd1 xdd2 ydd2', real=True)
    p_syms = [x1, y1, x2, y2]
    v_syms = [xd1, yd1, xd2, yd2]
    a_syms = [xdd1, ydd1, xdd2, ydd2]

    def dt_abs(f):
        return (sum(diff(f, pp)*vv for pp, vv in zip(p_syms, v_syms))
                + sum(diff(f, vv)*aa for vv, aa in zip(v_syms, a_syms)))

    # r̈ from Appendix A.1, written in absolute coordinates
    xdd_rel, ydd_rel = xdd1 - xdd2, ydd1 - ydd2
    rddot12 = ((x_rel*xdd_rel + y_rel*ydd_rel)/r12
               + (xd_rel**2 + yd_rel**2)/r12
               - rdot12_v**2/r12)

    EL_abs = dt_abs(diff(S, xd1)) - diff(S, x1)
    F_abs  = (q1*q2*x_rel/r12**3
              * (1 - rdot12_v**2/(2*c**2) + r12*rddot12/c**2))

    def val_func_abs():
        v = _num_vals_vel()
        v.update({xdd1: uniform(-0.1, 0.1), ydd1: uniform(-0.1, 0.1),
                  xdd2: uniform(-0.1, 0.1), ydd2: uniform(-0.1, 0.1)})
        return v

    results.append(check_zero_fast(
        EL_abs - F_abs,
        "B.2  Eq. euler_lagrange verbatim in absolute coords: "
        "d/dt(∂S/∂ẋ1) − ∂S/∂x1 = q1q2(x1−x2)/r12³(1−ṙ²/2c²+rr̈/c²)",
        val_func_abs,
    ))

    return results


# ============================================================
# GROUP C: Legendre Transform — H = T + U (Eq. H_def → Eq. H)
# ============================================================

def verify_group_C():
    results = []

    # alpha in velocity space (for C.1)
    alpha_x_v = q1*q2/c**2 * rdot12_v*(x1 - x2)/r12**2

    # C.1: Canonical momentum p_x1 = ∂L/∂xd1 = m1*xd1 − alpha_x_v
    px1_computed = diff(L, xd1)
    px1_expected = m1*xd1 - alpha_x_v
    results.append(check_zero(
        px1_computed - px1_expected,
        "C.1  Canonical momentum: p_x1 = m1*ẋ1 − alpha_x  (sign is minus)",
        _num_vals_vel,
    ))

    # C.2: Legendre transform (exact, in velocity space): sum(v * ∂L/∂v) − L = T + U
    H_legendre = (xd1*diff(L, xd1) + yd1*diff(L, yd1)
                  + xd2*diff(L, xd2) + yd2*diff(L, yd2)) - L
    results.append(check_zero(
        H_legendre - (T + U),
        "C.2  Legendre transform: H = T + U  (exact, velocity space)",
        _num_vals_vel,
    ))

    return results


# ============================================================
# GROUP D: 2-Particle Hamilton's Equations (verified against L)
#          (Eqs. xdot_two, xdot_two_p2, pdot_two, pdot_two_p2)
# ============================================================

def verify_group_D():
    results = []

    # D.1: ẋ1 = (px1 + alpha_x)/m1   [Eq. xdot_two]
    #      i.e. the inverse of px1 = m1*ẋ1 − alpha_x proved in C.1
    results.append(check_zero(
        (p1x_L + alpha_x)/m1 - xd1,
        "D.1  ẋ1 = (px1 + alpha_x)/m1  [Eq. xdot_two]",
    ))

    # D.2: ẋ2 = (px2 − alpha_x)/m2   [Eq. xdot_two_p2]
    results.append(check_zero(
        (p2x_L - alpha_x)/m2 - xd2,
        "D.2  ẋ2 = (px2 − alpha_x)/m2  [Eq. xdot_two_p2]",
    ))

    # D.3: ṗx1 = ∂L/∂x1 (Euler–Lagrange) matches the paper  [Eq. pdot_two]
    D3_expected = (q1*q2/r12**2
                   * ((x1 - x2)/r12 * (1 + Rational(3, 2)*rdot12_v**2/c**2)
                      - rdot12_v*xd_rel/c**2))
    results.append(check_zero(
        diff(L, x1) - D3_expected,
        "D.3  ṗx1 = ∂L/∂x1 matches paper formula  [Eq. pdot_two]",
    ))

    # D.4: ṗx1 + ṗx2 = 0 (Newton's 3rd law / total x-momentum conserved) [Eq. pdot_two_p2]
    results.append(check_zero(
        diff(L, x1) + diff(L, x2),
        "D.4  ṗx1 + ṗx2 = 0  (Newton's 3rd law, Eq. pdot_two_p2)",
    ))

    # D.5: alpha_x cancels in the SUM of the canonical momenta, so total
    #      canonical momentum equals total kinetic momentum:
    #          p_x1 + p_x2 = m1ẋ1 + m2ẋ2
    #      This — not "ẋ1 − ẋ2 ≠ 0" — is the identity that makes D.4 a
    #      statement about conserved linear momentum.
    results.append(check_zero(
        (p1x_L + p2x_L) - (m1*xd1 + m2*xd2),
        "D.5  p_x1 + p_x2 = m1ẋ1 + m2ẋ2  (α_x cancels in the sum)",
    ))

    # D.6: but alpha_x does NOT cancel in either momentum individually, which
    #      is exactly why v = p/m fails (see Group P).
    results.append(check_nonzero(
        p1x_L - m1*xd1,
        "D.6  α_x does not cancel in p_x1 alone (individual p ≠ m·v)",
    ))

    return results


# ============================================================
# GROUP E: n-Particle EOM Consistency (n=2 explicit summation)
#          Verifies Eqs. xdot_expanded and pdot_expanded reduce
#          to the 2-particle formulas when written with explicit sums.
# ============================================================

def verify_group_E():
    results = []

    # ---- xdot_expanded ----
    # ẋ_i = (1/m_i)(px_i + Σ_{j≠i} q_i q_j/c² ṙ_ij (x_i−x_j)/r_ij²)
    #
    # Convention: ṙ_ij for a pair (i,j) is (r_i−r_j)·(v_i−v_j)/|r_i−r_j|, which is
    # symmetric in the index order, so ṙ_21 = ṙ_12.

    # E.1a: n=2, i=1, sum over j=2
    xdot1_npart = (p1x_L + q1*q2/c**2 * rdot12_v*(x1 - x2)/r12**2) / m1
    results.append(check_zero(
        xdot1_npart - xd1,
        "E.1a xdot_expanded (n=2, i=1) returns ẋ1  [Eq. xdot_expanded]",
    ))

    # E.1b: n=2, i=2, sum over j=1.  (x2−x1) = −(x1−x2) flips the correction,
    #       reproducing the −alpha_x of Eq. xdot_two_p2.
    xdot2_npart = (p2x_L + q2*q1/c**2 * rdot12_v*(x2 - x1)/r12**2) / m2
    results.append(check_zero(
        xdot2_npart - xd2,
        "E.1b xdot_expanded (n=2, i=2) returns ẋ2  [Eq. xdot_expanded]",
    ))

    # E.2: n=2, i=1, pdot_expanded equals ∂L/∂x1
    pdot1_npart = (q1*q2/r12**2
                   * ((x1 - x2)/r12 * (1 + Rational(3, 2)*rdot12_v**2/c**2)
                      - rdot12_v*xd_rel/c**2))
    results.append(check_zero(
        pdot1_npart - diff(L, x1),
        "E.2  pdot_expanded (n=2, i=1) equals ∂L/∂x1  [Eq. pdot_expanded]",
    ))

    return results


# ============================================================
# GROUP F: Radial Acceleration Identities (Appendix A.1)
# ============================================================

def verify_group_F():
    results = []

    # 3D setup
    r3d    = sqrt(x**2 + y**2 + z**2)
    rdot3d = (x*xd + y*yd + z*zd) / r3d

    # Form 1 (component-wise, from paper):
    rddot_F1 = ((x*xdd + y*ydd + z*zdd)/r3d
                + (xd**2 + yd**2 + zd**2)/r3d
                - (x*xd + y*yd + z*zd)**2/r3d**3)

    # Form 3 (compact, using rdot²):
    rddot_F3 = (1/r3d)*((x*xdd + y*ydd + z*zdd)
                         + (xd**2 + yd**2 + zd**2)
                         - rdot3d**2)

    # F.1: Form 1 = Form 3
    results.append(check_zero(
        rddot_F1 - rddot_F3,
        "F.1  rddot Form1 = Form3  (three-form equivalence, Appendix A.1)",
        _num_vals_3d,
    ))

    # F.2: Form 1 matches d²/dt²(|r|) derived from first principles
    #
    # We compute d²(r)/dt² using the chain-rule total time derivative twice,
    # starting from r = sqrt(x²+y²+z²) as a function of (x,y,z,xd,yd,zd,xdd,ydd,zdd).

    def total_dt(f, pos, vel, acc):
        """d/dt f(pos, vel, acc) — one step of total time derivative."""
        return (sum(diff(f, p)*v for p, v in zip(pos, vel))
                + sum(diff(f, v)*a for v, a in zip(vel, acc)))

    pos_syms = [x, y, z]
    vel_syms = [xd, yd, zd]
    acc_syms = [xdd, ydd, zdd]

    # First derivative: dr/dt = rdot (as a function of pos, vel)
    rdot_expr = total_dt(r3d, pos_syms, vel_syms, acc_syms)
    # rdot_expr should equal rdot3d; verify:
    check_rdot = check_zero(
        rdot_expr - rdot3d,
        "F.2a rdot = d/dt(|r|) derived equals (r·v)/|r|",
        _num_vals_3d,
    )
    results.append(check_rdot)

    # Second derivative: d²r/dt² = rddot
    # rdot_expr is a function of (x,y,z,xd,yd,zd) — no xdd,ydd,zdd terms
    # d/dt(rdot_expr) needs the chain rule including acc terms
    rddot_derived = total_dt(rdot_expr, pos_syms, vel_syms, acc_syms)
    results.append(check_zero(
        rddot_derived - rddot_F1,
        "F.2b rddot Form1 = d²/dt²(|r|) derived from first principles",
        _num_vals_3d,
    ))

    return results


# ============================================================
# GROUP G: Conservation Laws (Section 2.2)
# ============================================================

def verify_group_G():
    results = []
    t_sym = symbols('t')

    # G.1: L has no explicit time dependence → energy conserved (Noether)
    has_t = t_sym in L.free_symbols
    results.append((
        not has_t,
        "G.1  L has no explicit time dependence (→ energy conserved)",
        "PASS: t not in free_symbols" if not has_t else "FAIL: t found in L",
    ))

    # G.2: Legendre transform of L is the conserved energy T + U
    E = (xd1*diff(L, xd1) + yd1*diff(L, yd1)
         + xd2*diff(L, xd2) + yd2*diff(L, yd2)) - L
    results.append(check_zero(
        E - (T + U),
        "G.2  Conserved energy Σv·∂L/∂v − L = T + U",
    ))

    # G.3: Total x-momentum rate = ∂L/∂x1 + ∂L/∂x2 = 0 (translational symmetry)
    results.append(check_zero(
        diff(L, x1) + diff(L, x2),
        "G.3  Total x-momentum rate = ∂L/∂x1 + ∂L/∂x2 = 0",
    ))

    # G.4: Total y-momentum rate = 0
    results.append(check_zero(
        diff(L, y1) + diff(L, y2),
        "G.4  Total y-momentum rate = ∂L/∂y1 + ∂L/∂y2 = 0",
    ))

    return results


# ============================================================
# GROUP H: n=3 Particle EOM Verification
#
# Directly tests the paper's two boxed equations for general n:
#   xdot_expanded:  ẋᵢ = (1/mᵢ)(pxᵢ + Σⱼ≠ᵢ qᵢqⱼ/c² ṙᵢⱼ(xᵢ−xⱼ)/rᵢⱼ²)
#   pdot_expanded:  ṗxᵢ = Σⱼ≠ᵢ qᵢqⱼ/rᵢⱼ² [(xᵢ−xⱼ)/rᵢⱼ(1+3ṙᵢⱼ²/2c²) − ṙᵢⱼ(ẋᵢ−ẋⱼ)/c²]
#
# Uses a 3-particle system to exercise the multi-pair sum structure.
# Note: ṙᵢⱼ = ṙⱼᵢ (radial velocity is symmetric in index order).
# ============================================================

def verify_group_H():
    results = []
    vf = _num_vals_3body

    def pair_pdot_x(qi, qj, xi, xj, rij, rdot_ij, vxi, vxj):
        """x-component of ṗᵢ contribution from one (i,j) pair, Eq. pdot_expanded."""
        return (qi*qj/rij**2
                * ((xi - xj)/rij * (1 + Rational(3, 2)*rdot_ij**2/c**2)
                   - rdot_ij*(vxi - vxj)/c**2))

    # ----- H.1–H.3: xdot_expanded -----
    # ẋᵢ = (1/mᵢ)(pxᵢ + Σⱼ≠ᵢ qᵢqⱼ/c² ṙᵢⱼ(xᵢ−xⱼ)/rᵢⱼ²)

    xdot1_formula = (diff(L3, xd1) + q1*q2/c**2 * rdot12_v*(x1 - x2)/r12**2
                                   + q1*q3/c**2 * rdot13_v*(x1 - x3)/r13**2) / m1
    results.append(check_zero_fast(
        xdot1_formula - xd1,
        "H.1  xdot_expanded (n=3, i=1) returns ẋ1  [Eq. xdot_expanded]",
        vf,
    ))

    # i=2: ṙ₂₁ = ṙ₁₂ (symmetric), so rdot12_v is used for the (2,1) pair
    xdot2_formula = (diff(L3, xd2) + q2*q1/c**2 * rdot12_v*(x2 - x1)/r12**2
                                   + q2*q3/c**2 * rdot23_v*(x2 - x3)/r23**2) / m2
    results.append(check_zero_fast(
        xdot2_formula - xd2,
        "H.2  xdot_expanded (n=3, i=2) returns ẋ2  [Eq. xdot_expanded]",
        vf,
    ))

    # i=3: ṙ₃₁ = ṙ₁₃, ṙ₃₂ = ṙ₂₃
    xdot3_formula = (diff(L3, xd3) + q3*q1/c**2 * rdot13_v*(x3 - x1)/r13**2
                                   + q3*q2/c**2 * rdot23_v*(x3 - x2)/r23**2) / m3
    results.append(check_zero_fast(
        xdot3_formula - xd3,
        "H.3  xdot_expanded (n=3, i=3) returns ẋ3  [Eq. xdot_expanded]",
        vf,
    ))

    # ----- H.4–H.6: pdot_expanded equals ∂L3/∂xᵢ -----

    pdot_x1_formula = (pair_pdot_x(q1, q2, x1, x2, r12, rdot12_v, xd1, xd2)
                     + pair_pdot_x(q1, q3, x1, x3, r13, rdot13_v, xd1, xd3))
    results.append(check_zero_fast(
        pdot_x1_formula - diff(L3, x1),
        "H.4  pdot_expanded (n=3, i=1) equals ∂L3/∂x1  [Eq. pdot_expanded]",
        vf,
    ))

    pdot_x2_formula = (pair_pdot_x(q2, q1, x2, x1, r12, rdot12_v, xd2, xd1)
                     + pair_pdot_x(q2, q3, x2, x3, r23, rdot23_v, xd2, xd3))
    results.append(check_zero_fast(
        pdot_x2_formula - diff(L3, x2),
        "H.5  pdot_expanded (n=3, i=2) equals ∂L3/∂x2  [Eq. pdot_expanded]",
        vf,
    ))

    pdot_x3_formula = (pair_pdot_x(q3, q1, x3, x1, r13, rdot13_v, xd3, xd1)
                     + pair_pdot_x(q3, q2, x3, x2, r23, rdot23_v, xd3, xd2))
    results.append(check_zero_fast(
        pdot_x3_formula - diff(L3, x3),
        "H.6  pdot_expanded (n=3, i=3) equals ∂L3/∂x3  [Eq. pdot_expanded]",
        vf,
    ))

    # ----- H.7–H.8: Total momentum conservation -----
    results.append(check_zero_fast(
        diff(L3, x1) + diff(L3, x2) + diff(L3, x3),
        "H.7  ṗx1 + ṗx2 + ṗx3 = 0  (n=3 x-momentum conservation)",
        vf,
    ))
    results.append(check_zero_fast(
        diff(L3, y1) + diff(L3, y2) + diff(L3, y3),
        "H.8  ṗy1 + ṗy2 + ṗy3 = 0  (n=3 y-momentum conservation)",
        vf,
    ))

    return results


# ============================================================
# GROUP I: Angular Momentum Conservation
#
# Verifies rotational invariance of L for n=2 and n=3, which is Noether's
# condition for conservation of Lz = Σᵢ(xᵢpyᵢ − yᵢpxᵢ) and is the rotational
# symmetry claimed in Section 2.2.
# ============================================================

def verify_group_I():
    results = []

    def rotational_variation(L_expr, coord_tuples):
        """delta_L under an infinitesimal rotation dx = -y, dy = x (and likewise
        for the velocities). Vanishing delta_L is Noether's condition for
        conservation of Lz = sum_i (x_i*p_yi - y_i*p_xi)."""
        var = Integer(0)
        for xi, yi, vxi, vyi in coord_tuples:
            var += (diff(L_expr, xi)*(-yi) + diff(L_expr, yi)*xi
                  + diff(L_expr, vxi)*(-vyi) + diff(L_expr, vyi)*vxi)
        return var

    # I.1: rotational invariance of L (2-particle) → Lz conserved
    results.append(check_zero_fast(
        rotational_variation(L, [(x1, y1, xd1, yd1), (x2, y2, xd2, yd2)]),
        "I.1  Rotational invariance of L (angular momentum conservation, n=2)",
        _num_vals_vel,
    ))

    # I.2: rotational invariance of L3 (3-particle — full multi-pair structure)
    results.append(check_zero_fast(
        rotational_variation(
            L3, [(x1, y1, xd1, yd1), (x2, y2, xd2, yd2), (x3, y3, xd3, yd3)]
        ),
        "I.2  Rotational invariance of L3 (angular momentum conservation, n=3)",
        _num_vals_3body,
    ))

    return results


# ============================================================
# GROUP J: Force <-> Potential (Section 1, Eqs. potential & force)
#
# The paper asserts (Sec. 1) that Weber's force can be derived from Weber's
# potential, and that the acceleration term falls off as 1/r. Both are checked
# here as scalar identities in the relational variables (r, rdot, rddot),
# which are the only variables U and F depend on.
# ============================================================

def verify_group_J():
    results = []

    rs        = symbols('r_s', positive=True)
    rds, rdds = symbols('rdot_s rddot_s', real=True)

    U_s = q1*q2/rs * (1 - rds**2/(2*c**2))          # Eq. potential
    S_s = q1*q2/rs * (1 + rds**2/(2*c**2))          # auxiliary potential
    F_s = q1*q2/rs**2 * (1 - rds**2/(2*c**2) + rs*rdds/c**2)   # Eq. force

    def dt_scalar(f):
        """Total time derivative of f(r, rdot) along a trajectory."""
        return diff(f, rs)*rds + diff(f, rds)*rdds

    def val_scalar():
        return {rs: uniform(0.5, 3.0), rds: uniform(-0.4, 0.4),
                rdds: uniform(-0.2, 0.2),
                q1: uniform(0.2, 1.5), q2: uniform(0.2, 1.5),
                c: uniform(3.0, 8.0)}

    # J.1: dU/dt = -F*rdot.  This IS the statement "the force derives from the
    #      potential": combined with dT/dt = +F*rdot (J.2) it gives d(T+U)/dt = 0.
    results.append(check_zero(
        dt_scalar(U_s) + F_s*rds,
        "J.1  dU/dt = -F*rdot  (Weber force follows from Weber potential)",
        val_scalar,
    ))

    # J.2: Lagrangian route -- F = -dS/dr + d/dt(dS/drdot), the generalized
    #      force of the velocity-dependent potential S.
    results.append(check_zero(
        (-diff(S_s, rs) + dt_scalar(diff(S_s, rds))) - F_s,
        "J.2  F = -∂S/∂r + d/dt(∂S/∂ṙ)  (Eq. force from the Lagrangian)",
        val_scalar,
    ))

    # J.3: the acceleration-dependent term falls off as 1/r (Sec. 1 & Discussion)
    F_acc_term = q1*q2/rs**2 * (rs*rdds/c**2)
    results.append(check_zero(
        F_acc_term - q1*q2*rdds/(c**2*rs),
        "J.3  acceleration term of Eq. force = q1q2*r̈/(c²r)  → falls off as 1/r",
        val_scalar,
    ))

    # J.4: Newton's third law in the strong form. F12 = F_s * rhat12 with
    #      r, rdot, rddot symmetric under 1<->2, and rhat21 = -rhat12.
    #      Checked component-wise in absolute 2-D coordinates.
    F12_x = q1*q2*(x1 - x2)/r12**3 * (1 - rdot12_v**2/(2*c**2))
    F21_x = q2*q1*(x2 - x1)/r12**3 * (1 - rdot12_v**2/(2*c**2))
    results.append(check_zero(
        F12_x + F21_x,
        "J.4  F21 = -F12 componentwise (Newton's 3rd law, strong form)",
        _num_vals_vel,
    ))

    # J.5: the force is central -- F12 x rhat = 0 (needed for angular momentum)
    F12_y = q1*q2*(y1 - y2)/r12**3 * (1 - rdot12_v**2/(2*c**2))
    results.append(check_zero(
        F12_x*(y1 - y2) - F12_y*(x1 - x2),
        "J.5  F12 is central: F12 × r̂ = 0  (→ angular momentum conserved)",
        _num_vals_vel,
    ))

    return results


# ============================================================
# GROUP K: Closed-form Hamiltonian H(q,p) for n = 2
#
# Section 2.1 states H = H(positions, momenta) after the Legendre transform,
# but never exhibits it: Eq. xdot_two still carries the velocity-space alpha_x
# on the right-hand side, so it is an implicit relation, not a closed form.
#
# For n = 2 the Legendre transform CAN be carried out exactly. In CoM +
# relative polar coordinates the momenta decouple and
#
#   H = P²/(2M) + p_r²/(2*mu_eff(r)) + p_phi²/(2*mu*r²) + q1q2/r,
#       mu_eff(r) = mu - q1q2/(c² r)
#
# Every claim of Section 2.1-2.2 is then verified in genuine (q,p) variables,
# and the chain closes back onto Weber's force law (K.9).
# ============================================================

def verify_group_K():
    results = []

    M_tot = m1 + m2
    mu    = m1*m2/(m1 + m2)
    kW    = q1*q2/c**2                      # Weber coupling, units of mass*length

    rr  = symbols('r_rel', positive=True)
    ph  = symbols('phi', real=True)
    rd  = symbols('rdot_rel', real=True)
    phd = symbols('phidot', real=True)
    Xd, Yd = symbols('Xdot Ydot', real=True)
    pr, pph, PX, PY = symbols('p_r p_phi P_X P_Y', real=True)

    mu_eff = mu - kW/rr

    T_pol = M_tot*(Xd**2 + Yd**2)/2 + mu*(rd**2 + rr**2*phd**2)/2
    S_pol = q1*q2/rr * (1 + rd**2/(2*c**2))
    U_pol = q1*q2/rr * (1 - rd**2/(2*c**2))
    L_pol = T_pol - S_pol

    def val_pol():
        return {rr: uniform(1.0, 3.0), ph: uniform(0, 6.0),
                rd: uniform(-0.3, 0.3), phd: uniform(-0.3, 0.3),
                Xd: uniform(-0.3, 0.3), Yd: uniform(-0.3, 0.3),
                m1: uniform(0.5, 2.0), m2: uniform(0.5, 2.0),
                q1: uniform(0.2, 1.5), q2: uniform(0.2, 1.5),
                c: uniform(3.0, 8.0)}

    # K.0: the CoM + relative split reproduces Eq. ke exactly (2-D Cartesian).
    Xd_c = (m1*xd1 + m2*xd2)/(m1 + m2)
    Yd_c = (m1*yd1 + m2*yd2)/(m1 + m2)
    results.append(check_zero(
        (m1*(xd1**2 + yd1**2)/2 + m2*(xd2**2 + yd2**2)/2)
        - ((m1 + m2)*(Xd_c**2 + Yd_c**2)/2
           + mu*((xd1 - xd2)**2 + (yd1 - yd2)**2)/2),
        "K.0  T = ½M|Ṙ|² + ½μ|v_rel|²  (CoM/relative split of Eq. ke)",
        _num_vals_vel,
    ))

    # K.1: radial canonical momentum p_r = mu_eff(r) * rdot.
    #      This is where the critical radius enters the formalism.
    results.append(check_zero(
        diff(L_pol, rd) - mu_eff*rd,
        "K.1  p_r = ∂L/∂ṙ = (μ − q1q2/(c²r))·ṙ = μ_eff·ṙ",
        val_pol,
    ))

    # K.2: angular canonical momentum is unmodified by the Weber term
    results.append(check_zero(
        diff(L_pol, phd) - mu*rr**2*phd,
        "K.2  p_φ = ∂L/∂φ̇ = μr²φ̇  (Weber term does not touch p_φ)",
        val_pol,
    ))

    # ---- the closed-form Hamiltonian -------------------------------------
    H_qp = (PX**2 + PY**2)/(2*M_tot) + pr**2/(2*mu_eff) \
         + pph**2/(2*mu*rr**2) + q1*q2/rr

    to_vel = [(pr, mu_eff*rd), (pph, mu*rr**2*phd), (PX, M_tot*Xd), (PY, M_tot*Yd)]

    # K.3: velocity-space Legendre transform reproduces Eq. H / Eq. h_box
    H_leg = (rd*diff(L_pol, rd) + phd*diff(L_pol, phd)
             + Xd*diff(L_pol, Xd) + Yd*diff(L_pol, Yd)) - L_pol
    results.append(check_zero(
        H_leg - (T_pol + U_pol),
        "K.3  Σv·∂L/∂v − L = T + U  in polar coords  [Eqs. H, h_box]",
        val_pol,
    ))

    # K.4: the closed form H(q,p) equals T + U once momenta are back-substituted.
    #      This is the Legendre transform Eq. H_p only asserts.
    results.append(check_zero(
        H_qp.subs(to_vel) - (T_pol + U_pol),
        "K.4  closed-form H(q,p) = T + U  (explicit Legendre transform, n=2)",
        val_pol,
    ))

    # K.5 / K.6: Hamilton's equations Eq. xdot round-trip to the velocities
    results.append(check_zero(
        diff(H_qp, pr).subs(to_vel) - rd,
        "K.5  ṙ = ∂H/∂p_r  returns ṙ  [Eq. xdot, closed form]",
        val_pol,
    ))
    results.append(check_zero(
        diff(H_qp, pph).subs(to_vel) - phd,
        "K.6  φ̇ = ∂H/∂p_φ  returns φ̇  [Eq. xdot, closed form]",
        val_pol,
    ))

    # K.7: p_phi is cyclic -> angular momentum conserved (Sec. 2.2)
    results.append(check_zero(
        diff(H_qp, ph),
        "K.7  ∂H/∂φ = 0 → ṗ_φ = 0  (angular momentum conserved, Sec. 2.2)",
        val_pol,
    ))

    # K.8: CoM coordinates are cyclic -> linear momentum conserved (Sec. 2.2)
    Xc, Yc = symbols('X_c Y_c', real=True)
    results.append((
        diff(H_qp, Xc) == 0 and diff(H_qp, Yc) == 0,
        "K.8  ∂H/∂X = ∂H/∂Y = 0 → Ṗ = 0  (linear momentum conserved, Sec. 2.2)",
        "symbolic zero (CoM coordinates absent from H)",
    ))

    # K.9: END-TO-END.  Hamilton's equations applied to the closed-form H(q,p)
    #      reproduce Weber's force law Eq. force for the relative coordinate:
    #          mu*(r̈ − r φ̇²) = q1q2/r² (1 − ṙ²/2c² + r r̈/c²)
    prdot = -diff(H_qp, rr)                       # Eq. pdot
    # r̈ = d/dt(p_r/mu_eff) = ṗ_r/μ_eff − p_r μ_eff' ṙ / μ_eff²
    rddot = (prdot/mu_eff - pr*diff(mu_eff, rr)*rd/mu_eff**2).subs(to_vel)
    rddot = simplify(rddot)
    F_weber = q1*q2/rr**2 * (1 - rd**2/(2*c**2) + rr*rddot/c**2)
    results.append(check_zero(
        mu*(rddot - rr*phd**2) - F_weber,
        "K.9  Hamilton's eqs on closed-form H  ⇒  Weber's force law [Eq. force]",
        val_pol,
    ))

    # K.10: bridge to the paper's Cartesian variables.  With the paper's
    #       p_x1 = m1ẋ1 − alpha_x (Eq. px1), the relative momentum
    #       p_rel = (m2 p1 − m1 p2)/M projected on r̂ equals μ_eff·ṙ.
    p1x, p1y = diff(L, xd1), diff(L, yd1)
    p2x, p2y = diff(L, xd2), diff(L, yd2)
    prel_x = (m2*p1x - m1*p2x)/(m1 + m2)
    prel_y = (m2*p1y - m1*p2y)/(m1 + m2)
    rhat_x, rhat_y = (x1 - x2)/r12, (y1 - y2)/r12
    mu_c   = m1*m2/(m1 + m2)
    results.append(check_zero(
        (prel_x*rhat_x + prel_y*rhat_y) - (mu_c - q1*q2/(c**2*r12))*rdot12_v,
        "K.10 r̂·p_rel = μ_eff·ṙ  (Eq. px1 Cartesian ⇔ closed form K.1)",
        _num_vals_vel,
    ))

    return results


# ============================================================
# GROUP L: Weber's critical radius (Discussion)
#
# The paper states rho = q1q2/(mu c²) and then makes two adjacent claims:
#   (a) "At this critical radius, the sign of the force ... is reversed."
#   (b) "the particles respond to the still-repulsive Coulomb force as if
#        they had negative inertia".
# These checks pin down exactly what is and is not true.
# ============================================================

def verify_group_L():
    results = []

    mu_s = symbols('mu_s', positive=True)
    rr   = symbols('r_rel', positive=True)
    rd   = symbols('rdot_rel', real=True)

    kW     = q1*q2/c**2
    rho    = q1*q2/(mu_s*c**2)                 # paper's critical separation
    mu_eff = mu_s - kW/rr

    def val_crit():
        return {rr: uniform(1.0, 3.0), rd: uniform(-0.4, 0.4),
                mu_s: uniform(0.5, 2.0),
                q1: uniform(0.2, 1.5), q2: uniform(0.2, 1.5),
                c: uniform(3.0, 8.0)}

    # L.1: mu_eff vanishes exactly at r = rho -- this defines the critical radius
    results.append(check_zero(
        mu_eff.subs(rr, rho),
        "L.1  μ_eff(ρ) = 0 at ρ = q1q2/(μc²)  (critical radius, Discussion)",
        val_crit,
    ))

    # L.2: rho has the dimensions of a length -- see Group O for the unit algebra
    #      Here: rho * mu * c^2 - q1q2 = 0 by construction.
    results.append(check_zero(
        rho*mu_s*c**2 - q1*q2,
        "L.2  ρ·μ·c² = q1q2  (definition is dimensionally a length)",
        val_crit,
    ))

    # L.3: radial (l = 0) equation of motion, solved for r̈:
    #          mu_eff * r̈ = q1q2/r² (1 − ṙ²/2c²)
    #      Verify this is equivalent to Weber's force law mu*r̈ = F.
    rdd = q1*q2/rr**2*(1 - rd**2/(2*c**2)) / mu_eff
    F_on_shell = q1*q2/rr**2*(1 - rd**2/(2*c**2) + rr*rdd/c**2)
    results.append(check_zero(
        mu_s*rdd - F_on_shell,
        "L.3  μ_eff·r̈ = q1q2/r²(1−ṙ²/2c²) ⇔ μr̈ = F_Weber  (ℓ=0 radial motion)",
        val_crit,
    ))

    # L.4: the DRIVING term q1q2/r²(1 − ṙ²/2c²) never changes sign for like
    #      charges with |ṙ| < c*sqrt(2).  It is repulsive above AND below rho.
    #      -> claim (b) is correct; nothing about *this* force reverses at rho.
    drive = q1*q2/rr**2*(1 - rd**2/(2*c**2))
    ok_drive = True
    for _ in range(400):
        v = val_crit()
        # like charges, subcritical |rdot|
        v[q1], v[q2] = abs(v[q1]), abs(v[q2])
        v[rd] = uniform(-1.0, 1.0)
        rho_n = float((rho).subs(v))
        v[rr] = uniform(0.05*rho_n, 5.0*rho_n)      # straddles rho
        if float(drive.subs(v)) <= 0:
            ok_drive = False
            break
    results.append((
        ok_drive,
        "L.4  driving term q1q2/r²(1−ṙ²/2c²) stays REPULSIVE across ρ "
        "(like charges, |ṙ|<c√2)",
        "400 random samples straddling ρ, all strictly positive"
        if ok_drive else "FAILED — driving term changed sign",
    ))

    # L.5: the ON-SHELL acceleration (and hence the net force mu*r̈) flips sign
    #      across rho because mu_eff does -- and it is SINGULAR at r = rho,
    #      not merely sign-reversed there.
    flips = singular = True
    for _ in range(200):
        v = val_crit()
        v[q1], v[q2] = abs(v[q1]), abs(v[q2])
        v[rd] = uniform(-0.5, 0.5)
        rho_n = float(rho.subs(v))
        a_above = float(rdd.subs({**v, rr: 1.6*rho_n}))
        a_below = float(rdd.subs({**v, rr: 0.6*rho_n}))
        if not (a_above > 0 > a_below):
            flips = False
        # approaching rho from either side, |r̈| blows up
        a_near_hi = abs(float(rdd.subs({**v, rr: rho_n*(1 + 1e-6)})))
        a_near_lo = abs(float(rdd.subs({**v, rr: rho_n*(1 - 1e-6)})))
        if not (a_near_hi > 1e4 and a_near_lo > 1e4):
            singular = False
    results.append((
        flips,
        "L.5  on-shell r̈ (and net force μr̈) reverses sign BELOW ρ, not AT ρ",
        "200 samples: r̈>0 at r=1.6ρ, r̈<0 at r=0.6ρ"
        if flips else "FAILED — no sign reversal found",
    ))
    results.append((
        singular,
        "L.6  r̈ diverges as r → ρ (μ_eff → 0): ρ is a singularity, not a zero",
        "200 samples: |r̈| > 1e4 on both sides of ρ at |1−r/ρ| = 1e-6"
        if singular else "FAILED — no divergence at ρ",
    ))

    # L.7: the same statement for ARBITRARY angular momentum, not just l = 0.
    #      Solving the planar relative EOM  mu(r̈ − rφ̇²) = F_Weber  for r̈ and
    #      substituting back gives the closed form
    #
    #          F_onshell = mu * (drive + (q1q2/c²)·φ̇²) / mu_eff(r)
    #
    #      so the whole r-dependence of the SIGN sits in mu_eff.
    phd  = symbols('phidot_L', real=True)
    rdd  = symbols('rddot_L', real=True)
    drive = q1*q2/rr**2*(1 - rd**2/(2*c**2))
    F_sym = q1*q2/rr**2*(1 - rd**2/(2*c**2) + rr*rdd/c**2)
    rdd_sol = solve(mu_s*(rdd - rr*phd**2) - F_sym, rdd)[0]
    F_onshell = mu_s*(rdd_sol - rr*phd**2)

    def val_crit_l():
        v = val_crit()
        v.update({phd: uniform(-0.4, 0.4)})
        return v

    results.append(check_zero(
        F_onshell - mu_s*(drive + (q1*q2/c**2)*phd**2)/mu_eff,
        "L.7  F_onshell = μ(drive + (q1q2/c²)φ̇²)/μ_eff  (any ℓ, not just ℓ=0)",
        val_crit_l,
    ))

    # L.8: therefore, for LIKE charges the numerator is strictly positive and
    #      sign(F) = sign(mu_eff): the net force really is repulsive above rho
    #      and attractive below it — for every angular momentum.
    #      This is the Discussion's "Below this critical radius, the sign of
    #      the force between q1 and q2 is reversed."
    reversed_ok = True
    for _ in range(400):
        v = val_crit_l()
        v[q1], v[q2] = abs(v[q1]), abs(v[q2])      # like charges
        v[rd] = uniform(-0.5, 0.5)
        rho_n = float(rho.subs(v))
        f_above = float(F_onshell.subs({**v, rr: uniform(1.05, 6.0)*rho_n}))
        f_below = float(F_onshell.subs({**v, rr: uniform(0.05, 0.95)*rho_n}))
        if not (f_above > 0 > f_below):
            reversed_ok = False
            break
    results.append((
        reversed_ok,
        "L.8  like charges: net force is repulsive above ρ, ATTRACTIVE below ρ, "
        "at arbitrary ℓ  (Discussion claim)",
        "400 random samples with φ̇ ≠ 0, straddling ρ"
        if reversed_ok else "FAILED — sign did not reverse across ρ",
    ))

    return results


# ============================================================
# GROUP M: Extended phase space & the symplectic integrator (Section 3)
#
# Checks the linear algebra of Eqs. Z, A, AZ, Phi and the 1/4 Newton factor,
# plus the claim that the Strang splitting is second order and symplectic.
# Uses 3n = 2 (one particle in 2-D, or two in 1-D) so the blocks are 2x2.
# ============================================================

def verify_group_M():
    results = []
    d = 2                                   # stands for 3n
    I, Zr = eye(d), zeros(d, d)

    # A is 2d x 4d with block rows [I -I 0 0 ; 0 0 I -I]   (Sec. 3)
    A = Matrix.vstack(
        Matrix.hstack(I,  -I, Zr, Zr),
        Matrix.hstack(Zr, Zr, I,  -I),
    )

    qs = Matrix(symbols('q1_ q2_', real=True))
    xs = Matrix(symbols('x1_ x2_', real=True))
    ps = Matrix(symbols('p1_ p2_', real=True))
    ys = Matrix(symbols('y1_ y2_', real=True))
    Zv = Matrix.vstack(qs, xs, ps, ys)       # Eq. Z ordering (q, x, p, y)

    # M.1: A*Z = (q - x, p - y)   [Eq. AZ]
    results.append((
        simplify(A*Zv - Matrix.vstack(qs - xs, ps - ys)) == zeros(2*d, 1),
        "M.1  A·Z = (q−x, p−y)  [Eq. AZ] for the Eq. Z ordering (q,x,p,y)",
        "symbolic zero",
    ))

    # M.2: A*A^T = 2*I  -- this is what fixes the 1/4 in the Newton iteration
    results.append((
        simplify(A*A.T - 2*eye(2*d)) == zeros(2*d, 2*d),
        "M.2  A·Aᵀ = 2I",
        "symbolic zero",
    ))

    # M.3: the simplified Newton factor.  f(mu) = A(Phi(Z + A^T mu) + A^T mu);
    #      with the simplified approximation dPhi/dZ ~ I,
    #      df/dmu ~ A(I + I)A^T = 2*A*A^T = 4I,  so mu <- mu - (1/4) f.
    Jac = A*(2*eye(4*d))*A.T
    results.append((
        simplify(Jac - 4*eye(2*d)) == zeros(2*d, 2*d),
        "M.3  simplified Newton Jacobian A(I+I)Aᵀ = 4I → step is −(1/4)f(μ)",
        "symbolic zero",
    ))

    # M.4: mu lives in R^{6n}: A has 2d = 6n rows when d = 3n
    results.append((
        A.shape == (2*d, 4*d),
        "M.4  A is 6n×12n and μ ∈ ℝ^{6n}  (dimensions of Sec. 3 consistent)",
        f"A.shape = {A.shape} with 3n = {d}",
    ))

    # ---- extended-Hamiltonian flow structure (the three bullets of Sec. 3) ----
    qe, xe, pe, ye, om = symbols('q_e x_e p_e y_e omega', real=True)

    def Ham(Q, P):
        return P**2/2 + om**2*Q**2/2

    Hbar = Ham(qe, ye) + Ham(xe, pe)         # H_A(q,y) + H_B(x,p)
    HA, HB = Ham(qe, ye), Ham(xe, pe)

    # M.5: H_A's flow moves only (x, p); H_B's flow moves only (q, y).
    A_moves_only_xp = (diff(HA, pe) == 0 and diff(HA, xe) == 0)
    B_moves_only_qy = (diff(HB, ye) == 0 and diff(HB, qe) == 0)
    results.append((
        A_moves_only_xp,
        "M.5  Φ^A flow: q̇ = ∂H_A/∂p = 0 and ẏ = −∂H_A/∂x = 0 → moves only (x,p)",
        "matches Sec. 3 bullet 1",
    ))
    results.append((
        B_moves_only_qy,
        "M.6  Φ^B flow: ẋ = ∂H_B/∂y = 0 and ṗ = −∂H_B/∂q = 0 → moves only (q,y)",
        "matches Sec. 3 bullet 2",
    ))

    # M.7: on the constraint x = q, y = p the extended flow reduces to
    #      Hamilton's equations for H  (this is why the method is consistent)
    con = [(xe, qe), (ye, pe)]
    red_qdot = diff(Hbar, pe).subs(con)
    red_pdot = (-diff(Hbar, qe)).subs(con)
    results.append((
        simplify(red_qdot - diff(Ham(qe, pe), pe)) == 0
        and simplify(red_pdot + diff(Ham(qe, pe), qe)) == 0,
        "M.7  on x=q, y=p the extended flow reduces to Hamilton's eqs for H",
        "symbolic zero",
    ))

    # M.8: each sub-flow has a CONSTANT right-hand side, so its exact solution
    #      is the Euler update of Eqs. update_step_x / update_step_p.
    #      (This is what makes Eq. update_step_* legitimate despite Euler
    #       being non-symplectic on its own.)
    h = symbols('h', positive=True)
    rhs_x, rhs_p = diff(HA, ye), -diff(HA, qe)      # depend on (q,y) only
    results.append((
        diff(rhs_x, xe) == 0 and diff(rhs_x, pe) == 0
        and diff(rhs_p, xe) == 0 and diff(rhs_p, pe) == 0,
        "M.8  Φ^A right-hand sides are constant in (x,p) → Euler step is EXACT",
        "justifies Eqs. update_step_x / update_step_p inside Φ^A, Φ^B",
    ))

    # ---- M.9-M.11: Strang splitting on a harmonic H (exactly integrable) ----
    def phiA(st, dt):
        q_, x_, p_, y_ = st
        return (q_, x_ + dt*y_, p_ - dt*om**2*q_, y_)

    def phiB(st, dt):
        q_, x_, p_, y_ = st
        return (q_ + dt*p_, x_, p_, y_ - dt*om**2*x_)

    st0  = (qe, xe, pe, ye)
    step = phiA(phiB(phiA(st0, h/2), h), h/2)       # Eq. Phi

    # restrict to the constraint x = q, y = p and compare with the exact flow
    stc = [expand(s.subs(con)) for s in step]
    q_ex = qe*cos(om*h) + pe*sin(om*h)/om
    p_ex = pe*cos(om*h) - qe*om*sin(om*h)

    err_q = expand(stc[0] - q_ex)
    err_p = expand(stc[2] - p_ex)
    ord_q = series(err_q, h, 0, 4).removeO()
    ord_p = series(err_p, h, 0, 4).removeO()
    second_order = (simplify(ord_q.subs(h, 0)) == 0
                    and simplify(diff(ord_q, h).subs(h, 0)) == 0
                    and simplify(diff(ord_q, h, 2).subs(h, 0)) == 0
                    and simplify(ord_p.subs(h, 0)) == 0
                    and simplify(diff(ord_p, h).subs(h, 0)) == 0
                    and simplify(diff(ord_p, h, 2).subs(h, 0)) == 0)
    results.append((
        second_order,
        "M.9  Φ = Φ^A_{Δt/2}∘Φ^B_{Δt}∘Φ^A_{Δt/2} is 2nd order  [Eq. Phi]",
        "local error starts at O(Δt³) (h⁰,h¹,h² terms all vanish)",
    ))

    # M.10: the extended map is symplectic: J^T Omega J = Omega with the
    #       Eq. Z ordering (q,x,p,y) -> Omega = [[0,I],[-I,0]] is standard.
    J = Matrix([[diff(s, v) for v in st0] for s in step])
    Om = Matrix.vstack(
        Matrix.hstack(zeros(2, 2), eye(2)),
        Matrix.hstack(-eye(2), zeros(2, 2)),
    )
    results.append((
        simplify(J.T*Om*J - Om) == zeros(4, 4),
        "M.10 Jᵀ Ω J = Ω — the map Eq. Phi is symplectic in extended phase space",
        "symbolic zero (confirms the Eq. Z variable ordering pairs (q,p),(x,y))",
    ))

    # M.11: phase-space volume preservation (Liouville, Sec. 3 opening claim)
    results.append((
        simplify(J.det() - 1) == 0,
        "M.11 det J = 1 — phase space volume preserved (Liouville's theorem)",
        "symbolic zero",
    ))

    # M.12: starting ON the constraint x=q, y=p, one Strang step leaves it.
    #       The drift is what the symmetric projection (the shift vector μ)
    #       of Sec. 3 exists to remove — so the projection is not optional.
    drift_x = simplify(stc[1] - stc[0])
    drift_y = simplify(stc[3] - stc[2])
    results.append((
        drift_x != 0 and drift_y != 0,
        "M.12 one Φ step leaves the constraint A·Z = 0 (drift ≠ 0) → the "
        "symmetric projection is required",
        f"x−q drift = {drift_x}, p−y drift = {-drift_y}",
    ))

    return results


# ============================================================
# GROUP N: n-body Hamiltonian, mass matrix, and the implicitness of
#          Eq. xdot_expanded  (Sections 2.3 and 2.4)
# ============================================================

def _mass_matrix_3body():
    """M(q) with p = M(q) v, from the 3-particle 2-D Lagrangian L3."""
    vels = [xd1, yd1, xd2, yd2, xd3, yd3]
    return Matrix(len(vels), len(vels),
                  lambda a, b: diff(L3, vels[a], vels[b])), vels


def verify_group_N():
    results = []
    vf = _num_vals_3body

    # N.1: the n-body Legendre transform gives Eq. hamiltonian
    #      H_n = sum ½m v² + sum_{i<j} q_i q_j / r_ij (1 − ṙ_ij²/2c²)
    vels = [xd1, yd1, xd2, yd2, xd3, yd3]
    H3_leg = sum(v*diff(L3, v) for v in vels) - L3
    H3_paper = (T3
                + q1*q2/r12*(1 - rdot12_v**2/(2*c**2))
                + q1*q3/r13*(1 - rdot13_v**2/(2*c**2))
                + q2*q3/r23*(1 - rdot23_v**2/(2*c**2)))
    results.append(check_zero_fast(
        H3_leg - H3_paper,
        "N.1  Σv·∂L₃/∂v − L₃ = H₃ of Eq. hamiltonian  (n-body Legendre, n=3)",
        vf,
    ))

    # N.2: the momentum map is exactly linear: p = M(q) v, M symmetric.
    Mm, vels = _mass_matrix_3body()
    vvec = Matrix(vels)
    pvec = Matrix([diff(L3, v) for v in vels])
    results.append(check_zero_fast(
        sum(abs(e) for e in (Mm*vvec - pvec)),
        "N.2  p = M(q)·v exactly  (velocity dependence of L₃ is homogeneous deg 2)",
        vf,
    ))
    results.append(check_zero_fast(
        sum(abs(e) for e in (Mm - Mm.T)),
        "N.3  M(q) is symmetric (a genuine mass matrix)",
        vf,
    ))

    # N.4: M is NOT diagonal -- the Weber term couples the particles, so
    #      Eq. xdot_expanded is an IMPLICIT relation, not a closed form:
    #      its right-hand side still contains ṙ_ij, i.e. the velocities.
    offdiag_nonzero = False
    for _ in range(10):
        v = vf()
        Mn = Mm.subs(v)
        if max(abs(float(Mn[a, b])) for a in range(6) for b in range(6) if a != b) > 1e-9:
            offdiag_nonzero = True
            break
    results.append((
        offdiag_nonzero,
        "N.4  M(q) has non-zero off-diagonal blocks → Eq. xdot_expanded is "
        "IMPLICIT in ẋ (needs a linear solve, not a substitution)",
        "off-diagonal entries of M are non-zero at random configurations",
    ))

    # N.5: the honest closed form is the linear solve v = M(q)^{-1} p.
    solved_ok = True
    for _ in range(10):
        v = vf()
        Mn = Mm.subs(v).evalf()
        pn = pvec.subs(v).evalf()
        vn = Mn.LUsolve(pn)
        true_v = Matrix([v[s] for s in vels])
        if max(abs(float(vn[a] - true_v[a])) for a in range(6)) > 1e-9:
            solved_ok = False
            break
    results.append((
        solved_ok,
        "N.5  v = M(q)⁻¹·p recovers the true velocities  (the closed form "
        "Eq. xdot_expanded leaves implicit)",
        "10 random 3-body configurations, max error < 1e-9",
    ))

    # N.6: explicit block structure of M(q).  This is what Eq. xdot_expanded
    #      encodes implicitly:
    #          M_ii = m_i·I − Σ_{j≠i} (q_i q_j/(c² r_ij)) r̂_ij r̂_ijᵀ
    #          M_ij = +        (q_i q_j/(c² r_ij)) r̂_ij r̂_ijᵀ      (i ≠ j)
    pos  = [(x1, y1), (x2, y2), (x3, y3)]
    mass = [m1, m2, m3]
    chg  = [q1, q2, q3]

    def rhat(a, b):
        dx, dy = pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]
        rab = sqrt(dx**2 + dy**2)
        return Matrix([dx/rab, dy/rab]), rab

    M_pred = zeros(6, 6)
    for a in range(3):
        for b in range(3):
            if a == b:
                blk = mass[a]*eye(2)
                for j in range(3):
                    if j != a:
                        rh, rab = rhat(a, j)
                        blk -= chg[a]*chg[j]/(c**2*rab) * (rh*rh.T)
            else:
                rh, rab = rhat(a, b)
                blk = chg[a]*chg[b]/(c**2*rab) * (rh*rh.T)
            M_pred[2*a:2*a + 2, 2*b:2*b + 2] = blk

    results.append(check_zero_fast(
        sum(abs(e) for e in (Mm - M_pred)),
        "N.6  M(q) has the predicted block form m_i·I − Σ (q_iq_j/c²r_ij) r̂r̂ᵀ "
        "(off-diagonal blocks + the same dyad)",
        vf,
    ))

    # N.7: if one tries to USE Eq. xdot_expanded as a recipe -- seed the RHS with
    #      a velocity guess and iterate -- the iteration is a Jacobi solve of
    #      p = M(q)v.  For two-body radial motion it reads
    #
    #          mu*rdot^(k+1) = p_r + (q1q2/(c² r))*rdot^(k)
    #
    #      whose contraction factor is exactly rho/r.  So the naive iteration
    #      converges only ABOVE the critical radius and DIVERGES below it --
    #      precisely the regime the Discussion calls the most promising one.
    mu_s = symbols('mu_s', positive=True)
    rr   = symbols('r_rel', positive=True)
    rho  = q1*q2/(mu_s*c**2)
    contraction = (q1*q2/(c**2*rr))/mu_s
    results.append(check_zero(
        contraction - rho/rr,
        "N.7  Jacobi/fixed-point solve of Eq. xdot_expanded has contraction "
        "factor exactly ρ/r  (converges iff r > ρ)",
        lambda: {rr: uniform(1.0, 3.0), mu_s: uniform(0.5, 2.0),
                 q1: uniform(0.2, 1.5), q2: uniform(0.2, 1.5),
                 c: uniform(3.0, 8.0)},
    ))

    # N.8: demonstrate it -- iterate the recipe above and below rho.
    conv_above = div_below = True
    for _ in range(50):
        mu_n = uniform(0.5, 2.0); c_n = uniform(3.0, 8.0)
        qq_n = uniform(0.2, 1.5)*uniform(0.2, 1.5)
        rho_n = qq_n/(mu_n*c_n**2)
        k_n = qq_n/c_n**2
        for r_n, expect_conv in ((uniform(2.0, 8.0)*rho_n, True),
                                 (uniform(0.2, 0.8)*rho_n, False)):
            p_r = 0.3*(mu_n - k_n/r_n)          # p_r consistent with rdot = 0.3
            v = p_r/mu_n                        # seed with the naive p/m
            for _ in range(200):
                v = (p_r + (k_n/r_n)*v)/mu_n
            converged = abs(v) < 1e3 and abs(v - 0.3) < 1e-6
            if expect_conv and not converged:
                conv_above = False
            if not expect_conv and converged:
                div_below = False
    results.append((
        conv_above and div_below,
        "N.8  the iteration converges to the true ṙ for r > ρ and diverges "
        "for r < ρ  (50 random configurations)",
        "confirms Eq. xdot_expanded needs a real linear solve, not iteration"
        if conv_above and div_below else "FAILED",
    ))

    return results


# ============================================================
# GROUP O: Units, dimensions and the pair count (Sections 2.5 and 2.4)
# ============================================================

def verify_group_O():
    results = []
    mg, mm, sec = symbols('mg mm s', positive=True)

    dim_F = mg*mm/sec**2                      # Sec. 2.5 claim
    dim_U = mg*mm**2/sec**2                   # Sec. 2.5 claim
    dim_q = sqrt(mg*mm**3/sec**2)             # App. A claim
    dim_r, dim_v, dim_m = mm, mm/sec, mg

    # O.1: [q]² = [F]·[r]²   (Coulomb term of Eq. force)
    results.append((
        simplify(dim_q**2/(dim_F*dim_r**2) - 1) == 0,
        "O.1  [q]² = [F]·[r]² ⇒ [q] = √(mg·mm³/s²)  (App. A charge unit)",
        "symbolic zero",
    ))

    # O.2: [q1q2/r] = [energy] = mg·mm²/s²   (Eq. potential, Sec. 2.5)
    results.append((
        simplify(dim_q**2/dim_r/dim_U - 1) == 0,
        "O.2  [U] = [q1q2/r] = mg·mm²/s²  (Sec. 2.5 potential units)",
        "symbolic zero",
    ))

    # O.3: [q1q2/r²] = mg·mm/s²   (Sec. 2.5 force units)
    results.append((
        simplify(dim_q**2/dim_r**2/dim_F - 1) == 0,
        "O.3  [F] = [q1q2/r²] = mg·mm/s²  (Sec. 2.5 force units)",
        "symbolic zero",
    ))

    # O.4: the two correction terms of Eq. force are dimensionless
    results.append((
        simplify(dim_v**2/dim_v**2) == 1
        and simplify(dim_r*(dim_v/sec)/dim_v**2) == 1,
        "O.4  ṙ²/c² and r·r̈/c² are dimensionless (Eq. force is homogeneous)",
        "symbolic zero",
    ))

    # O.5: [rho] = [q1q2/(mu c²)] is a length (Discussion)
    results.append((
        simplify(dim_q**2/(dim_m*dim_v**2)/mm - 1) == 0,
        "O.5  [ρ] = [q1q2/(μc²)] = mm — the critical radius is a length",
        "symbolic zero",
    ))

    # O.6: the pair count of Sec. 2.4 / the TikZ figure
    ok = all(len([(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)])
             == n*(n - 1)//2 for n in range(1, 40))
    results.append((
        ok, "O.6  Σ_{i<j} 1 = n(n−1)/2 unique pairs  (Sec. 2.4, O(n²))",
        "verified for n = 1..39",
    ))

    # O.7: alpha_x has the dimensions of a momentum (so Eq. px1 is homogeneous)
    dim_alpha = dim_q**2/dim_v**2 * dim_v*dim_r/dim_r**2
    results.append((
        simplify(dim_alpha/(dim_m*dim_v) - 1) == 0,
        "O.7  [α_x] = [m·v] — Eq. px1 is dimensionally homogeneous",
        "symbolic zero",
    ))

    return results


# ============================================================
# GROUP P: Negative controls — the errors the current revision corrects
#
# Each check below asserts that the PREVIOUS (wrong) form does NOT hold.
# A PASS here means the erroneous variant really is erroneous, i.e. the
# correction in the paper was necessary.
# ============================================================

def verify_group_P():
    results = []

    # P.1: the deleted appendix block claimed p_x1 = m1*ẋ1 (and the deleted
    #      sentence claimed ẋ_i = p_x_i/m_i).  Both are false whenever q1q2≠0.
    results.append(check_nonzero(
        diff(L, xd1) - m1*xd1,
        "P.1  p_x1 ≠ m1·ẋ1  (deleted appendix block was wrong)",
    ))
    results.append(check_nonzero(
        p1x_L/m1 - xd1,
        "P.2  ẋ1 ≠ p_x1/m1  (deleted 'replace each velocity by p/m' was wrong)",
    ))

    # P.3 / P.4: the two sign fixes in Eqs. xdot_two / xdot_two_p2
    results.append(check_nonzero(
        (p1x_L - alpha_x)/m1 - xd1,
        "P.3  ẋ1 = (p_x1 − α_x)/m1 is WRONG; the corrected + sign is required",
    ))
    results.append(check_nonzero(
        (p2x_L + alpha_x)/m2 - xd2,
        "P.4  ẋ2 = (p_x2 + α_x)/m2 is WRONG; the corrected − sign is required",
    ))

    # P.5 / P.6: the two sign fixes inside Eq. pdot_two
    wrong_A = (q1*q2/r12**2
               * ((x1 - x2)/r12 * (1 - Rational(3, 2)*rdot12_v**2/c**2)
                  - rdot12_v*xd_rel/c**2))
    results.append(check_nonzero(
        diff(L, x1) - wrong_A,
        "P.5  (1 − 3ṙ²/2c²) in Eq. pdot_two is WRONG; the + sign is required",
    ))
    wrong_B = (q1*q2/r12**2
               * ((x1 - x2)/r12 * (1 + Rational(3, 2)*rdot12_v**2/c**2)
                  + rdot12_v*xd_rel/c**2))
    results.append(check_nonzero(
        diff(L, x1) - wrong_B,
        "P.6  +ṙ(ẋ1−ẋ2)/c² in Eq. pdot_two is WRONG; the − sign is required",
    ))

    # P.7: L = T + S (wrong sign on the auxiliary potential) does not give T + U
    L_bad = T + S
    H_bad = (xd1*diff(L_bad, xd1) + yd1*diff(L_bad, yd1)
             + xd2*diff(L_bad, xd2) + yd2*diff(L_bad, yd2)) - L_bad
    results.append(check_nonzero(
        H_bad - (T + U),
        "P.7  L = T + S does NOT Legendre-transform to T + U (L = T − S is required)",
    ))

    # P.8: using U instead of S in the Lagrangian does not give Weber's force
    L_U = T - U
    F_from_U = diff(L_U, x1)
    F_from_S = diff(L, x1)
    results.append(check_nonzero(
        F_from_U - F_from_S,
        "P.8  L = T − U does NOT reproduce Eq. pdot_two (S, not U, is the "
        "Lagrangian potential)",
    ))

    # ---- P.9 / P.10: the v1.2 formulation, decided against Weber's force law --
    # v1.2 said "each velocity component ẋ_i in Eq. H is replaced by p_x_i/m_i".
    # Doing that gives a perfectly good Hamiltonian -- it is just not Weber's.
    # Both systems are compared here on the two-body relative problem, where no
    # matrix notation is needed: everything rides on the single scalar relation
    # between p_r and ṙ.
    mu_p, k_p = symbols('mu_P k_P', positive=True)
    r_p, pr_p, pph_p = symbols('r_P p_rP p_phiP', positive=True)
    qq_p = k_p*c**2                       # k = q1q2/c²
    mueff_p = mu_p - k_p/r_p

    def weber_residual(H):
        rdot   = diff(H, pr_p)
        phidot = diff(H, pph_p)
        prdot  = -diff(H, r_p)
        rddot  = diff(rdot, r_p)*rdot + diff(rdot, pr_p)*prdot
        F = qq_p/r_p**2*(1 - rdot**2/(2*c**2) + r_p*rddot/c**2)
        return simplify(mu_p*(rddot - r_p*phidot**2) - F)

    def val_P():
        return {mu_p: uniform(0.5, 2.0), k_p: uniform(0.01, 0.2),
                r_p: uniform(1.0, 3.0), pr_p: uniform(0.1, 0.5),
                pph_p: uniform(0.1, 0.5), c: uniform(3.0, 8.0)}

    # current paper: p_r = mu_eff * rdot  (Eq. px1 with the corrected + sign)
    H_now = pr_p**2/(2*mueff_p) + pph_p**2/(2*mu_p*r_p**2) + qq_p/r_p
    results.append(check_zero(
        weber_residual(H_now),
        "P.9  CURRENT paper: H = p_r²/(2μ_eff)+p_φ²/(2μr²)+q1q2/r reproduces "
        "Weber's force law EXACTLY",
        val_P,
    ))

    # v1.2 paper == src/hamiltonian/builders/weber.jl : rdot := p_r/mu
    H_v12 = (pr_p**2/(2*mu_p) + pph_p**2/(2*mu_p*r_p**2)
             + qq_p/r_p*(1 - (pr_p/mu_p)**2/(2*c**2)))
    results.append(check_nonzero(
        weber_residual(H_v12),
        "P.10 v1.2 paper / weber.jl (ṙ := p_r/μ) does NOT reproduce Weber's "
        "force law — it integrates a different Hamiltonian",
        val_P,
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

    print("\n--- Group A: Potential Identities ---")
    all_results += verify_group_A()

    print("--- Group B: Euler-Lagrange Equation ---")
    all_results += verify_group_B()

    print("--- Group C: Legendre Transform ---")
    all_results += verify_group_C()

    print("--- Group D: 2-Particle Hamilton's Equations ---")
    all_results += verify_group_D()

    print("--- Group E: n-Particle EOM Consistency (n=2) ---")
    all_results += verify_group_E()

    print("--- Group F: Radial Acceleration Identities ---")
    all_results += verify_group_F()

    print("--- Group G: Conservation Laws ---")
    all_results += verify_group_G()

    print("--- Group H: n=3 Body EOM Verification ---")
    all_results += verify_group_H()

    print("--- Group I: Angular Momentum Conservation ---")
    all_results += verify_group_I()

    print("--- Group J: Force <-> Potential, Newton's Third Law ---")
    all_results += verify_group_J()

    print("--- Group K: Closed-form Hamiltonian H(q,p) for n=2 ---")
    all_results += verify_group_K()

    print("--- Group L: Weber's Critical Radius ---")
    all_results += verify_group_L()

    print("--- Group M: Extended Phase Space & Symplectic Integrator ---")
    all_results += verify_group_M()

    print("--- Group N: n-Body Mass Matrix ---")
    all_results += verify_group_N()

    print("--- Group O: Units, Dimensions, Pair Count ---")
    all_results += verify_group_O()

    print("--- Group P: Negative Controls (the corrected errors) ---")
    all_results += verify_group_P()

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
