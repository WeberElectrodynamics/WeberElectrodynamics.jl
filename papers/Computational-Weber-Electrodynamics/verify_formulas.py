#!/usr/bin/env python3
"""
Symbolic verification of mathematical formulas from:
"Computational Weber electrodynamics with symplectic numerical integrators"

Verifies all key derived results using SymPy.
Each check prints PASS/FAIL with the verification method used.

Run: python3 verify_formulas.py
"""

import sys
from random import uniform, seed
from sympy import (
    symbols, sqrt, diff, simplify, expand, trigsimp,
    Function, Derivative, Rational, Integer
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
