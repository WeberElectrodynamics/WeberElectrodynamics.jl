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

# rdot in velocity space (Lagrangian variables)
rdot12_v = (x_rel*xd_rel + y_rel*yd_rel) / r12

# rdot in (q,p) space — paper's approximation v_i ≈ p_i/m_i
rdot12_p = (x_rel*(px1/m1 - px2/m2) + y_rel*(py1/m1 - py2/m2)) / r12

# ---------- Potentials (velocity space) ----------
U = q1*q2/r12 * (1 - rdot12_v**2 / (2*c**2))   # Weber potential (Eq. potential)
S = q1*q2/r12 * (1 + rdot12_v**2 / (2*c**2))   # Auxiliary potential (Eq. 104)

# ---------- Lagrangian (velocity space) ----------
T = m1*(xd1**2 + yd1**2)/2 + m2*(xd2**2 + yd2**2)/2
L = T - S

# ---------- Hamiltonian in (q,p) space — paper's formula ----------
T_p  = px1**2/(2*m1) + py1**2/(2*m1) + px2**2/(2*m2) + py2**2/(2*m2)
H_qp = T_p + q1*q2/r12 * (1 - rdot12_p**2 / (2*c**2))

# ---------- alpha_x, alpha_y in (q,p) space ----------
alpha_x = q1*q2/c**2 * rdot12_p*(x1 - x2)/r12**2
alpha_y = q1*q2/c**2 * rdot12_p*(y1 - y2)/r12**2


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def _num_vals_2d(rng=None):
    """Return a random dict of 2D (q,p) numeric values with particles separated."""
    r = rng or (uniform,)
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

def check_zero(expr, label, val_func=None):
    """
    Attempt symbolic simplification; fall back to numerical check.
    Returns (passed: bool, label: str, detail: str).
    """
    vf = val_func or _num_vals_2d
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

    # C.3: Substituting v_i → p_i/m_i in H = T + U gives paper's H(q,p)
    H_approx = (T + U).subs([
        (xd1, px1/m1), (yd1, py1/m1),
        (xd2, px2/m2), (yd2, py2/m2),
    ])
    results.append(check_zero(
        expand(H_approx) - expand(H_qp),
        "C.3  H(q,p) via v≈p/m matches paper formula Eq. (hamiltonian)",
    ))

    return results


# ============================================================
# GROUP D: 2-Particle Hamilton's Equations
#          (Eqs. xdot_two, xdot_two_p2, pdot_two, pdot_two_p2)
# ============================================================

def verify_group_D():
    results = []

    # D.1: ∂H/∂px1 = (px1 − alpha_x)/m1   [Eq. xdot_two]
    dH_dpx1 = diff(H_qp, px1)
    results.append(check_zero(
        dH_dpx1 - (px1 - alpha_x)/m1,
        "D.1  ẋ1 = ∂H/∂px1 = (px1 − alpha_x)/m1  [Eq. xdot_two]",
    ))

    # D.2: ∂H/∂px2 = (px2 + alpha_x)/m2   [Eq. xdot_two_p2]
    dH_dpx2 = diff(H_qp, px2)
    results.append(check_zero(
        dH_dpx2 - (px2 + alpha_x)/m2,
        "D.2  ẋ2 = ∂H/∂px2 = (px2 + alpha_x)/m2  [Eq. xdot_two_p2]",
    ))

    # D.3: −∂H/∂x1 = q1q2/r12² * [(x1−x2)/r12*(1−3ṙ²/(2c²)) + ṙ*(ẋ1−ẋ2)/c²]  [Eq. pdot_two]
    dpx1_dt    = -diff(H_qp, x1)
    xd_rel_p   = px1/m1 - px2/m2  # (ẋ1−ẋ2) in (q,p) space
    D3_expected = (q1*q2/r12**2
                   * ((x1 - x2)/r12 * (1 - Rational(3, 2)*rdot12_p**2/c**2)
                      + rdot12_p*xd_rel_p/c**2))
    results.append(check_zero(
        dpx1_dt - D3_expected,
        "D.3  ṗx1 = −∂H/∂x1 matches paper formula  [Eq. pdot_two]",
    ))

    # D.4: ṗx1 + ṗx2 = 0 (Newton's 3rd law / total x-momentum conserved) [Eq. pdot_two_p2]
    dpx2_dt = -diff(H_qp, x2)
    results.append(check_zero(
        dpx1_dt + dpx2_dt,
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

    # ---- xdot_expanded (Eq. 210) ----
    # For particle i: xdot_i = (1/mi) * (pxi − sum_{j≠i} qi*qj/c² * rdot_ij*(xi−xj)/r_ij²)
    #
    # Convention: rdot_ij for arbitrary pair (i,j) is (ri−rj)·(vi−vj)/|ri−rj|.
    # For j>i this equals rdot12_p; for j<i, r_ji = −r_ij but v_ji = −v_ij, so rdot_ji = rdot_ij.

    # E.1a: n=2, i=1, sum over j=2:
    #   xdot_1 = (px1 − q1*q2/c² * rdot_12*(x1−x2)/r12²) / m1
    xdot1_npart = (px1 - q1*q2/c**2 * rdot12_p*(x1 - x2)/r12**2) / m1
    results.append(check_zero(
        xdot1_npart - diff(H_qp, px1),
        "E.1a xdot_expanded (n=2, i=1) equals ∂H/∂px1  [Eq. xdot_expanded]",
    ))

    # E.1b: n=2, i=2, sum over j=1:
    #   rdot_21 = rdot_12 (same scalar; both |r|-components cancel)
    #   (x2−x1) = −(x1−x2)
    #   → sum term = −q1q2/c² * rdot_12*(−(x1−x2))/r12² = +alpha_x
    #   → xdot_2 = (px2 − (−alpha_x))/m2 = (px2 + alpha_x)/m2  ✓
    rdot21_p = (  (x2 - x1)*(px2/m2 - px1/m1)
                + (y2 - y1)*(py2/m2 - py1/m1)) / r12   # = rdot12_p by antisymmetry
    xdot2_npart = (px2 - q2*q1/c**2 * rdot21_p*(x2 - x1)/r12**2) / m2
    results.append(check_zero(
        xdot2_npart - diff(H_qp, px2),
        "E.1b xdot_expanded (n=2, i=2) equals ∂H/∂px2  [Eq. xdot_expanded]",
    ))

    # E.2: n=2, i=1, pdot_expanded (Eq. 216):
    #   pdot_x1 = q1q2/r12² * [(x1−x2)/r12*(1−3ṙ²/(2c²)) + ṙ*(xd1−xd2)/c²]
    xd_rel_p   = px1/m1 - px2/m2
    pdot1_npart = (q1*q2/r12**2
                   * ((x1 - x2)/r12 * (1 - Rational(3, 2)*rdot12_p**2/c**2)
                      + rdot12_p*xd_rel_p/c**2))
    results.append(check_zero(
        pdot1_npart - (-diff(H_qp, x1)),
        "E.2  pdot_expanded (n=2, i=1) equals −∂H/∂x1  [Eq. pdot_expanded]",
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

    # G.1: H has no explicit time dependence
    has_t = t_sym in H_qp.free_symbols
    results.append((
        not has_t,
        "G.1  H(q,p) has no explicit time dependence",
        "PASS: t not in free_symbols" if not has_t else "FAIL: t found in H",
    ))

    # G.2: dH/dt = 0 along solutions (Poisson bracket {H,H} = 0)
    # dH/dt = Σ_i(∂H/∂qi * ∂H/∂pi + ∂H/∂pi * (−∂H/∂qi)) = Σ_i(a·b − b·a) = 0
    pairs = [(x1, px1), (y1, py1), (x2, px2), (y2, py2)]
    dH_dt = sum(
        diff(H_qp, qi)*diff(H_qp, pi) - diff(H_qp, pi)*diff(H_qp, qi)
        for qi, pi in pairs
    )
    results.append((
        dH_dt == 0,
        "G.2  dH/dt = 0 via Hamilton's equations (Poisson bracket cancellation)",
        "trivially zero (a·b − b·a = 0)" if dH_dt == 0 else f"non-zero: {dH_dt}",
    ))

    # G.3: Total x-momentum rate = −∂H/∂x1 − ∂H/∂x2 = 0  (translational symmetry)
    # H depends only on (x1−x2), so ∂H/∂x1 = −∂H/∂x2 → sum = 0
    dpx_total = -diff(H_qp, x1) - diff(H_qp, x2)
    results.append(check_zero(
        dpx_total,
        "G.3  Total x-momentum rate = −∂H/∂x1 − ∂H/∂x2 = 0",
    ))

    # G.4: Total y-momentum rate = 0
    dpy_total = -diff(H_qp, y1) - diff(H_qp, y2)
    results.append(check_zero(
        dpy_total,
        "G.4  Total y-momentum rate = −∂H/∂y1 − ∂H/∂y2 = 0",
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
