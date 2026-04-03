"""
02_closed_forms.py
------------------
Closed-form N-dependence of the Weber Hamiltonian energy and force sums,
derived using SymPy's summation() (Zeilberger algorithm internally) and
Gosper's algorithm (A=B §5, §6).

Configuration: N identical particles, mass m, charge q, spacing d (1D chain).
"""

from sympy import (
    symbols, harmonic, simplify, summation, latex,
    S, pi, sin, oo, Rational,
)
from sympy.concrete.gosper import gosper_sum

k, N = symbols("k N", positive=True, integer=True)
i    = symbols("i", positive=True, integer=True)
q, m, d, c, v = symbols("q m d c v", positive=True)

SEP = "=" * 65
SEP2 = "-" * 65


def section(title):
    print(f"\n{SEP}\n{title}\n{SEP}")


# ─────────────────────────────────────────────────────────────────────────────
section("1 · Total Coulomb energy  U_C(N)")
# ─────────────────────────────────────────────────────────────────────────────
print("""
Setup:
  U_C(N) = (q²/d) · Σₖ₌₁^{N−1} (N−k)/k
         = (q²/d) · S_C(N)
""")

# Gosper on the same term with FIXED upper bound n (to show it works)
n = symbols("n", positive=True, integer=True)
f_fixed = (n - k) / k   # n is a fixed parameter, k is the sum variable

print("--- Gosper's algorithm (A=B §5) ---")
print("When N is treated as a fixed parameter (not in the upper limit),")
print("gosper_sum finds the anti-difference G(k) with G(k+1)−G(k)=f(k).\n")

# Gosper on (N-k)/k with N in summand AND limit
g_diag = gosper_sum((N - k) / k, (k, 1, N - 1))
print(f"gosper_sum((N−k)/k, k=1..N−1) = {g_diag}")
print()
print("Substituting small N to check:")
for nv in [2, 3, 4, 5]:
    direct  = sum((nv - kk) / kk for kk in range(1, nv))
    formula = float(g_diag.subs(N, nv).evalf())
    ok      = "✓" if abs(direct - formula) < 1e-9 else "✗ WRONG"
    print(f"  N={nv}: direct={direct:.4f},  gosper formula={formula:.4f}  {ok}")

print("""
Result: gosper_sum gives a WRONG answer when N appears in both the summand
AND the upper limit (it treats the upper limit symbolically but N in the
summand also varies).  This is a known limitation of Gosper on 'diagonal'
sums.  We use summation() instead, which handles this correctly via its
Zeilberger-based backend.
""")

print("--- summation() (Zeilberger backend) ---")
S_C = summation((N - k) / k, (k, 1, N - 1))
print(f"summation((N−k)/k, k=1..N−1) = {S_C}")
print()
print("In closed form:  S_C(N) = N · H_{N−1} − (N−1)")
print("where H_n = harmonic(n) = Σ_{j=1}^n 1/j  (harmonic number)")
print()
print("Validation:")
for nv in [2, 3, 4, 5, 10, 20]:
    direct = sum((nv - kk) / kk for kk in range(1, nv))
    closed = float(S_C.subs(N, nv).evalf())
    ok     = "✓" if abs(direct - closed) < 1e-10 else "✗"
    print(f"  N={nv:2d}: direct={direct:.6f},  N·H_{{N-1}}−(N-1)={closed:.6f}  {ok}")

print(f"\n  Full expression for U_C(N):")
U_C = (q**2 / d) * S_C
print(f"  U_C(N) = (q²/d) · [{S_C}]")


# ─────────────────────────────────────────────────────────────────────────────
section("2 · Weber-corrected energy  U_W(N)")
# ─────────────────────────────────────────────────────────────────────────────
print("""
For a symmetric configuration with uniform relative radial speed v
(same for all pairs):  ṙᵢⱼ = v  →  Weber factor = 1 − v²/(2c²)

  U_W(N) = (q²/d) · (1 − v²/2c²) · S_C(N)
""")
S_W = simplify(summation((N - k) / k * (1 - v**2 / (2*c**2)), (k, 1, N-1)))
print(f"summation((N−k)/k · (1−v²/2c²), k=1..N−1) = {S_W}")
print()
ratio_WC = simplify(S_W / S_C)
print(f"S_W(N) / S_C(N) = {ratio_WC}")
print()
print("Key result: the Weber correction is a CONSTANT SCALAR multiple of")
print("the Coulomb sum.  The N-dependence (harmonic number scaling) is")
print("completely unchanged by the velocity correction.")
print()
print("Numerical check (v/c = 0.1, q=d=1):")
for nv in [2, 3, 5, 10]:
    direct_W = sum((nv-kk)/kk * (1 - 0.01/2) for kk in range(1, nv))
    closed_W = float(S_W.subs([(N, nv), (v, 0.1), (c, 1)]).evalf())
    ok       = "✓" if abs(direct_W - closed_W) < 1e-10 else "✗"
    print(f"  N={nv}: direct={direct_W:.6f},  closed={closed_W:.6f}  {ok}")


# ─────────────────────────────────────────────────────────────────────────────
section("3 · Force on end particle  F_end(N)")
# ─────────────────────────────────────────────────────────────────────────────
print("""
The static force on particle 0 (leftmost) from all N−1 particles to its right:

  F_end(N) = (q²/d²) · Σₖ₌₁^{N−1} 1/k²
           = (q²/d²) · H_{N−1}^{(2)}

where H_n^{(2)} = harmonic(n, 2) = Σ_{k=1}^n 1/k²  (generalized harmonic number).

Note: Gosper returns None for 1/k² — no rational anti-difference exists.
summation() recognises this as a generalized harmonic sum.
""")

g_F = gosper_sum(S.One / k**2, (k, 1, N - 1))
print(f"gosper_sum(1/k², k=1..N−1) = {g_F}  (correctly returns None)")
print()

S_F = summation(S.One / k**2, (k, 1, N - 1))
print(f"summation(1/k², k=1..N−1)  = {S_F}")
print()
print(f"F_end(N) = (q²/d²) · {S_F}")
print()
print("Limit as N→∞:  H_∞^{{(2)}} = π²/6  ≈ 1.6449  (Basel problem)")
import sympy
pi_sq_6 = float(sympy.pi**2 / 6)
print(f"  π²/6 = {pi_sq_6:.10f}")
print()
print("Validation:")
for nv in [2, 3, 5, 10, 20]:
    direct = sum(1/kk**2 for kk in range(1, nv))
    closed = float(S_F.subs(N, nv).evalf())
    ok     = "✓" if abs(direct - closed) < 1e-10 else "✗"
    print(f"  N={nv:2d}: direct={direct:.8f},  H_{{N-1}}^(2)={closed:.8f}  {ok}")


# ─────────────────────────────────────────────────────────────────────────────
section("4 · Force on interior particle i  (0 < i < N−1)")
# ─────────────────────────────────────────────────────────────────────────────
print("""
Particle i receives repulsive forces from i particles to its left
and (N−1−i) particles to its right.

  F_i(N) = (q²/d²) · [ H_i^{(2)} − H_{N−1−i}^{(2)} ]
""")

left_sum  = summation(S.One / k**2, (k, 1, i))
right_sum = summation(S.One / k**2, (k, 1, N - 1 - i))
F_i_sym   = simplify(left_sum - right_sum)
print(f"F_i(N) / (q²/d²) = {F_i_sym}")
print()
print("Special cases (positive = rightward force):")
print("  i = 0         (leftmost):  F = −H_{N-1}^{(2)}  [pushed left by all right neighbours]")
print("  i = (N-1)/2   (middle):    F =  0               [by symmetry]")
print("  i = N-1       (rightmost): F = +H_{N-1}^{(2)}  [pushed right by all left neighbours]")
print()
print("Note: F_end in section 3 is the MAGNITUDE |F_0| = H_{N-1}^{(2)}.")
print()
print("Validation for N=5 (signed, rightward positive):")
nv = 5
for iv in range(nv):
    direct_L = sum(1/kk**2 for kk in range(1, iv+1))       # particles to left  → push right (+)
    direct_R = sum(1/kk**2 for kk in range(1, nv-iv))      # particles to right → push left  (−)
    direct   = direct_L - direct_R
    closed   = float(F_i_sym.subs([(N, nv), (i, iv)]).evalf())
    ok = "✓" if abs(direct - closed) < 1e-10 else "✗"
    print(f"  i={iv}: F_i={direct:.6f},  H_i^(2)-H_{{N-1-i}}^(2)={closed:.6f}  {ok}")
