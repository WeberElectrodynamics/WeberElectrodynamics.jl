"""
01_ratio_tests.py
-----------------
Hypergeometric ratio test (A=B §3.2) applied to the pair-sum terms that
appear in the N-body Weber Hamilton's equations (paper Eqs. 211, 217).

A sequence f(k) is HYPERGEOMETRIC iff the ratio r(k) = f(k+1)/f(k)
is a rational function of k.
"""

from sympy import symbols, simplify, sin, pi, S

k, N = symbols("k N", positive=True, integer=True)
v, c = symbols("v c", positive=True)

SEP = "-" * 60


def ratio_test(name, f_expr, k_sym=k):
    ratio = simplify(f_expr.subs(k_sym, k_sym + 1) / f_expr)
    is_hyper = ratio.is_rational_function(k_sym)
    label = "HYPERGEOMETRIC  ✓" if is_hyper else "NOT hypergeometric  ✗"
    print(f"\n{name}")
    print(f"  f(k)          = {f_expr}")
    print(f"  f(k+1)/f(k)   = {ratio}")
    print(f"  Rational in k?  {is_hyper}   →  {label}")
    return ratio, is_hyper


print(SEP)
print("Hypergeometric Ratio Test — N-body Weber Hamilton's equations")
print(SEP)

print("""
Hamilton's equations for the N-body Weber Hamiltonian (paper Eqs. 211, 217):

  ẋᵢ   = (1/mᵢ)[ pᵢ − Σⱼ≠ᵢ (qᵢqⱼ/c²) ṙᵢⱼ(xᵢ−xⱼ)/rᵢⱼ² ]

  ṗₓᵢ  = Σⱼ≠ᵢ (qᵢqⱼ/rᵢⱼ²)[ (xᵢ−xⱼ)/rᵢⱼ (1 − 3ṙᵢⱼ²/2c²)
                            + ṙᵢⱼ(ẋᵢ−ẋⱼ)/c² ]

For N IDENTICAL particles (mass m, charge q) in a 1D chain with spacing d:
  rᵢⱼ = |i−j|·d,  re-index by k = |i−j|.

The double pair-sum collapses to a single sum over k = 1..N-1 where:
  • (N−k) pairs have relative distance k
  • each contributes a k-dependent term
""")

print(SEP)
print("Test 1 — Coulomb energy term: f_C(k) = (N−k)/k")
print("  (appears in total potential energy sum)")
f_C = (N - k) / k
ratio_test("Coulomb energy pair term  f_C(k) = (N−k)/k", f_C)

print(SEP)
print("Test 2 — Force term: f_F(k) = 1/k²")
print("  (appears in ṗ equation, static limit ṙᵢⱼ = 0)")
f_F = S.One / k**2
ratio_test("Force pair term  f_F(k) = 1/k²", f_F)

print(SEP)
print("Test 3 — Weber-corrected energy: f_W(k) = (N−k)/k · (1 − v²/2c²)")
print("  (uniform relative radial speed v for all pairs)")
f_W = (N - k) / k * (1 - v**2 / (2 * c**2))
ratio_test("Weber energy pair term  f_W(k) = (N−k)/k·(1−v²/2c²)", f_W)

print(SEP)
print("Test 4 — Ring configuration: f_ring(k) = csc(πk/N)")
print("  (N particles on circle radius R, pair distance 2R·sin(πk/N))")
f_ring = 1 / sin(pi * k / N)
ratio_test("Ring csc term  f_ring(k) = csc(πk/N)", f_ring)

print(SEP)
print("""
SUMMARY
-------
Term                  | f(k)                    | Hypergeometric?
Coulomb energy        | (N−k)/k                 | YES  ✓
Force (1/k²)          | 1/k²                    | YES  ✓
Weber correction      | (N−k)/k · (1−v²/2c²)   | YES  ✓
Ring csc              | csc(πk/N)               | NO   ✗

The 1D-chain pair sums in Hamilton's equations ARE hypergeometric in k.
The ring configuration is NOT → Gosper/Zeilberger do not apply there.
""")
