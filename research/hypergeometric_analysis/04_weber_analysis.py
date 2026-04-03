"""
04_weber_analysis.py
--------------------
Detailed analysis of the Weber velocity-correction term in the N-body
Hamilton's equations.

Weber potential between a pair (i, j):
  U_ij = (q²/r_ij) · (1 − ṙ_ij²/2c²)

For the symmetric 1D chain of N identical particles, three physically
distinct cases are examined.
"""

from sympy import symbols, harmonic, simplify, summation, Rational, S, sqrt

k, N = symbols("k N", positive=True, integer=True)
v, c, q, d, m = symbols("v c q d m", positive=True)

SEP = "=" * 65

print(SEP)
print("Weber Velocity-Correction Analysis — N-body 1D Chain")
print(SEP)

print("""
Weber pair potential:  U_ij = (q²/r_ij) · (1 − ṙ_ij²/2c²)
where ṙ_ij = (r⃗_ij · v⃗_ij)/r_ij  is the radial component of relative velocity.

For the 1D chain:  r_ij = |i−j|·d,  ṙ_ij = (sign of approach)·|v_i − v_j|·(along chain)
""")

# Base Coulomb sum
S_C = summation((N - k) / k, (k, 1, N - 1))  # = N*H_{N-1} - (N-1)

print(SEP)
print("Case A: Collinear uniform motion (all particles same velocity v₀)")
print(SEP)
print("""
All particles move at the same velocity v₀ along the chain axis.
Relative velocity of any pair: v_rel = v₀ − v₀ = 0
Radial relative velocity: ṙ_ij = 0  for ALL pairs
Weber factor = 1 − 0²/(2c²) = 1

  U_W^A(N) = U_C(N)  [pure Coulomb — Weber correction is zero]

Physical meaning: uniform translation is a symmetry; Weber's law reduces to
Coulomb's for systems in uniform motion.
""")

# Momentum equation (Eq. 211) at rest (p_i = 0 initially, dot r = 0)
print("Position equation (Eq. 211) with ṙ_ij = 0:")
print("  ẋ_i = p_{x_i}/m  (velocity-coupling sum vanishes)")
print()
print("Momentum equation (Eq. 217) with ṙ_ij = 0:")
print("  ṗ_{x_i} = Σⱼ≠ᵢ (q²/r_ij²) · (x_i−x_j)/r_ij")
print("          = (q²/d²) · Σₖ₌₁^{N-1−i} 1/k² − (q²/d²) · Σₖ₌₁^i 1/k²")
print("  Static force formula: F_i = (q²/d²)[H_i^{(2)} − H_{N−1−i}^{(2)}]
  (left neighbours push right +, right neighbours push left −)")

print()
print(SEP)
print("Case B: Head-on symmetric approach (outer pair moving inward at speed v)")
print(SEP)
print("""
Particle 0 moves right at +v, particle N−1 moves left at −v, rest stationary.
For particle 0 vs particle k:  v_rel = v − 0 = v  (along chain axis)
Radial relative velocity: ṙ_{0k} = v  (particle 0 approaching particle k)
Weber factor for pair (0, k) = 1 − v²/(2c²)

For ALL OTHER pairs (i,j with i,j ≥ 1): ṙ_ij = 0 (both stationary).

Force on particle 0 (from all k = 1..N−1):
  ṗ_{x_0} = (q²/d²) · Σₖ₌₁^{N−1} (1/k²) · [1 − 3v²/(2c²) + v·0/c²]
           = (q²/d²) · (1 − 3v²/2c²) · H_{N−1}^{(2)}
""")
S_F = summation(S.One / k**2, (k, 1, N - 1))
# Weber factor for ṗ from Eq. 217 with ṙ_ij = v:
# (1 - 3*v^2/(2*c^2)) factor on the positional term
weber_force_factor = 1 - Rational(3, 2) * v**2 / c**2
print(f"  ṗ_{{x_0}} = (q²/d²) · (1 − 3v²/2c²) · H_{{N−1}}^{{(2)}}")
print()
print("Numerical check (v/c = 0.1, q=d=1):")
for nv in [3, 5, 10]:
    hf2 = float(S_F.subs(N, nv).evalf())
    # force factor: 1 - 3*(0.1)^2/2 = 1 - 0.015 = 0.985
    factor_val = 1 - 3*(0.1)**2/2
    print(f"  N={nv:2d}: F₀ = {factor_val:.4f} × {hf2:.6f} = {factor_val*hf2:.6f}")

print()
print(SEP)
print("Case C: Uniform relative radial speed v (all pairs, symmetric IC)")
print(SEP)
print("""
For configurations with a fixed relative radial speed v the same for all pairs,
the energy sum factors cleanly:

  U_W(N) = (q²/d) · (1 − v²/2c²) · Σₖ₌₁^{N−1} (N−k)/k

The Weber factor is a CONSTANT SCALAR — does not depend on k or N.
""")
S_W = simplify(summation((N - k) / k * (1 - v**2 / (2*c**2)), (k, 1, N-1)))
ratio_WC = simplify(S_W / S_C)
print(f"summation result: {S_W}")
print(f"S_W(N) / S_C(N) = {ratio_WC}")
print()

print("Weber correction fraction  (ΔU_Weber/U_Coulomb) = v²/2c²:")
print(f"{'v/c':>6} | {'ΔU/U_C':>10} | {'U_W/U_C':>10}")
print("-" * 32)
for beta in [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5]:
    correction = beta**2 / 2
    ratio_num  = 1 - correction
    print(f"{beta:>6.3f} | {correction:>10.6f} | {ratio_num:>10.6f}")

print()
print(SEP)
print("Key conclusion")
print(SEP)
print("""
For the symmetric N-body Weber chain:

1. STATIC configuration (all velocities zero):
   U_W = U_C,  F_i = (q²/d²)[H_i^{(2)} − H_{N−1−i}^{(2)}]
   A=B closed forms apply fully.

2. UNIFORM RELATIVE SPEED (Case C):
   U_W = (1 − v²/2c²) · U_C
   The Weber correction is N-INDEPENDENT — a constant scalar on U_C.
   The harmonic-number N-scaling is PRESERVED.

3. CASE-B (approach, Eq. 217):
   ṗ_{x_0} has (1 − 3v²/2c²) factor instead of (1 − v²/2c²).
   Still a constant scalar: N-dependence unchanged.

The A=B hypergeometric machinery gives EXACT closed forms for the
N-dependence in all three cases.  Weber's velocity term modifies only
the prefactor, never the harmonic-number structure.
""")
