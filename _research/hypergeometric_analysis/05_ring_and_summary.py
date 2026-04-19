"""
05_ring_and_summary.py
----------------------
Ring configuration analysis (non-hypergeometric) and full summary table
of A=B results for the N-body Weber Hamiltonian.
"""

from sympy import symbols, harmonic, sin, pi, simplify, Rational, S, N as N_evalf
import math

k, N = symbols("k N", positive=True, integer=True)
q, d, c, v, R = symbols("q d c v R", positive=True)

SEP  = "=" * 65
SEP2 = "-" * 65


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("Ring Configuration: N charges equally spaced on circle radius R")
print(SEP)

print("""
Pair distance for relative index k:  r_k = 2R · sin(πk/N)

Total Coulomb energy:
  U_ring(N) = (Nq²/4R) · Σₖ₌₁^{N−1} csc(πk/N)

The csc sum is NOT hypergeometric in k → Gosper/Zeilberger do NOT apply.
""")

# Ratio test (already done in 01, just print result)
f_ring = 1 / sin(pi * k / N)
ratio_ring = simplify(f_ring.subs(k, k + 1) / f_ring)
print(f"Ratio f_ring(k+1)/f_ring(k) = {ratio_ring}")
print(f"Rational in k?  {ratio_ring.is_rational_function(k)}   → NOT hypergeometric\n")

# Numerical values of the ring sum
def ring_csc_sum(n_val):
    return sum(1.0 / math.sin(math.pi * kk / n_val) for kk in range(1, n_val))

print("Numerical values of  T(N) = Σₖ₌₁^{N−1} csc(πk/N):")
print(f"{'N':>4} | {'T(N)':>12} | {'T(N)/(N/π·lnN)':>16} | {'U_ring·4R/Nq²':>14}")
print(SEP2)
for nv in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]:
    t = ring_csc_sum(nv)
    if nv > 1:
        asymp = (nv / math.pi) * math.log(nv) if nv > 1 else 0
        ratio_a = t / asymp if asymp > 0 else float('nan')
    else:
        ratio_a = float('nan')
    print(f"{nv:>4} | {t:>12.6f} | {ratio_a:>16.6f} | {t:>14.6f}")

print(f"""
Large-N asymptotic:  T(N) ~ (N/π) · ln(N)
This is NOT a closed form — only an asymptotic expansion.

Special cases with known exact values:
  N=2: T(2) = csc(π/2) = 1
  N=4: T(4) = csc(π/4) + csc(π/2) + csc(3π/4) = 1 + 2√2 ≈ {1 + 2*math.sqrt(2):.6f}
  N=6: T(6) = 2·csc(π/6) + csc(π/3) + csc(π/2) + csc(2π/3)
            = 4 + 2/√3 + 1 + 2/√3 = 5 + 4/√3 ≈ {5 + 4/math.sqrt(3):.6f}

These involve algebraic numbers (√2, √3) but no general elementary formula.
""")


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("Full N-body Chain: Energy and force table  (q=d=1)")
print(SEP)

def H(n_val):
    """Harmonic number H_n = Σ 1/k, k=1..n"""
    return float(harmonic(n_val).evalf())

def H2(n_val):
    """Generalized harmonic H_n^{(2)} = Σ 1/k², k=1..n"""
    return float(harmonic(n_val, 2).evalf())

def S_C(n_val):
    """N·H_{N-1} − (N-1)"""
    return n_val * H(n_val - 1) - (n_val - 1)

print(f"\n{'N':>4} | {'U_C(N)':>12} | {'U_C/N':>10} | {'F_end(N)':>12} | {'F_end/π²·6':>12}")
print(SEP2)
pi2_6 = math.pi**2 / 6
for nv in range(2, 21):
    uc     = S_C(nv)
    f_end  = H2(nv - 1)
    print(f"{nv:>4} | {uc:>12.6f} | {uc/nv:>10.6f} | {f_end:>12.8f} | {f_end/pi2_6:>12.6f}")

print(f"\n  U_C(N) / N  →  ln(N) + γ  (Euler-Mascheroni γ ≈ 0.5772)")
print(f"  F_end(N)    →  π²/6 ≈ {pi2_6:.8f}  (Basel problem)")


# ─────────────────────────────────────────────────────────────────────────────
print()
print(SEP)
print("COMPLETE SUMMARY OF A=B RESULTS")
print(SEP)
print("""
Configuration: N identical particles, mass m, charge q, spacing d (1D chain).
All results are exact closed forms in terms of harmonic numbers.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[1]  TOTAL COULOMB ENERGY
     U_C(N) = (q²/d) · [N · H_{N−1} − (N−1)]
     Method: summation() with Zeilberger backend
     Growth: O(N ln N)

[2]  WEBER-CORRECTED ENERGY (uniform radial speed v for all pairs)
     U_W(N) = (q²/d) · (1 − v²/2c²) · [N · H_{N−1} − (N−1)]
     = (1 − v²/2c²) · U_C(N)
     Method: scalar multiplication — same A=B structure

[3]  FORCE ON END PARTICLE (static, particle 0)
     F_end(N) = (q²/d²) · H_{N−1}^{(2)}
     Limit: F_end → (q²/d²) · π²/6  as N → ∞  (Basel)
     Note: Gosper returns None for 1/k²;
           summation() gives generalized harmonic number.

[4]  FORCE ON INTERIOR PARTICLE i  (0 < i < N−1)
     F_i(N)   = (q²/d²) · [H_i^{(2)} − H_{N−1−i}^{(2)}]
     F_{(N-1)/2} = 0  (middle particle, by symmetry)

[5]  ENERGY RECURRENCE (Zeilberger certificate)
     E(N+1) − E(N) = H_N,   E(2) = 1
     where E(N) = N · H_{N−1} − (N−1)
     Physical: adding the N-th particle increases energy by H_N.

[6]  RING CONFIGURATION
     U_ring(N) ∝ Σ_{k=1}^{N-1} csc(πk/N)   (NOT hypergeometric)
     ~ (N/π) · ln N  for large N
     No Gosper/Zeilberger closed form exists.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

NOTATION
  H_n       = harmonic(n)    = Σ_{k=1}^n 1/k
  H_n^{(2)} = harmonic(n, 2) = Σ_{k=1}^n 1/k²
  γ         ≈ 0.5772157  (Euler-Mascheroni constant)
  H_n       ~ ln(n) + γ    for large n

SCOPE NOTE
  These are EXACT N-scaling formulas for energy and force at fixed
  symmetric configurations — NOT solutions to the full nonlinear dynamics.
  The general N-body problem has no closed-form solution (Poincaré-Bruns).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")
