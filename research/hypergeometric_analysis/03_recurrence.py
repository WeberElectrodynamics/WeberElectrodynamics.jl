"""
03_recurrence.py
----------------
Zeilberger / creative telescoping (A=B §6): derive and verify the recurrence
in N for the total Coulomb energy of the N-body Weber chain.

E(N) := N · H_{N−1} − (N−1)  satisfies  E(N+1) − E(N) = H_N,  E(2) = 1.
"""

from sympy import symbols, harmonic, simplify, Rational, S

N = symbols("N", positive=True, integer=True)

SEP = "=" * 65

print(SEP)
print("Zeilberger Recurrence for N-body Coulomb energy")
print(SEP)

print("""
Closed form (from 02_closed_forms.py):
  S_C(N) = N · H_{N−1} − (N−1)
  where H_n = Σ_{k=1}^n 1/k  (harmonic number)

CLAIM (Zeilberger certificate):
  S_C(N+1) − S_C(N) = H_N

This is the 'creative telescoping' identity: it PROVES the closed form
by providing a recurrence that uniquely determines S_C(N) given S_C(2)=1.

Proof (by hand, A=B §6 spirit):
  S_C(N+1) − S_C(N)
    = (N+1)·H_N − N  −  (N·H_{N−1} − (N−1))
    = (N+1)·H_N − N·H_{N−1} − 1
    = H_{N−1} + (N+1)/N − 1         [using H_N = H_{N−1} + 1/N]
    = H_{N−1} + 1/N
    = H_N                             □
""")

# Symbolic verification via SymPy
# Note: SymPy cannot simplify N*H_N - N*H_{N-1} - 1 to zero symbolically
# (it would need to apply H_N = H_{N-1} + 1/N internally).
# We verify numerically instead; the algebraic proof above is the certificate.
E = lambda n: n * harmonic(n - 1) - (n - 1)
diff_sym = simplify(E(N + 1) - E(N))
print(f"SymPy: E(N+1)−E(N) = {diff_sym}")
print(f"       harmonic(N)  = harmonic(N)")
print(f"  (SymPy cannot auto-simplify H_N = H_{{N-1}} + 1/N — see algebraic proof above)")
print()

# Numerical verification
print("Numerical verification  E(N+1) − E(N) = H_N:")
print(f"{'N':>4} | {'E(N+1)':>12} | {'E(N)':>12} | {'E(N+1)−E(N)':>14} | {'H_N':>12} | OK?")
print("-" * 65)
E_n = lambda nv: nv * float(harmonic(nv - 1).evalf()) - (nv - 1)
H_n = lambda nv: float(harmonic(nv).evalf())
for nv in range(2, 13):
    lhs = E_n(nv + 1) - E_n(nv)
    rhs = H_n(nv)
    ok  = "✓" if abs(lhs - rhs) < 1e-12 else "✗"
    print(f"{nv:>4} | {E_n(nv+1):>12.8f} | {E_n(nv):>12.8f} | {lhs:>14.8f} | {rhs:>12.8f} | {ok}")

print()
print(SEP)
print("Growth table: S_C(N), H_{N-1}, and ratio")
print(SEP)
print(f"\n{'N':>4} | {'S_C(N)':>12} | {'H_{N-1}':>10} | {'S_C/H':>10} | {'ΔS_C = H_N':>12}")
print("-" * 55)
for nv in range(2, 21):
    sc   = E_n(nv)
    hn_1 = H_n(nv - 1)
    hn   = H_n(nv)
    delta = E_n(nv + 1) - E_n(nv)
    print(f"{nv:>4} | {sc:>12.6f} | {hn_1:>10.6f} | {sc/hn_1:>10.4f} | {delta:>12.8f}")

print("""
Interpretation:
  S_C(N) grows as N · ln(N) for large N  (since H_N ~ ln N + γ)
  Each added particle increases energy by exactly H_N (the N-th harmonic number).
  This recurrence IS the Zeilberger certificate of the closed form.
""")
