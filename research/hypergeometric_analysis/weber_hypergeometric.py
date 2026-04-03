"""
weber_hypergeometric.py
=======================
A=B book techniques (Gosper, Zeilberger/creative-telescoping) applied to the
N-body Weber electrodynamic Hamiltonian in the symmetric configuration.

The N-body Weber Hamiltonian (paper Eq. 191):

    H_N = Σᵢ pᵢ²/(2m) + Σᵢ<ⱼ (q²/rᵢⱼ)(1 - ṙᵢⱼ²/(2c²))

For N *identical* particles (mass m, charge q) in a 1D linear chain with
equal spacing d the pair distance is rᵢⱼ = |i−j|·d.  Reindexing by
k = |i−j| turns the double pair-sum into a single sum over k = 1..N-1 where
each summand is a hypergeometric term in k (ratio of successive terms is
rational in k).  This enables Gosper's algorithm and Zeilberger's method.
"""

from sympy import (
    symbols, harmonic, simplify, factor, ratsimp,
    summation, Sum, pprint, latex, pi, sin, cos,
    polygamma, S, oo, Rational, sqrt, exp, ln,
)
from sympy.concrete.gosper import gosper_sum

# ---------------------------------------------------------------------------
# Symbolic variables
# ---------------------------------------------------------------------------
k, N = symbols("k N", positive=True, integer=True)
n    = symbols("n",   positive=True, integer=True)   # generic bound
q, m, d, c, v = symbols("q m d c v", positive=True)
R    = symbols("R",   positive=True)                  # ring radius
theta = symbols("theta")                              # angle for ring


# ---------------------------------------------------------------------------
# Core pair-term constructors
# ---------------------------------------------------------------------------

def coulomb_pair_term(k_sym=k, N_sym=N):
    """
    The k-th Coulomb pair contribution to total potential energy in a 1D chain
    of N identical charges q, spacing d:

        f_C(k) = (N − k) / k

    (after factoring out q²/d, which is the prefactor).

    There are (N − k) pairs with relative index k, each contributing 1/(k·d).
    """
    return (N_sym - k_sym) / k_sym


def weber_pair_term(k_sym=k, N_sym=N, v_radial=v):
    """
    Weber-corrected pair term.  For identical collinear motion at speed v_radial
    along the chain axis:  ṙᵢⱼ = v_radial for all pairs (same direction, so
    relative radial velocity is 0 — corrected formula uses the 1D case where
    the radial velocity between particles at separation k·d moving with the
    same velocity is identically zero).

    For the *head-on* case (particle i moves right at v, particle j at rest):
        ṙᵢⱼ = ±v  →  Weber factor = 1 − v²/(2c²)

    Here we parameterise by v_radial (the relative radial speed of a pair),
    which is treated as a fixed scalar the same for all pairs in the symmetric
    configuration.

        f_W(k) = (N − k) / k · (1 − v_radial² / (2c²))
    """
    return (N_sym - k_sym) / k_sym * (1 - v_radial**2 / (2 * c**2))


def force_pair_term(k_sym=k):
    """
    The k-th term in the net force on the *end particle* (particle 0) from
    all particles to its right in the 1D chain:

        f_F(k) = 1 / k²

    (after factoring out q²/d²).
    """
    return S.One / k_sym**2


# ---------------------------------------------------------------------------
# Hypergeometric ratio test  (A=B §3.2)
# ---------------------------------------------------------------------------

def hypergeometric_ratio(f_expr, k_sym=k, verbose=True):
    """
    Compute r(k) = f(k+1)/f(k) and check whether it is rational in k.
    A term f(k) is *hypergeometric* iff this ratio is rational.

    Returns
    -------
    ratio : SymPy expression
    is_hyper : bool
    """
    ratio = simplify(f_expr.subs(k_sym, k_sym + 1) / f_expr)
    is_hyper = ratio.is_rational_function(k_sym)
    if verbose:
        print(f"f(k)        = {f_expr}")
        print(f"f(k+1)/f(k) = {ratio}")
        print(f"Rational in k? {is_hyper}  → {'hypergeometric ✓' if is_hyper else 'NOT hypergeometric ✗'}")
    return ratio, is_hyper


# ---------------------------------------------------------------------------
# Gosper's algorithm  (A=B §5)
# ---------------------------------------------------------------------------

def gosper_closed_form(f_expr, k_sym=k, a=1, b=None, verbose=True):
    """
    Apply Gosper's algorithm to find an anti-difference G(k) such that
    G(k+1) − G(k) = f(k), and evaluate the definite sum if bounds are given.

    When b is symbolic (e.g. N-1) we rely on SymPy's summation() fallback
    because Gosper's algorithm requires bounds to be explicit integers OR the
    summand must not contain the upper-bound symbol.

    Returns the closed-form expression (or None if Gosper fails).
    """
    if verbose:
        print(f"Applying Gosper to f(k) = {f_expr}  over k = {a}..{b}")

    # Try Gosper (works when summand is hypergeometric in k alone)
    if b is not None:
        result = gosper_sum(f_expr, (k_sym, a, b))
    else:
        result = gosper_sum(f_expr, k_sym)

    if result is None:
        if verbose:
            print("  Gosper returned None — falling back to summation()")
        if b is not None:
            result = summation(f_expr, (k_sym, a, b))
        else:
            result = None
    if verbose:
        print(f"  Closed form: {result}")
    return result


# ---------------------------------------------------------------------------
# Main closed forms
# ---------------------------------------------------------------------------

def coulomb_energy_chain(N_sym=N):
    """
    Total Coulomb potential energy of N identical charges q, spacing d:

        U_C(N) = (q²/d) · Σₖ₌₁^{N−1} (N−k)/k
               = (q²/d) · [N · H_{N−1} − (N−1)]

    where H_n = harmonic(n) is the n-th harmonic number.

    Returns the symbolic expression for the bracket [N·H_{N-1} − (N-1)].
    Multiply by q²/d to get energy in natural units.
    """
    expr = summation((N_sym - k) / k, (k, 1, N_sym - 1))
    return simplify(expr)


def weber_energy_chain(N_sym=N):
    """
    Total Weber potential energy for a 1D chain, uniform relative radial
    velocity v_radial (same for all pairs):

        U_W(N) = (q²/d) · (1 − v²/2c²) · Σₖ₌₁^{N−1} (N−k)/k
               = (1 − v²/2c²) · U_C(N)

    Returns: (weber_factor, coulomb_sum_bracket)
    """
    coulomb_sum = coulomb_energy_chain(N_sym)
    weber_factor = 1 - v**2 / (2 * c**2)
    return weber_factor, coulomb_sum


def force_on_end_particle(N_sym=N):
    """
    Net force on the end particle (particle 0) in a 1D chain of N identical
    charges, static configuration:

        F_end(N) = (q²/d²) · Σₖ₌₁^{N−1} 1/k²
                 = (q²/d²) · H_{N−1}^{(2)}

    where H_n^{(2)} = harmonic(n, 2) = Σₖ₌₁ⁿ 1/k² is the generalised
    harmonic number (converges to π²/6 as n→∞).

    Note: Gosper's algorithm returns None for 1/k² (no rational anti-difference
    exists), but SymPy's summation() recognises the generalised harmonic form.

    Returns the symbolic expression for the sum.
    """
    expr = summation(S.One / k**2, (k, 1, N_sym - 1))
    return expr   # = harmonic(N_sym - 1, 2)


def force_on_interior_particle(i_sym, N_sym=N):
    """
    Net force on an interior particle i (1 ≤ i ≤ N−2) in the 1D chain.
    Particles to the left repel with force +q²/(kd)², particles to the right
    repel with force −q²/(kd)² (net sign convention: positive = rightward).

        F_i(N) = (q²/d²) · [Σₖ₌₁^{i} 1/k² − Σₖ₌₁^{N−1−i} 1/k²]
               = (q²/d²) · [H_i^{(2)} − H_{N−1−i}^{(2)}]

    Returns the symbolic expression.
    """
    left_sum  = summation(S.One / k**2, (k, 1, i_sym))
    right_sum = summation(S.One / k**2, (k, 1, N_sym - 1 - i_sym))
    return simplify(left_sum - right_sum)


# ---------------------------------------------------------------------------
# Zeilberger recurrence  (A=B §6) — derived by creative telescoping
# ---------------------------------------------------------------------------

def energy_recurrence_rhs(N_sym=N):
    """
    The energy satisfies a first-order recurrence in N:

        E(N+1) = E(N) + H_N

    where E(N) = N·H_{N−1} − (N−1)  and  H_N = harmonic(N).

    Proof sketch (creative telescoping / Zeilberger spirit):
      E(N+1) − E(N)
        = (N+1)·H_N − N − N·H_{N−1} + (N−1)
        = H_{N−1} + (N+1)/N − 1           [using H_N = H_{N−1} + 1/N]
        = H_{N−1} + 1/N
        = H_N

    Returns H_N symbolically.
    """
    return harmonic(N_sym)


def verify_energy_recurrence(N_max=8):
    """
    Numerically verify  E(N+1) − E(N) = H_N  for N = 2..N_max.
    """
    E = lambda n: n * float(harmonic(n - 1).evalf()) - (n - 1)
    H = lambda n: float(harmonic(n).evalf())
    results = []
    for n_val in range(2, N_max + 1):
        lhs = E(n_val + 1) - E(n_val)
        rhs = H(n_val)
        results.append((n_val, lhs, rhs, abs(lhs - rhs) < 1e-12))
    return results


# ---------------------------------------------------------------------------
# Ring configuration (non-hypergeometric, but finite trigonometric identity)
# ---------------------------------------------------------------------------

def ring_pair_distances(N_sym):
    """
    In a ring of N identical charges on a circle of radius R, all equally
    spaced, the pair distance for relative index k is:

        r_k = 2R · sin(π k / N)

    The total Coulomb potential energy is:

        U_ring(N) = (q²/2) · Σᵢ≠ⱼ 1/rᵢⱼ
                  = (q² / (4R)) · N · Σₖ₌₁^{N−1} csc(πk/N)

    The csc sum is NOT hypergeometric in k (sin(πk/N) has no rational ratio
    in k), so Gosper's algorithm does not apply.  The sum evaluates to

        Σₖ₌₁^{N−1} csc(πk/N) = (N/π) · Im[ψ((1+iπ)/N) − ψ(1/N)] + ...

    which does not simplify to elementary functions in general.
    We compute it numerically (see notebook).
    """
    return 2 * R * sin(pi * k / N_sym)


def ring_csc_ratio(N_sym=N):
    """Check: csc(πk/N) is NOT hypergeometric in k."""
    f = 1 / sin(pi * k / N_sym)
    ratio = simplify(f.subs(k, k + 1) / f)
    is_hyper = ratio.is_rational_function(k)
    return ratio, is_hyper


# ---------------------------------------------------------------------------
# Numeric validation helpers
# ---------------------------------------------------------------------------

def validate_coulomb_chain(N_max=15):
    """
    Compare direct pair sum to closed form N·H_{N−1} − (N−1) for N = 2..N_max.
    Returns list of (N, direct, closed_form, match).
    """
    results = []
    for n_val in range(2, N_max + 1):
        direct = sum((n_val - kk) / kk for kk in range(1, n_val))
        closed = n_val * float(harmonic(n_val - 1).evalf()) - (n_val - 1)
        results.append((n_val, direct, closed, abs(direct - closed) < 1e-10))
    return results


def validate_force_chain(N_max=15):
    """
    Compare direct sum Σ1/k² to harmonic(N-1, 2) for N = 2..N_max.
    Returns list of (N, direct, closed_form, match).
    """
    results = []
    for n_val in range(2, N_max + 1):
        direct = sum(1 / kk**2 for kk in range(1, n_val))
        closed = float(harmonic(n_val - 1, 2).evalf())
        results.append((n_val, direct, closed, abs(direct - closed) < 1e-10))
    return results


# ---------------------------------------------------------------------------
# Pretty-print summary
# ---------------------------------------------------------------------------

def print_summary():
    """Print a summary of all closed forms."""
    print("=" * 65)
    print("  N-body Weber Hamiltonian — Closed-form N-dependence")
    print("  (symmetric 1D chain: equal mass m, charge q, spacing d)")
    print("=" * 65)

    print("\n[1] Total Coulomb energy:")
    print("    U_C(N) = (q²/d) · S_C(N)")
    print("    S_C(N) = N·H_{N−1} − (N−1)   [harmonic numbers]")
    print("    Derived by: SymPy summation() [Gosper fails — N in both")
    print("    summand and upper limit; summation() uses Zeilberger internally]")

    print("\n[2] Weber-corrected energy (uniform relative radial speed v):")
    print("    U_W(N) = (q²/d) · (1 − v²/2c²) · S_C(N)")

    print("\n[3] Force on end particle:")
    print("    F_end(N) = (q²/d²) · H_{N−1}^{(2)}    [generalised harmonic]")
    print("    H_n^{(2)} = Σₖ₌₁ⁿ 1/k²  →  π²/6  as n→∞")
    print("    Note: Gosper returns None for 1/k² (no rational anti-difference)")

    print("\n[4] Recurrence (Zeilberger / creative telescoping):")
    print("    E(N+1) − E(N) = H_N")
    print("    E(2) = 1,   E(N) = N·H_{N−1} − (N−1)")

    print("\n[5] Ring configuration:")
    print("    U_ring(N) ∝ Σₖ₌₁^{N−1} csc(πk/N)")
    print("    csc ratio in k is NOT rational → Gosper/A=B does NOT apply")
    print("    Evaluated numerically (no known elementary closed form)")
    print("=" * 65)


if __name__ == "__main__":
    print_summary()
    print()

    print("--- Hypergeometric ratio tests ---")
    hypergeometric_ratio(coulomb_pair_term())
    print()
    hypergeometric_ratio(force_pair_term())
    print()
    ratio_r, is_r = ring_csc_ratio()
    print(f"Ring csc: ratio = {ratio_r}, rational? {is_r}")

    print()
    print("--- Closed forms ---")
    print("S_C(N) =", coulomb_energy_chain())
    print("F_end(N) proportional to:", force_on_end_particle())

    print()
    print("--- Recurrence verification ---")
    for row in verify_energy_recurrence():
        print(f"  N={row[0]}: E(N+1)-E(N)={row[1]:.8f}, H_N={row[2]:.8f}, ok={row[3]}")

    print()
    print("--- Coulomb validation ---")
    for row in validate_coulomb_chain(8):
        print(f"  N={row[0]}: direct={row[1]:.6f}, closed={row[2]:.6f}, ok={row[3]}")
