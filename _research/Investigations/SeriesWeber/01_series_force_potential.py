"""
01 — Symbolic derivation of the SERIES Weber force, potential, Lagrangian,
     conserved energy, and velocity-dependent (Finsler) effective mass.

Series ansatz (Assis 1992, Can. J. Phys. 70, 330, Eqs. [6]-[7]):

    U(r, rdot) = (q/r) * g(rdot^2/c^2),   g(s) = 1 - a1 s - a2 s^2 - a3 s^3 - ...

with (Assis notation) a1 = alpha, a2 = beta, a3 = gamma, ...
Weber  : alpha=1/2, beta=gamma=...=0
Phipps : g(s) = sqrt(1-s)  ->  alpha=1/2, beta=1/8, gamma=1/16, ...

Everything is verified with SymPy.  No physics is assumed beyond
    F = -hat r dU/dr      (Assis's work/energy route, d(rdot^2)/dr = 2 rddot)
    mu rddot_vec = F hat r   (Newton for the relative coordinate)
"""
import sympy as sp

r, c, q, mu = sp.symbols('r c q mu', positive=True)
rd, rdd = sp.symbols('rdot rddot', real=True)          # rdot = dr/dt, rddot = d^2 r/dt^2
a1, a2, a3 = sp.symbols('alpha beta gamma', real=True)  # Assis potential coefficients
w = sp.symbols('w', real=True)                          # placeholder for rdot^2

print("="*72)
print("SERIES WEBER ELECTRODYNAMICS — symbolic derivations")
print("="*72)

# ---------------------------------------------------------------- potential
g_of = lambda s: 1 - a1*s - a2*s**2 - a3*s**3      # generating function g(s)
U = q/r * g_of(w/c**2)                              # potential as function of (r, w=rdot^2)

# ---- (1) Force law via F = -dU/dr  with  d(rdot^2)/dr = 2 rddot ------------
# total derivative of U along the trajectory wrt r:  dU/dr = U_r + U_w * (dw/dr),
# and dw/dr = d(rdot^2)/dr = 2 rddot.
dUdr = sp.diff(U, r) + sp.diff(U, w)*(2*rdd)
F = sp.simplify((-dUdr).subs(w, rd**2))
print("\n[1] Series force law  F(r,rdot,rddot):")
sp.pprint(sp.collect(sp.expand(F*r**2/q), rdd))
print("    (printed as r^2 F / q)")

# Assis's published closed form, Eq. [7]:
F_assis = q/r**2 * (1
    - a1*(rd**2 - 2*r*rdd)/c**2
    - a2*(rd**4 - 4*r*rd**2*rdd)/c**4
    - a3*(rd**6 - 6*r*rd**4*rdd)/c**6)
print("\n    F - F_Assis[7]  simplifies to:", sp.simplify(F - F_assis), " (must be 0)")

# General n-th term check:  term_n = a_n * (q/r^2) * (rdot^{2n} - 2 n r rdot^{2n-2} rddot)/c^{2n}
# (verified for concrete n=1..5; symbolic-exponent simplification is unreliable in SymPy)
an = sp.symbols('a_n', real=True)
rdp = sp.symbols('rdot', positive=True)        # positive so rdot**(2n-2) simplifies cleanly
allmatch = True
for nn in range(1, 6):
    term_pot = q/r * (-an) * (w/c**2)**nn
    tf = sp.simplify((-(sp.diff(term_pot, r) + sp.diff(term_pot, w)*2*rdd)).subs(w, rdp**2))
    # Assis [7]: higher terms carry explicit minus with positive coeff an=alpha,beta,gamma
    tc = -an*q/r**2 * (rdp**(2*nn) - 2*nn*r*rdp**(2*nn-2)*rdd)/c**(2*nn)
    allmatch &= (sp.simplify(tf - tc) == 0)
print("    general n-th force term = a_n q/r^2 (rdot^2n - 2n r rdot^(2n-2) rddot)/c^2n "
      "for n=1..5 :", allmatch)

# ---- (2) Effective (velocity-dependent) mass  m_eff = mu - dF/d(rddot) -----
# Newton (l=0):  mu rddot = F(r,rdot,rddot)  =>  (mu - dF/d rddot) rddot = F|_{rddot=0}
B = sp.simplify(sp.diff(F, rdd))            # coefficient of rddot in F
m_eff = sp.simplify(mu - B)
print("\n[2] Effective inertial mass  m_eff(r,rdot) = mu - dF/d(rddot):")
sp.pprint(m_eff)
gp = sp.diff(g_of(w/c**2), w).subs(w, rd**2)        # g'(s) * (1/c^2)
print("    check  m_eff == mu + (2 q / (r c^2)) g'(s):",
      sp.simplify(m_eff - (mu + 2*q/(r*c**2)*sp.diff(g_of(sp.symbols('s')), sp.symbols('s')).subs(sp.symbols('s'), rd**2/c**2))) == 0)
print("    Weber limit (alpha=1/2,beta=gamma=0):",
      sp.simplify(m_eff.subs({a1: sp.Rational(1,2), a2: 0, a3: 0})),
      " =>  mu(1 - q/(mu r c^2)) = mu(1 - rho/r)")

# ---- (3) Conserved energy  E = 1/2 mu rdot^2 + U(r,rdot)  (l=0) ------------
# d/dt(1/2 mu rdot^2) = mu rddot rdot = F rdot = -(dU/dr) rdot = -dU/dt.
E = sp.Rational(1,2)*mu*rd**2 + (q/r)*g_of(rd**2/c**2)
print("\n[3] Conserved energy (l=0):  E = 1/2 mu rdot^2 + (q/r) g(rdot^2/c^2)")
print("    Weber:", sp.simplify(E.subs({a1: sp.Rational(1,2), a2:0, a3:0})))

# ---- (4) Series Weber Lagrangian: the 'odd-denominator' resummation --------
# Claim: L = 1/2 mu rdot^2 - (q/r) h(s),   h(s) = 1 + sum_{n>=1} a_n/(2n-1) s^n,
# reproduces the same force AND the same energy E.  (a_n = alpha,beta,gamma.)
h_of = lambda s: 1 + (a1/1)*s + (a2/3)*s**2 + (a3/5)*s**3   # NOTE denominators 1,3,5
L = sp.Rational(1,2)*mu*rd**2 - (q/r)*h_of(rd**2/c**2)
# Euler-Lagrange (radial, l=0):  d/dt(dL/drdot) - dL/dr = 0  ->  m_eff rddot = ...
Lrd = sp.diff(L, rd)
# replace explicit r,rdot time-derivatives:  d/dt -> rdot d/dr + rddot d/drdot
EL = sp.diff(Lrd, r)*rd + sp.diff(Lrd, rd)*rdd - sp.diff(L, r)
EL = sp.simplify(EL)                       # = 0 is the EOM:  m_eff rddot - (force terms) = 0
# Newton form from the force route:  m_eff rddot - F|_{rddot=0} = 0
Newton = sp.simplify(m_eff*rdd - F.subs(rdd, 0))
print("\n[4] Lagrangian L = 1/2 mu rdot^2 - (q/r)(1 + alpha s + beta/3 s^2 + gamma/5 s^3),  s=rdot^2/c^2")
print("    Euler-Lagrange EOM minus Newton/force EOM:", sp.simplify(EL - Newton), " (must be 0)")
# Energy from Lagrangian:  E_L = rdot dL/drdot - L  must equal E above.
E_L = sp.simplify(rd*Lrd - L)
print("    Jacobi energy rdot*L_rdot - L equals 1/2 mu rdot^2 + (q/r) g :",
      sp.simplify(E_L - E) == 0)
print("    => the SAME g appears in the energy; denominators 1,3,5,7,... are the signature.")

# ---- (5) Phipps closed form sanity: h and g for g(s)=sqrt(1-s) -------------
s = sp.symbols('s', real=True)
g_ph = sp.sqrt(1 - s)
# Assis coefficients from -Taylor of sqrt(1-s):  g = 1 - 1/2 s - 1/8 s^2 - 1/16 s^3 - ...
ser = sp.series(g_ph, s, 0, 5).removeO()
print("\n[5] Phipps g(s)=sqrt(1-s) Taylor:", ser, " => alpha=1/2,beta=1/8,gamma=1/16")
# h for Phipps: solve h - 2 s h' = g  with h(s)=sqrt(1+ s)?  check closed form h=sqrt(1+s)? no.
# The ODE h - 2 s h' = g has solution h(s) = -sum a_n/(2n-1) s^n (+ C sqrt s); for Phipps:
h_ph_series = sum(sp.Rational(1, 2*k-1)*sp.binomial(sp.Rational(1,2), k)*(-1)**0*0 for k in range(0)) # placeholder
print("    Phipps Lagrangian density h(s) series (1 + 1/2 s + 1/8/3 s^2 + ... ):")
# build h from g_ph Taylor: a_n = -coeff, h_n = a_n/(2n-1) for n>=1, h_0=1
gp_ser = sp.series(g_ph, s, 0, 6).removeO()
poly = sp.Poly(gp_ser, s)
hs = 1
for k in range(1, 6):
    a_k = -poly.coeff_monomial(s**k)          # Assis coeff (alpha,beta,...) = -taylor coeff
    hs += sp.Rational(a_k.p, a_k.q)/(2*k-1)*s**k if hasattr(a_k,'p') else a_k/(2*k-1)*s**k
print("    h(s) =", sp.nsimplify(hs))

print("\nDone.  All boxed checks above must report 0 / True.")
