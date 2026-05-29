"""
03 — Head-on (l=0) collision speed: the series REPAIRS Weber's superluminal value.

For l=0 a particle reaches the collision locus r->0 only if the velocity factor
g(s) -> 0 there (so that (1/r) g stays finite against the +1/r Coulomb barrier).
The collision speed is therefore  s* = smallest positive root of g(s)=0,
    rdot_collision / c = sqrt(s*).

Bare Weber       g = 1 - s/2            -> s*=2     -> rdot = sqrt(2) c  (SUPERLUMINAL)
Phipps (full)    g = sqrt(1-s)          -> s*=1     -> rdot =        c  (light barrier)
Truncations of the sqrt-series interpolate; we show s* -> 1 as order grows.
"""
import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from math import comb

os.makedirs(os.path.join(os.path.dirname(__file__), "figs"), exist_ok=True)

def sqrt_taylor_coeffs(N):
    """coeffs c_k of sqrt(1-s) = sum c_k s^k, k=0..N  (c_0=1, rest negative)."""
    # binom(1/2, k) * (-1)^k
    from math import gamma
    coeffs = []
    for k in range(N+1):
        # binom(1/2,k) = prod_{i=0}^{k-1} (1/2 - i) / k!
        num = 1.0
        for i in range(k):
            num *= (0.5 - i)
        from math import factorial
        b = num/factorial(k)
        coeffs.append(b*((-1)**k))
    return np.array(coeffs)

print("order N  | series g_N(s) = partial sum of sqrt(1-s)")
print(" (N=1 is bare Weber)        smallest positive root s*    rdot*/c = sqrt(s*)")
print("-"*70)
orders = list(range(1, 21))
roots = []
for N in orders:
    c = sqrt_taylor_coeffs(N)            # g_N(s) = sum c_k s^k
    r = np.roots(c[::-1])                # numpy wants highest-degree first
    real_pos = sorted([x.real for x in r if abs(x.imag) < 1e-8 and x.real > 1e-9])
    s_star = real_pos[0] if real_pos else np.nan
    roots.append(s_star)
    tag = "  <-- Weber" if N == 1 else ("  <-- Phipps coeffs through gamma" if N == 3 else "")
    if N <= 8 or N % 4 == 0:
        print(f"  N={N:2d}     s* = {s_star:8.5f}     rdot*/c = {np.sqrt(s_star):8.5f}{tag}")

print("-"*70)
print(f"  Phipps (full resummation)   s* = 1.00000     rdot*/c = 1.00000  (exactly c)")

fig, ax = plt.subplots(figsize=(7.0, 4.6))
ax.plot(orders, np.sqrt(roots), "o-", color="C0", label=r"truncated series $\dot r_*/c=\sqrt{s_*}$")
ax.axhline(np.sqrt(2), color="k", ls="--", lw=1.2, label=r"bare Weber $\sqrt{2}\,c$ (superluminal)")
ax.axhline(1.0, color="C3", ls="-", lw=1.4, label=r"Phipps / light barrier $c$")
ax.set_xlabel("series order N (highest power of $v^2/c^2$ kept)")
ax.set_ylabel(r"head-on collision speed  $\dot r_* / c$")
ax.set_title("Series resummation drives the l=0 collision speed from $\\sqrt{2}\\,c$ down to $c$")
ax.set_xticks(orders[::2]); ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "figs", "03_collision_speed.png")
fig.savefig(out, dpi=130); print("wrote", out)
