"""
04 — Topology of the energy level sets inside the nucleus (l != 0).

Conserved energy (mu=q=c=1, so r_c=1):
    E(r, rdot) = 1/2 rdot^2 + l^2/(2 r^2) + (1/r) g(rdot^2).

We contour E at a fixed level h for two models and read off the TOPOLOGY of the
orbit's phase curve {E=h}:

  * bare Weber  g = 1 - s/2 :  g -> -inf as s->inf, so the (1/r)g term can cancel
    the centrifugal +l^2/2r^2 barrier  ->  level set is OPEN, dives to r=0 (collision).
  * Phipps      g = sqrt(1-s), |g|<=1 :  E >= l^2/(2 r^2)  ->  hard floor
    r >= l/sqrt(2h).  The level set is a CLOSED loop  ->  a genuine bound (periodic)
    orbit.  The terminal spiral of bare Weber becomes a compact invariant curve.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

here = os.path.dirname(__file__)
os.makedirs(os.path.join(here, "figs"), exist_ok=True)

ell = 0.30
# common IC used in 05: at rest (rdot=0) at r0=0.6  ->  h = l^2/(2 r0^2) + 1/r0
r0 = 0.60
h = ell**2/(2*r0**2) + 1.0/r0
print(f"l = {ell},  level h = {h:.4f}  (rest at r0={r0})")

def E_weber(r, rd):  return 0.5*rd**2 + ell**2/(2*r**2) + (1.0/r)*(1 - 0.5*rd**2)
def E_phipps(r, rd):
    s = rd**2
    g = np.sqrt(np.clip(1 - s, 0, None))
    return 0.5*rd**2 + ell**2/(2*r**2) + (1.0/r)*g

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))

# ---- bare Weber: rdot ranges widely (collision at sqrt2 c for l=0, inf for l!=0)
r = np.linspace(0.02, 1.05, 600)
rd = np.linspace(-2.2, 2.2, 600)
R, RD = np.meshgrid(r, rd)
axes[0].contour(R, RD, E_weber(R, RD), levels=[h], colors="k", linewidths=2)
axes[0].contourf(R, RD, E_weber(R, RD), levels=[h, 1e9], colors=["#ffe9e9"])
axes[0].axvline(1.0, color="gray", ls=":", label="$r_c=1$")
axes[0].set_title("bare Weber  $g=1-s/2$\nlevel set is OPEN — dives to collision $r\\to0$")
axes[0].set_xlabel("$r$"); axes[0].set_ylabel("$\\dot r/c$")
axes[0].set_xlim(0, 1.05); axes[0].set_ylim(-2.2, 2.2); axes[0].legend(loc="upper right", fontsize=8)

# ---- Phipps: |rdot|<c, closed loop
r2 = np.linspace(0.02, 1.05, 600)
rd2 = np.linspace(-0.999, 0.999, 600)
R2, RD2 = np.meshgrid(r2, rd2)
axes[1].contour(R2, RD2, E_phipps(R2, RD2), levels=[h], colors="C3", linewidths=2)
rmin = ell/np.sqrt(2*(h-0.5))      # inner edge: loop touches rdot=c here
rfloor = ell/np.sqrt(2*h)          # looser absolute bound from E >= l^2/2r^2
axes[1].axvline(rmin, color="C0", ls="--",
                label=f"light barrier $r_{{min}}=\\ell/\\sqrt{{2(h-1/2)}}={rmin:.3f}$ ($\\dot r=c$)")
axes[1].axvline(1.0, color="gray", ls=":")
axes[1].set_title("Phipps  $g=\\sqrt{1-s}$\nlevel set is a CLOSED loop — bound periodic orbit")
axes[1].set_xlabel("$r$"); axes[1].set_ylabel("$\\dot r/c$")
axes[1].set_xlim(0, 1.05); axes[1].set_ylim(-1.0, 1.0); axes[1].legend(loc="upper right", fontsize=8)

fig.suptitle("Resummation changes the topology of the interior energy level set "
             "(l != 0): non-compact  ->  compact", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = os.path.join(here, "figs", "04_phase_topology.png")
fig.savefig(out, dpi=130); print("wrote", out)

print(f"\nPhipps inner light barrier (rdot=c):  r_min = l/sqrt(2(h-1/2)) = {rmin:.4f}")
print(f"Phipps absolute floor (E>=l^2/2r^2):  r >= l/sqrt(2h)        = {rfloor:.4f}  (no collision)")
print("bare Weber: no floor — (1/r)g -> -inf is unbounded below, collision survives.")
