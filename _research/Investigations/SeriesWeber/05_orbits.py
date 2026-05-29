"""
05 — Direct orbit integration inside the nucleus, l != 0.

Planar relative EOM (mu=q=c=1, critical radius r_c=1):
    rddot = [ l^2/r^3 + (1/r^2) g(s) ] / m_eff(r,rdot),   s = rdot^2,
    m_eff(r,rdot) = 1 + (2/r) g'(s),       phidot = l / r^2.

Same sub-critical IC (rest at r0=0.6, l=0.3) for both models.

  * bare Weber  g = 1 - s/2 :  the velocity factor is unbounded below, so the orbit
    SPIRALS INTO THE ORIGIN (terminal collision, infinite speed/winding) — exactly
    Frauenfelder & Weber's result: no periodic interior orbits.

  * Phipps      g = sqrt(1-s) :  |g|<=1.  The orbit can NOT reach r=0.  Instead it is
    confined to an annulus [r_min, r_max] and reaches an INNER LIGHT BARRIER at
        r_min = l / sqrt(2(h - 1/2))   where  rdot = c,
    a finite-speed (=c) reflection, C^0-continuable — the direct analogue of
    Frauenfelder's l=0 sqrt(2)c bounce, but now at finite separation and at l != 0.
    Reflecting rdot -> -rdot at the s=1 surface (energy-conserving) traces a bounded
    precessing rosette: a genuine collision-free bound state of the Weber nucleus.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

here = os.path.dirname(__file__)
os.makedirs(os.path.join(here, "figs"), exist_ok=True)

ell = 0.30
r0, vr0 = 0.60, 0.0
h = ell**2/(2*r0**2) + 1.0/r0
print(f"l={ell}, h={h:.4f}")

# ---- bare Weber: integrate to collision ----------------------------------
def rhs_weber(t, y):
    r, vr, phi = y
    s = vr*vr
    meff = 1.0 - 1.0/r                      # g'=-1/2
    num = ell**2/r**3 + (1 - 0.5*s)/r**2
    return [vr, num/meff, ell/r**2]
def hit_origin(t, y): return y[0] - 1e-4
hit_origin.terminal = True; hit_origin.direction = -1
solW = solve_ivp(rhs_weber, (0, 4.0), [r0, vr0, 0.0], events=hit_origin,
                 rtol=1e-11, atol=1e-13, max_step=2e-4)
rW, vrW, phiW = solW.y
EW = 0.5*vrW**2 + ell**2/(2*rW**2) + (1/rW)*(1-0.5*vrW**2)
print(f"Weber : r_min={rW.min():.2e} winding={(phiW[-1]-phiW[0])/(2*np.pi):.1f}rev "
      f"max|vr|={np.abs(vrW).max():.1f} reached_origin={len(solW.t_events[0])>0} "
      f"Edrift={np.ptp(EW):.1e}")

# ---- Phipps: integrate with light-barrier bounce at s = 1 ----------------
SCAP = 1.0 - 1e-7
def rhs_phipps(t, y):
    r, vr, phi = y
    s = min(vr*vr, SCAP)
    root = np.sqrt(1.0 - s)
    meff = 1.0 - 1.0/(r*root)               # g'=-1/(2 root)
    num = ell**2/r**3 + root/r**2
    return [vr, num/meff, ell/r**2]
def hit_barrier(t, y): return (y[1]*y[1]) - SCAP
hit_barrier.terminal = True; hit_barrier.direction = +1
def hit_rmax(t, y): return y[1]            # vr=0 (outer turning), not terminal
hit_rmax.terminal = False

state = [r0, vr0, 0.0]
T = []; R = []; VR = []; PHI = []; t_now = 0.0; nb = 0
for seg in range(400):                      # many radial half-periods
    sol = solve_ivp(rhs_phipps, (t_now, t_now+20), state, events=hit_barrier,
                    rtol=1e-10, atol=1e-12, max_step=5e-3, dense_output=False)
    T.append(sol.t); R.append(sol.y[0]); VR.append(sol.y[1]); PHI.append(sol.y[2])
    t_now = sol.t[-1]
    state = [sol.y[0,-1], -sol.y[1,-1], sol.y[2,-1]]   # reflect rdot at light barrier
    if len(sol.t_events[0]) == 0: break                # finished window w/o bounce
    nb += 1
    if t_now > 40: break
rP = np.concatenate(R); vrP = np.concatenate(VR); phiP = np.concatenate(PHI)
EP = 0.5*vrP**2 + ell**2/(2*rP**2) + (1/rP)*np.sqrt(np.clip(1-vrP**2,0,None))
rmin_pred = ell/np.sqrt(2*(h-0.5))
print(f"Phipps: r in [{rP.min():.4f},{rP.max():.4f}]  r_min_pred={rmin_pred:.4f}  "
      f"max|vr|={np.abs(vrP).max():.4f}  bounces={nb}  "
      f"winding={(phiP[-1]-phiP[0])/(2*np.pi):.1f}rev  Edrift={np.ptp(EP):.1e}")

# ---- figure --------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
xW, yW = rW*np.cos(phiW), rW*np.sin(phiW)
axes[0].plot(xW, yW, "k", lw=0.5)
axes[0].plot(xW[0], yW[0], "o", color="C2", ms=6, label="start")
axes[0].plot(0, 0, "r+", ms=10, label="collision r=0")
th = np.linspace(0, 2*np.pi, 200)
axes[0].plot(np.cos(th), np.sin(th), color="gray", ls=":", lw=1, label="$r_c=1$")
axes[0].set_aspect("equal"); axes[0].legend(fontsize=8)
axes[0].set_title("bare Weber: terminal spiral into the origin")

xP, yP = rP*np.cos(phiP), rP*np.sin(phiP)
axes[1].plot(xP, yP, "C3", lw=0.4)
axes[1].plot(xP[0], yP[0], "o", color="C2", ms=6, label="start")
for rr, cc, lab in [(rmin_pred, "C0", f"inner light barrier $r_{{min}}={rmin_pred:.3f}$"),
                    (rP.max(), "C1", f"outer turning $r_{{max}}={rP.max():.3f}$"),
                    (1.0, "gray", "$r_c=1$")]:
    axes[1].plot(rr*np.cos(th), rr*np.sin(th), color=cc, ls=":", lw=1.2, label=lab)
axes[1].plot(0, 0, "k+", ms=8)
axes[1].set_aspect("equal"); axes[1].legend(fontsize=8)
axes[1].set_title("Phipps: bounded rosette in an annulus (no collision)")
for ax in axes: ax.set_xlabel("x"); ax.set_ylabel("y")
fig.suptitle(f"Same sub-critical IC (l={ell}, rest at r={r0}): collapse vs bound state",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = os.path.join(here, "figs", "05_orbits.png")
fig.savefig(out, dpi=130); print("wrote", out)
