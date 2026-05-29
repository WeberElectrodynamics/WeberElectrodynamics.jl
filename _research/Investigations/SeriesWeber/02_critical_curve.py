"""
02 — The critical RADIUS becomes a critical CURVE.

Frauenfelder & Weber: bare Weber has a velocity-INDEPENDENT critical radius
r_c = 1/c^2 (here, with mu=q=c=1, r_c = 1).  The metric g_rr = (r-r_c)/r changes
signature on the vertical line r = r_c in the (r, rdot) phase plane.

For the SERIES force the fibre-Hessian (effective mass) is velocity dependent,
    m_eff(r, rdot) = mu + (2 q / (r c^2)) g'(rdot^2/c^2),
so the signature-change locus m_eff = 0 is the CURVE
    r_crit(rdot) = -(2 q/(mu c^2)) g'(rdot^2/c^2)        (Finsler signature change).
With mu=q=c=1:  r_crit(rdot) = -2 g'(rdot^2).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

os.makedirs(os.path.join(os.path.dirname(__file__), "figs"), exist_ok=True)

# models: g'(s) with s = rdot^2  (mu=q=c=1, so r_c=1)
models = {
    "Weber (alpha=1/2)":            lambda s: -0.5 + 0*s,
    "+ beta s^2 (4th order)":       lambda s: -0.5 - 0.25*s,
    "+ gamma s^3 (6th order)":      lambda s: -0.5 - 0.25*s - (3/16)*s**2,
    "Phipps  g=sqrt(1-s)":          lambda s: -0.5/np.sqrt(np.clip(1-s, 1e-9, None)),
}
rcrit = {name: (lambda s, gp=gp: -2*gp(s)) for name, gp in models.items()}

rd = np.linspace(0.0, 0.999, 400)          # radial speed (units of c)
fig, ax = plt.subplots(figsize=(7.2, 5.2))
colors = {"Weber (alpha=1/2)": "k", "+ beta s^2 (4th order)": "C0",
          "+ gamma s^3 (6th order)": "C1", "Phipps  g=sqrt(1-s)": "C3"}
for name, rc in rcrit.items():
    s = rd**2
    rr = rc(s)
    ax.plot(rr, rd, color=colors[name], lw=2, label=name)
ax.axvline(1.0, color="gray", ls=":", lw=1)
ax.set_xlabel(r"critical radius $r_{\rm crit}$   (units $r_c=1/c^2$)")
ax.set_ylabel(r"radial speed  $\dot r / c$")
ax.set_title("Signature-change locus: Weber line vs series curves\n"
             "(left of each curve = Lorentzian / negative-mass / binding region)")
ax.set_xlim(0.7, 5.0); ax.set_ylim(0, 1.0)
ax.legend(loc="upper right", fontsize=9)
ax.text(1.02, 0.05, "Weber: vertical line\n(velocity-independent)", fontsize=8, color="k")
fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "figs", "02_critical_curve.png")
fig.savefig(out, dpi=130)
print("wrote", out)

# numeric anchors
for name, rc in rcrit.items():
    print(f"  {name:28s}  r_crit(0)={rc(0.0):.4f}   r_crit(0.5c)={rc(0.25):.4f}   "
          f"r_crit(0.9c)={rc(0.81):.4f}")
print("\nPhipps:  r_crit(rdot) = 1/sqrt(1 - rdot^2/c^2)  ->  diverges as rdot->c "
      "(relativistic gamma-dilation of the Weber radius).")
