# Transformed Weber Hamiltonians

## Scope

This document catalogs canonical transformations and coordinate changes applied to the two-body Weber Hamiltonian. Each transformation provides a different perspective on the dynamics: some simplify the Hamiltonian, others regularize singularities, and yet others reveal hidden geometric structure. The goal is not merely simplification but deeper insight into the physical content of Weber electrodynamics.

All results were derived symbolically using the exploration script `scripts/weber_transforms.py`.

## Starting Point: The Weber Hamiltonian

For the reduced two-body problem in polar coordinates $(r, \phi)$ with angular momentum $\ell = p_\phi$ conserved, the Weber Hamiltonian takes two equivalent forms.

### Standard form

$$H = \frac{p_r^2}{2\mu}\left(1 - \frac{\rho}{r}\right) + \frac{\ell^2}{2\mu r^2} + \frac{k}{r}$$

where $k = q_1 q_2$, $\mu$ is the reduced mass, and $\rho = k/(\mu c^2)$ is the critical radius.

### Metric form (Frauenfelder--Weber 2024)

Absorbing the velocity-dependent potential into the kinetic energy:

$$H = \frac{r\,p_r^2}{2\mu(r - \rho)} + \frac{\ell^2}{2\mu r^2} + \frac{k}{r}$$

This is the Hamiltonian for a particle on the Weber plane with Riemannian/Lorentzian metric $g_{rr} = (r-\rho)/r$, $g_{\phi\phi} = r^2$, and a pure Coulomb potential $k/r$. The metric changes signature at $r = \rho$: Riemannian for $r > \rho$, Lorentzian for $r < \rho$. See [CriticalRadiusAndLikeChargeAttraction.md](CriticalRadiusAndLikeChargeAttraction.md) for details.

### Perturbative decomposition

$$H = H_0 + \rho\,H_1$$

where $H_0 = p_r^2/(2\mu) + \ell^2/(2\mu r^2) + k/r$ is the Kepler Hamiltonian and $H_1 = -p_r^2/(2\mu r)$ is the Weber correction.

---

## Part I: Coordinate Transformations

### 1. Radial Flattening Coordinate

**Motivation.** The metric form has non-standard kinetic energy $g^{rr} p_r^2/(2\mu)$ with $g^{rr} = r/(r-\rho)$. Can we find a coordinate in which the kinetic energy becomes standard?

**Transformation.** Define $r^*$ by

$$\frac{dr^*}{dr} = \frac{1}{\sqrt{g^{rr}}} = \sqrt{\frac{r - \rho}{r}}$$

This is a canonical point transformation with $p_r = P \cdot \sqrt{(r-\rho)/r}$.

**Result.** The Hamiltonian becomes

$$\boxed{H = \frac{P^2}{2\mu} + \frac{\ell^2}{2\mu r^2} + \frac{k}{r}}$$

which is *exactly* the Kepler/Coulomb Hamiltonian. All Weber physics is encoded in the nonlinear map $r(r^*)$.

**The flattening map.** Integrating $f'(r) = \sqrt{(r-\rho)/r}$:

$$r^*(r) = \sqrt{r(r-\rho)} - \rho\,\operatorname{arcosh}\!\left(\frac{\sqrt{r}}{\sqrt{\rho}}\right), \qquad r > \rho$$

For $r < \rho$ (sub-critical regime):

$$r^*(r) = \frac{\sqrt{r}\,(r - \rho)}{\sqrt{\rho - r}} + \rho\arcsin\!\left(\frac{\sqrt{r}}{\sqrt{\rho}}\right) \cdot i$$

The map is real for $r > \rho$ and acquires an imaginary part for $r < \rho$, reflecting the Lorentzian signature change.

**Physical interpretation.** The Weber two-body problem *is* the Kepler problem in disguise, but in a nonlinearly distorted radial coordinate. The distortion compresses (or stretches) the radial direction near $\rho$, creating an apparent "metric" effect. The critical radius $r = \rho$ maps to $r^* = 0$, the boundary of the coordinate domain.

**Note.** This differs from the Schwarzschild tortoise coordinate ($dr^*/dr = r/(r-\rho)$, which flattens the null-cone structure). Our coordinate flattens the kinetic energy.

### 2. Isothermal (Conformal) Coordinates on the Weber Plane

**Motivation.** Find coordinates in which the full 2D Weber metric is conformally flat.

**Transformation.** For $r > \rho$, define $\sigma$ by

$$\frac{d\sigma}{dr} = \frac{\sqrt{r - \rho}}{r^{3/2}}$$

so that the Weber metric becomes $ds^2 = r^2(d\sigma^2 + d\phi^2)$.

**Result.**

$$\sigma(r) = 2\operatorname{arcosh}\!\left(\frac{\sqrt{r}}{\sqrt{\rho}}\right) - \frac{2\sqrt{r-\rho}}{\sqrt{r}}$$

The Hamiltonian in isothermal coordinates is

$$H = \frac{p_\sigma^2 + p_\phi^2}{2\mu\, r(\sigma)^2} + \frac{k}{r(\sigma)}$$

where $r(\sigma)$ is the inverse of $\sigma(r)$ above.

**Physical interpretation.** Isothermal coordinates make the Weber plane conformally equivalent to flat space. All anisotropy from the Weber term is absorbed into the conformal factor $\Omega = r(\sigma)$. The kinetic energy is isotropic in $(\sigma, \phi)$ — the "speed of light" is the same in all directions.

### 3. Kruskal--Szekeres Analogy

**Motivation.** The Weber metric $g_{rr} = 1 - \rho/r$ has the same form as the Schwarzschild metric. The coordinate singularity at $r = \rho$ can be removed by the Kruskal--Szekeres (KS) construction.

**Transformation.** Define new coordinates $(T, X)$ by

$$T^2 - X^2 = \left(\frac{r}{\rho} - 1\right)e^{r/\rho}$$

with the conformal factor $F = 4\rho^3 e^{-r/\rho}/r$.

**Result.**

$$H = \frac{r(p_T^2 - p_X^2)\,e^{r/\rho}}{8\mu\rho^3} + \frac{\ell^2}{2\mu r^2} + \frac{k}{r}$$

where $r = r(T,X)$ is defined implicitly.

**Key property.** The conformal factor $F(\rho) = 4\rho^2/e$ is finite and nonzero. The coordinate singularity at $r = \rho$ is completely removed. The kinetic term has Lorentzian signature $(p_T^2 - p_X^2)$, reflecting the signature change: inside $\rho$, the $T$-direction is spacelike and the $X$-direction is timelike (roles swap).

**Physical interpretation.** The KS construction provides a maximally extended Weber phase space that covers both the distant and molecular regimes in a single coordinate chart. The critical radius $r = \rho$ is a regular point of the extended dynamics, not a true singularity.

### 4. Eddington--Finkelstein-like Coordinates

**Motivation.** Mix time and space in extended phase space to remove the coordinate singularity at $r = \rho$ while preserving a simpler structure than KS.

**Transformation.** In extended phase space $(r, t, p_r, p_t)$ with $p_t = -E$, define ingoing EF coordinate $v = t + r^*(r)$.

**Result.**

$$H_{\text{EF}} = \frac{P^2 r}{2\mu(r-\rho)} - \frac{2PP_v r^2}{2\mu(r-\rho)^2} + \frac{P_v^2 r^3}{2\mu(r-\rho)^3} + P_v + \frac{\ell^2}{2\mu r^2} + \frac{k}{r}$$

The constraint $H_{\text{EF}} = 0$ defines the physical trajectories. The cross-term $PP_v$ mixes the spatial and temporal momenta, which is the mechanism by which the singularity at $r = \rho$ is resolved.

### 5. Levi-Civita Regularization in the Metric Frame

**Motivation.** Apply the standard collision regularization $r = u^2$ to the metric-form Hamiltonian.

**Transformation.** $r = u^2$, $p_r = P/(2u)$ (canonical). Sundman time scaling $K = u^2(H - E)$.

**Result.** The LC transform regularizes $r = 0$ (Coulomb singularity) but **not** the critical radius $r = \rho$. The factor $(u^2 - \rho)$ persists in the denominator at $u = \sqrt{\rho}$.

**Physical interpretation.** The Coulomb and critical-radius singularities are structurally different. The Coulomb singularity is a coordinate artifact removable by LC/KS. The critical-radius singularity is a genuine metric signature change — it requires the KS-type construction (Section 3) to resolve, not the collision-type construction.

---

## Part II: Generating Function Experiments

### 6. Type-1: $F_1(r, R)$ --- Phase Space Exchange

**$F_1 = rR$ (identity/exchange).** The relations $p_r = R$, $P = -r$ exchange the roles of position and momentum:

$$H(R, P) = \frac{R^2}{2\mu}\left(1 + \frac{\rho}{P}\right) - \frac{k}{P} + \frac{\ell^2}{2\mu P^2}$$

This is the *momentum-space Weber Hamiltonian*. The Weber correction $\rho/P$ is a $1/P$ perturbation of the kinetic energy, whereas in position space it was a $1/r$ perturbation. For large momenta (high energy), the Weber correction vanishes — the system approaches Kepler.

**$F_1 = R \cdot r^*(r)$ (tortoise as generating function).** This implements the Schwarzschild-type tortoise coordinate via an $F_1$ generating function. The new momentum $R = p_r(1 - \rho/r)$ is the "tortoise momentum" and $P = -r^*$ is minus the tortoise coordinate.

### 7. Type-2: $F_2(r, P)$ --- Hamilton-Jacobi Framework

**Action generator (formal).** From $H = E$:

$$p_r^2 = \frac{2\mu E\,r^2 - 2\mu k\,r - \ell^2}{r(r - \rho)}$$

The integrand $\sqrt{p_r^2}$ has:
- Branch points at the roots of the numerator (turning points of the orbit)
- Poles at $r = 0$ and $r = \rho$ (from the Weber metric)

The action variable $J = (1/2\pi)\oint p_r\,dr$ is an elliptic-type integral with the additional pole at $r = \rho$ that is absent from Kepler.

**Perturbative expansion.** The Weber $p_r^2$ can be expanded:

$$p_{r,\text{Weber}}^2 \approx p_{r,\text{Kepler}}^2 \cdot \left(1 + \frac{\rho}{r} + \frac{\rho^2}{r^2} + \cdots\right)$$

giving a systematic series of corrections to the Kepler action.

**Sigmoid regularization of the sub-critical domain.** The $F_2$ generating function

$$F_2 = P \cdot \ln\!\left(\frac{r}{\rho - r}\right)$$

maps $r \in (0, \rho)$ to $R \in (-\infty, +\infty)$ via $r = \rho\,\sigma(R)$ where $\sigma$ is the sigmoid function. The momentum transforms as $p_r = P\rho/(r(\rho-r))$. This compactifies the finite sub-critical interval to the full real line — a regularization of the domain rather than the singularity.

### 8. Type-3: $F_3(p_r, R)$ --- Momentum-Space Transformations

The fundamental obstacle for $F_3$-type transformations: the Weber term $p_r^2 \rho/(2\mu r)$ couples momentum and position multiplicatively. No $F_3(p_r, R)$ can fully decouple this into a standard kinetic energy $P^2/(2\mu)$ plus a position-only potential, because the $F_3$ relations $r = -\partial F_3/\partial p_r$ and $P = -\partial F_3/\partial R$ cannot simultaneously eliminate the $p_r^2/r$ coupling.

An $F_3$ of the form $F_3 = -p_r R - \rho p_r \ln(R/\rho)$ implements the inverse tortoise map $r = R + \rho\ln(R/\rho)$, but the resulting expressions are transcendental and no simpler than direct point transformation methods.

### 9. Type-4: $F_4(p_r, P)$ --- Momentum-to-Momentum

**Phase-space exchange** ($F_4 = p_r P$) gives the same result as $F_1 = rR$: the momentum-space Hamiltonian.

**Lorentz boost.** A linear symplectic transformation with boost parameter $\xi$:

$$r_{\text{new}} = r\cosh\xi - \frac{p_r}{\mu}\sinh\xi, \qquad p_{\text{new}} = -\mu r\sinh\xi + p_r\cosh\xi$$

mixes position and momentum. This is the phase-space analogue of a Lorentz boost, motivated by the Lorentzian signature inside $r < \rho$. By choosing $\xi$ appropriately, one can eliminate specific terms in the Hamiltonian at leading order in $\rho$.

### 10. Automated $F_2$ Scan

A systematic scan of $F_2 = P \cdot r^a(r-\rho)^b$ for $a, b \in \{-1, -\frac{1}{2}, 0, \frac{1}{2}, 1, 2\}$ found 10 invertible cases. The simplest results:

| $a$ | $b$ | New $r$ | Transformed $H$ |
|-----|-----|---------|-----------------|
| 2 | 0 | $r = p_r/(2P)$ | $\frac{p_r^2}{2\mu} + \frac{2P\ell^2}{p_r^2\mu} + \frac{2Pk}{p_r} - \frac{P\rho}{mu}$ |
| 0 | 2 | $r = \rho + p_r/(2P)$ | Centers sub-critical domain at origin |
| 1 | 1 | $r = (P\rho + p_r)/(2P)$ | Shifts by $\rho/2$ |

The shift $F_2 = P(r - \rho/2)$ is particularly natural for the sub-critical regime, placing the midpoint of $[0, \rho]$ at the origin.

---

## Part III: Perturbation Theory and Action-Angle Variables

### 11. Deprit--Lie Perturbation Theory

**Decomposition.** Writing $H = H_0 + \rho H_1$ with $H_0$ the Kepler Hamiltonian, the Weber correction on the Kepler orbit becomes:

$$H_1 = -\frac{p_r^2}{2\mu r} = -\frac{E}{r} + \frac{\ell^2}{2\mu r^3} + \frac{k}{r^2}$$

using $p_r^2 = 2\mu(E - k/r) - \ell^2/r^2$ on the Kepler shell.

**Orbit average.** Using the standard Kepler averages $\langle 1/r \rangle = 1/a$, $\langle 1/r^2 \rangle = 1/(a^2(1-e^2))$, $\langle 1/r^3 \rangle = 1/(a^3(1-e^2)^{3/2})$:

$$\langle H_1 \rangle = -\frac{E}{a} + \frac{\ell^2}{2\mu a^3(1-e^2)^{3/2}} + \frac{k}{a^2(1-e^2)}$$

Substituting the Kepler relations $E = -k/(2a)$ and $\ell^2 = \mu k a(1-e^2)$:

$$\langle H_1 \rangle = \frac{k}{2a^2(1-e^2)^{3/2}}\left(1 - e^2 + 2\sqrt{1-e^2} + (1-e^2)^{3/2}\right)$$

This is the *secular Weber Hamiltonian* — the orbit-averaged perturbation that governs slow evolution of orbital elements.

### 12. Delaunay Variables

In Delaunay actions $(L, G)$ with $L = \sqrt{-\mu k a}$, $G = \ell = L\sqrt{1-e^2}$:

$$H_0(L) = -\frac{\mu k^2}{2L^2}$$

$$\langle H_1 \rangle(L, G) = \frac{k^3\mu^2(2L^2 - G^2 - GL)}{2G^2 L^4}$$

**Apsidal precession.** The precession rate of the argument of perihelion is:

$$\frac{dg}{dt} = \rho \cdot \frac{\partial \langle H_1 \rangle}{\partial G} = \rho \cdot \frac{k^3\mu^2(G - 4L)}{2G^3 L^3}$$

The precession per orbit is:

$$\Delta g = 2\pi\rho \cdot \frac{k\mu(G - 4L)}{2G^3}$$

For nearly circular orbits ($G \approx L$, $e \approx 0$): $\Delta g \approx -3\pi\rho k\mu/L^3$, a retrograde precession.

### 13. McGehee Blow-Up Coordinates

**$r = s^2$ (Levi-Civita-type).** The Hamiltonian becomes:

$$H = \frac{P^2}{8\mu s^2}\left(1 - \frac{\rho}{s^2}\right) + \frac{\ell^2}{2\mu s^4} + \frac{k}{s^2}$$

The Weber factor $(1 - \rho/s^2)$ diverges at $s = 0$, which is qualitatively different from Kepler. This reflects the finite-speed collision ($|\dot{r}| \to \sqrt{2}c$ at $r = 0$) — the kinetic energy blows up as $1/s^2$ rather than being regular.

**$r = 1/s^2$ (blow-up).** Transforms $r = 0$ to $s = \infty$:

$$H = \frac{P^2 s^6}{8\mu}\left(1 - \rho s^2\right) + \frac{\ell^2 s^4}{2\mu} + ks^2$$

The leading term at large $s$ is $-P^2\rho s^8/(8\mu)$ — a negative-definite term that dominates near collision. In Kepler, the leading term is $+P^2 s^6/(8\mu)$ (positive). This sign difference is the McGehee-coordinate manifestation of the negative effective mass in the sub-critical regime.

### 14. Hamilton--Jacobi Separation

The Hamilton--Jacobi equation for Weber is:

$$(W')^2 = \frac{2\mu E\,r^2 - 2\mu k\,r - \ell^2}{r(r - \rho)}$$

**Separability diagnosis.** In polar coordinates, $\phi$ separates trivially (yielding $\ell$). The radial part is a single ODE with an extra pole at $r = \rho$ compared to Kepler. True non-trivial separation (as in parabolic coordinates for Kepler + Stark) would require a hidden constant of motion analogous to the Runge--Lenz vector. For Weber, the Runge--Lenz vector is *not* conserved — it precesses at the rate computed in Section 12.

**Parabolic coordinates.** The Weber Hamiltonian in Cartesian 2D involves the term $(x p_x + y p_y)^2/r^3$, which couples $x$ and $y$ through the radial velocity squared. This is quadratic in momenta (unlike the Stark effect which is linear in coordinates), and is not separable in parabolic coordinates.

### 15. Action-Angle Variables for Sub-Critical Oscillation

For the molecular state ($r_0 < \rho$, $\ell = 0$, $E = k/r_0$):

$$p_r^2 = \frac{2\mu k(r_0 - r)}{r_0(\rho - r)}$$

The action integral is

$$J = \frac{1}{\pi}\int_0^{r_0} \sqrt{\frac{2\mu k(r_0 - r)}{r_0(\rho - r)}}\,dr$$

This is an elliptic integral with a branch point at $r = r_0$ and a pole at $r = \rho > r_0$.

**Small-amplitude limit** ($r_0 \ll \rho$): In this regime $\rho - r \approx \rho$, so the integrand simplifies:

$$J \approx \frac{1}{\pi}\sqrt{\frac{2\mu k}{\rho}}\int_0^{r_0} \sqrt{r_0 - r}\,dr = \frac{2}{3\pi}\sqrt{\frac{2\mu k}{\rho}}\,r_0^{3/2}$$

Inverting: $r_0 \propto J^{2/3}$, hence $E = k/r_0 \propto J^{-2/3}$.

**Frequency.** $\omega = dE/dJ \propto J^{-5/3}$, which is anharmonic (frequency depends on amplitude). This is consistent with Weber's period formula $T = 2\sqrt{2}\,r_0/c$ which gives $\omega = \pi c/(\sqrt{2}\,r_0)$.

The power-law $E \propto J^{-2/3}$ should be compared with:
- Kepler: $E \propto J^{-2}$ (attractive case)
- Harmonic oscillator: $E \propto J$ (linear)
- Particle in a box: $E \propto J^2$ (quadratic)

The Weber sub-critical oscillation has a weaker dependence on $J$ than Kepler, reflecting the fact that the effective potential well becomes deeper (not shallower) as $r \to 0$.

### 16. Jacobi--Maupertuis Metric and Gaussian Curvature

On the energy surface $H = E$, trajectories are geodesics of the Jacobi--Maupertuis metric:

$$ds^2_{\text{JM}} = 2\mu(E - k/r)\left[\frac{r - \rho}{r}\,dr^2 + r^2\,d\phi^2\right]$$

The metric components are:

$$g_{rr}^{\text{JM}} = \frac{2\mu(r-\rho)(Er - k)}{r^2}, \qquad g_{\phi\phi}^{\text{JM}} = 2\mu r(Er - k)$$

**Kepler limit.** For $\rho = 0$, the Gaussian curvature simplifies to:

$$K_{\text{Kepler}} = \frac{Ek}{4\mu(Er - k)^3}$$

For $E < 0$ (bound) and $k < 0$ (attractive): $K > 0$ (focusing) outside the turning point, $K < 0$ (defocusing) inside.

**Critical radius.** At $r = \rho$: $g_{rr}^{\text{JM}} = 0$, so the metric is degenerate and the Gaussian curvature diverges. This is the geometric manifestation of the critical-radius barrier — the JM metric has a metric singularity at $r = \rho$, directly analogous to the event horizon in Schwarzschild spacetime.

**Sub-critical regime.** For $r < \rho$: $g_{rr}^{\text{JM}} < 0$, the JM metric becomes Lorentzian, and the standard notion of Gaussian curvature does not apply. The geometry is pseudo-Riemannian, and the "geodesics" have the character of timelike curves in the Lorentzian sense.

---

## Summary: Map of Transformations

```
                    Flattening coordinate
    Weber (metric) ────────────────────────> Kepler (exactly!)
         │                                       │
         │ LC (r=u²)                              │ Delaunay
         ▼                                        ▼
    LC-regularized                          Action-angle
    (removes r=0,                           (orbit averages,
     keeps r=ρ)                              precession rate)
         │
         │ KS analogy
         ▼
    Kruskal-Szekeres
    (removes r=ρ,
     Lorentzian kinetic)

    Weber (standard) ──── F1=rR ────> Momentum-space Weber
         │
         │ F2 sigmoid
         ▼
    Sub-critical on ℝ
    (r ∈ (0,ρ) → R ∈ ℝ)
```

### Key Results

1. **The Weber problem is exactly Kepler** in the flattening coordinate $r^* = \sqrt{r(r-\rho)} - \rho\,\operatorname{arcosh}(\sqrt{r/\rho})$. All Weber dynamics arises from the nonlinear map $r \leftrightarrow r^*$.

2. **The critical radius is a coordinate singularity**, removable by KS-type coordinates. The metric signature change from Riemannian to Lorentzian is the invariant content.

3. **Apsidal precession** in Delaunay variables: $\Delta g = 2\pi\rho \cdot k\mu(G-4L)/(2G^3)$ per orbit, retrograde for nearly circular orbits.

4. **Sub-critical action-angle**: $E \propto J^{-2/3}$, anharmonic with $\omega \propto J^{-5/3}$. The action integral is elliptic, with a pole at $r = \rho$ absent from Kepler.

5. **The JM metric** has a metric singularity at $r = \rho$ (divergent Gaussian curvature) and becomes Lorentzian for $r < \rho$. This is the full geometric encoding of Weber dynamics.

6. **No $F_3$ or $F_4$ generating function** can fully decouple the Weber kinetic term $p_r^2/r$ — the multiplicative coupling of momentum and position is intrinsic. Only point transformations (which are special cases of $F_1$ or $F_2$) can achieve full decoupling, and only by introducing implicit/transcendental coordinate maps.

---

## References

- Frauenfelder, U., Weber, J. "A mathematical description of the Weber nucleus as a classical and quantum mechanical system." *Anal. Math. Phys.* **14**:31 (2024).
- Goldstein, H., Poole, C., Safko, J. *Classical Mechanics*. 3rd ed. Addison-Wesley, 2002. Chapters 9--10 (canonical transformations, Hamilton-Jacobi theory).
- McGehee, R. "Triple collision in the collinear three-body problem." *Invent. Math.* **27** (1974), 191--227.
- Deprit, A. "Canonical transformations depending on a small parameter." *Celest. Mech.* **1** (1969), 12--30.

See also: [WeberElectrodynamics.md](../theory/WeberElectrodynamics.md), [CriticalRadiusAndLikeChargeAttraction.md](CriticalRadiusAndLikeChargeAttraction.md), [Regularization.md](../theory/Regularization.md).
