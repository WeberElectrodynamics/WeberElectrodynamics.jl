# Non-Regularizability of the Angular Momentum Collision in Weber Electrodynamics

## Scope

This document investigates whether the spiraling collision singularity arising when two like charges with angular momentum $\ell \neq 0$ collide inside the critical radius can be regularized by any smooth coordinate-time transformation. Seven regularization approaches are tested analytically and numerically. All fail: the obstruction is topological (infinite winding number) and fundamental. The document records what was tried, why each approach fails, and the lessons learned.

For the basic theory of the critical radius and the two collision types, see [CriticalRadiusAndLikeChargeAttraction.md](CriticalRadiusAndLikeChargeAttraction.md). For the standard regularization techniques used here, see [Regularization.md](../theory/Regularization.md).

## The Problem

### Setup

Two like charges inside the critical radius $\rho$ on the Weber plane, in dimensionless units ($k = ee'/\mu = 1$, $c = 1$, $\rho = 1$). The Hamiltonian is

$$H = \frac{r\,p_r^2}{2(r - 1)} + \frac{\ell^2}{2r^2} + \frac{1}{r},$$

and the energy equation (Frauenfelder & Weber 2024, Theorem 2.1) is

$$\dot{r}^2 = \frac{\ell^2 + 2r - 2hr^2}{r(1 - r)}.$$

### Two collision types

| | $\ell = 0$ (head-on) | $\ell \neq 0$ (spiraling) |
|---|---|---|
| Speed at $r = 0$ | $\|\dot{r}\| \to \sqrt{2}\,c$ (finite) | $\|\dot{r}\| \sim \ell/\sqrt{r} \to \infty$ |
| Angular velocity | $\dot\phi = 0$ | $\dot\phi = \ell/r^2 \to \infty$ |
| Total winding | 0 | $\phi \to \infty$ (infinite revolutions) |
| Collision time | finite | finite |
| Regularizable? | **yes** (Sundman/LC) | **no** (this document) |

### Power-law asymptotics ($\ell \neq 0$)

As $r \to 0$ with $\ell \neq 0$:

$$\dot{r} \sim \frac{\ell}{\sqrt{\rho\,r}} \sim r^{-1/2}, \qquad \dot\phi = \frac{\ell}{r^2} \sim r^{-2}, \qquad r \sim (T_{\text{coll}} - t)^{2/3}.$$

The $r^{-1/2}$ radial divergence is the same as in the Kepler problem. The critical difference is that in Kepler, angular momentum prevents the particle from reaching $r = 0$ (centrifugal barrier). In Weber inside $\rho$, the Lorentzian metric *inverts* the centrifugal barrier: angular momentum drives the particle *inward*, and the collision occurs in finite time.

Numerical verification (Experiment 0): fitting $|\dot{r}|$ vs $r$ in the range $10^{-6} < r < 0.05$ yields exponents $-0.36$ to $-0.48$ (approaching the theoretical $-0.5$ as $\ell$ increases); fitting $r$ vs $(T_{\text{coll}} - t)$ yields exponents $0.67$--$0.70$ (approaching $2/3$).

## The Topological Obstruction

The total angle accumulated as the particle spirals to the origin is

$$\phi = \int_0^{T_{\text{coll}}} \frac{\ell}{r(t)^2}\,dt.$$

Substituting $dt = -\sqrt{\rho\,r}\,dr/\ell$ (from the leading-order energy equation):

$$\phi \sim \int_0^{r_0} \frac{\ell}{r^2} \cdot \frac{\sqrt{\rho\,r}}{\ell}\,dr = \sqrt{\rho}\int_0^{r_0} r^{-3/2}\,dr = \infty.$$

The integral diverges: the particle makes *infinitely many revolutions* before reaching $r = 0$.

This is a **topological invariant**. The winding number of a curve around the origin is preserved by all smooth (and in particular all conformal) diffeomorphisms of the punctured plane. No smooth coordinate change can reduce an infinite winding number to a finite one.

**Consequence.** In any coordinate system where the origin corresponds to a single point, the trajectory must wind infinitely many times around that point. For the solution to be extendable through the collision as a smooth curve, it would need to arrive at the collision point with a well-defined tangent direction — but infinite winding means the tangent direction is undefined.

## Survey of Regularization Approaches

### Experiment 0: Baseline

**Method.** Direct integration of the polar ODE with `scipy.integrate.solve_ivp` (Radau), event detection at $r_{\min} = 10^{-8}$.

**Results.** For $r_0 = 0.5$, $\ell = 0.1$: collision at $t = 0.561$, $\sim 3179$ revolutions, final $|\dot{r}| = 10^3$. For $\ell = 0$: collision at $t = 0.574$, final $|\dot{r}| = 1.414 \approx \sqrt{2}$ (finite). Confirms all analytical predictions.


### Experiment 1: Standard Sundman ($dt = r\,d\tau$)

**Transform.** Fictitious time $\tau$ with $dt = r\,d\tau$. The radial velocity transforms as

$$\frac{dr}{d\tau} = r \cdot \dot{r} = \frac{r^2 p_r}{r - 1} \to 0 \quad \text{as } r \to 0.$$

The angular velocity transforms as

$$\frac{d\phi}{d\tau} = r \cdot \dot\phi = \frac{\ell}{r} \to \infty \quad \text{as } r \to 0.$$

**Result.** Radial dynamics regularized ($dr/d\tau \to 0$). Angular dynamics **not** regularized ($d\phi/d\tau = \ell/r \to \infty$). At $r = 10^{-8}$: $|dr/d\tau| = 10^{-5}$ vs $|d\phi/d\tau| = 10^{7}$.

**Failure mode.** The Sundman factor $g = r$ slows the radial singularity ($r^{-1/2}$ becomes $r^{1/2}$) but only partially compensates the angular one ($r^{-2}$ becomes $r^{-1}$).


### Experiment 2: Higher-order Sundman ($dt = r^\alpha\,d\tau$)

**Transform.** For general $\alpha$:

$$\frac{d\phi}{d\tau} = r^\alpha \cdot \frac{\ell}{r^2} = \ell\,r^{\alpha - 2}.$$

| $\alpha$ | $d\phi/d\tau$ as $r \to 0$ | Collision in finite $\tau$? | Winding |
|---|---|---|---|
| $< 2$ | $\to \infty$ | yes | infinite |
| $= 2$ | $= \ell$ (constant) | no ($\tau \to \infty$) | $\phi = \ell\tau \to \infty$ |
| $> 2$ | $\to 0$ | no ($\tau \to \infty$) | $\phi$ still $\to \infty$ |

**Result.** At $\alpha = 2$, all instantaneous derivatives are finite: the ODE is smooth. But $r = 0$ is a *fixed point* of $dr/d\tau = r^3 p_r/(r-1)$, so the collision takes infinite fictitious time. The total angle $\phi = \ell\,\tau \to \infty$ as $\tau \to \infty$. This is **over-regularization**: the singularity is pushed to $\tau = \infty$ rather than removed.

For $\alpha > 2$, $d\phi/d\tau \to 0$ but $\tau \to \infty$ even faster, and $\phi$ still diverges.

**Conclusion.** No Sundman exponent $\alpha$ can eliminate the infinite winding. For $\alpha < 2$, the angular velocity diverges in finite fictitious time. For $\alpha \geq 2$, derivatives are bounded, but the collision takes infinite fictitious time and the total angle is still infinite.


### Experiment 3: Levi-Civita + Sundman

**Transform.** Complex squaring $z = u^2$ with Sundman $dt = |u|^2\,d\tau$. For pure Coulomb, this yields the perfectly regular Kamiltonian

$$K_{\text{Coulomb}} = \frac{|U|^2}{8} + \kappa - E|u|^2.$$

For Weber with $\ell \neq 0$, the angular momentum term $\ell^2/(2r^2)$ transforms to

$$\frac{\ell^2}{2|u|^2}$$

after multiplication by $r = |u|^2$ (the Sundman factor). This is **still singular** at $u = 0$.

**Result.** At $|u|^2 = 10^{-6}$: the angular momentum force $|\partial K/\partial u| \sim \ell^2/|u|^3 = 9.6 \times 10^6$ (diverges). LC momenta $|U| = 2 \times 10^5$ (diverges). The LC-lifted trajectory spirals into $u = 0$ with diverging forces.

For $\ell = 0$: LC works as expected. The $\ell^2/(2|u|^2)$ term vanishes, and all quantities remain bounded.

**Failure mode.** In Kepler, $\ell \neq 0$ creates a centrifugal barrier at $r = 0$ (equivalently $u = 0$), so the $\ell^2/(2|u|^2)$ term is never tested. In Weber inside $\rho$, the inverted barrier allows reaching $u = 0$, and the LC transform cannot absorb the angular singularity.


### Experiment 4: McGehee blow-up

**Transform.** $s = \sqrt{r}$, McGehee velocity $w = s\,\dot{s}$. This does not regularize but reveals the phase-space structure near collision.

**Results.** For $\ell = 0$: $w \to -\sqrt{2}/2 \approx -0.707$ (bounded). The collision manifold $\{s = 0\}$ is invariant, and the flow extends smoothly through it. This is why the head-on collision is regularizable.

For $\ell \neq 0$: both $w$ and the angular component $\omega = \ell/s$ diverge as $s \to 0$. The collision manifold is **not** an invariant set — trajectories hit it with infinite angular velocity. No coordinate can make the manifold invariant because the angular winding is infinite.

| $\ell$ | Final $w = s\dot{s}$ | Final $\omega = \ell/s$ | Manifold invariant? |
|---|---|---|---|
| 0 | $-0.707$ | 0 | yes |
| 0.05 | $-250$ | 500 | no |
| 0.1 | $-500$ | 1000 | no |
| 0.3 | $-1500$ | 3000 | no |


### Experiment 5: Partial radial regularization

**Method.** Apply Sundman $dt = r\,d\tau$ to the radial ODE only, tracking the angle $\phi$ by quadrature:

$$\frac{d\phi}{d\tau} = \frac{\ell}{r(\tau)}.$$

**Result.** The radial trajectory $r(\tau)$ is smooth: $dr/d\tau \to 0$ as $r \to 0$ (regularized). The radial ODE can be integrated through or very close to $r = 0$ without difficulty.

The angular quadrature $\phi(\tau) = \int \ell/r\,d\tau$ diverges: near the collision, $r \sim (\tau_{\text{coll}} - \tau)^\beta$ for some $\beta$, giving $d\phi/d\tau \sim (\tau_{\text{coll}} - \tau)^{-\beta}$ which is non-integrable for $\beta \geq 1$.

**Practical utility.** Moderate. One can integrate the radial dynamics smoothly through the collision in regularized time, accepting that the angle is undefined at the collision point. This is analogous to the collision bounce approach used for $\ell = 0$ but formulated in regularized coordinates.


### Experiments 6--7: Conformal maps

**Experiment 6 (Birkhoff inversion).** $w = 1/z$ maps the collision at $z = 0$ to spatial infinity $|w| \to \infty$. The infinite spiral at the origin becomes an infinite spiral at infinity. The winding number is preserved (conformal maps preserve winding numbers). Velocity: $|dw/dt| \sim \ell/r^3 \to \infty$.

**Experiment 7 (Logarithmic map).** $w = \log z = \log r + i\phi$ "unwinds" the spiral: the trajectory goes to $(-\infty, +\infty)$ in the $w$-plane without winding. But both velocity components diverge: $dw_1/dt = \dot{r}/r \sim r^{-3/2} \to \infty$, $dw_2/dt = \dot\phi = \ell/r^2 \to \infty$. Even composing with compactification ($\arctan \circ \log$) leaves velocity divergence: $r^{-3/2}/\log^2 r \to \infty$.

**Conclusion.** The log map removes the winding but converts it to an escape-to-infinity singularity with divergent speed. No conformal or smooth coordinate change can simultaneously eliminate both the topological obstruction (winding) and the metric obstruction (speed divergence).


## Summary Table

| # | Approach | Radial at $r = 0$ | Angular at $r = 0$ | Finite $\tau$ to collision? | Through-collision? | Verdict |
|---|---|---|---|---|---|---|
| 0 | Baseline (physical time) | $\dot{r} \to \infty$ | $\dot\phi \to \infty$ | N/A (finite $t$) | no | -- |
| 1 | Sundman $dt = r\,d\tau$ | $\to 0$ | $\ell/r \to \infty$ | yes | radial only | partial |
| 2 | Sundman $dt = r^2\,d\tau$ | $\to 0$ | $\ell$ (const) | **no** ($\tau \to \infty$) | never reaches | over-regularized |
| 3 | Levi-Civita + Sundman | smooth | $\ell^2/|u|^2 \to \infty$ | yes | no | fails |
| 4 | McGehee blow-up | $w$ bounded ($\ell=0$) | $\ell/s \to \infty$ | N/A | $\ell=0$ only | diagnostic |
| 5 | Partial radial | $\to 0$ | log divergence | yes | radial yes, angular no | moderate utility |
| 6 | Birkhoff $w = 1/z$ | maps to $\infty$ | winding preserved | N/A | no | fails |
| 7 | Log map $w = \log z$ | $\to -\infty$ | unwound but $\dot{w} \to \infty$ | N/A | no | fails |

## Why the Head-On Collision ($\ell = 0$) IS Regularizable

For $\ell = 0$, the dynamics reduce to a 1D radial problem with no angular component. The collision speed is finite ($\sqrt{2}\,c$), and there is no winding. The standard Sundman transform $dt = r\,d\tau$ reduces the velocity to zero in fictitious time, and the trajectory passes through $r = 0$ smoothly:

- **McGehee velocity** $w = s\dot{s} \to -\sqrt{2}/2$ (bounded), collision manifold is invariant
- **Levi-Civita** works because the $\ell^2/(2|u|^2)$ singular term is absent
- **Collision bounce** (reflection $q_{\text{rel}} \to -q_{\text{rel}}$) is the practical implementation used in the WeberElectrodynamics package

The regularizability of the $\ell = 0$ case versus the non-regularizability of $\ell \neq 0$ is a clean dichotomy with no intermediate cases.

## Connection to Quantum Theory

Frauenfelder & Weber (2024, Theorem A) show that the Weber--Schr\"{o}dinger equation on the interval $(0, \rho)$ is a singular Sturm--Liouville problem that is limit circle at both endpoints. The proof reveals a classical--quantum parallel:

| Angular momentum | Classical | Quantum |
|---|---|---|
| $\ell = 0$ | collision (regularizable) | non-oscillating (natural BCs exist) |
| $\ell \neq 0$ | spiraling (not regularizable) | wildly oscillating (no natural BCs) |

The classical regularity that permits regularization at $\ell = 0$ is mirrored by the existence of natural boundary conditions in the quantum theory. Conversely, the classical non-regularizability at $\ell \neq 0$ corresponds to wildly oscillating quantum solutions for which no canonical boundary condition exists.

## Practical Implications

1. **Collision bounce** remains the best practical approach for the $\ell = 0$ sub-critical like-charge oscillation in the WeberElectrodynamics package (see [CollisionBounceRegularization.md](../exploratory/CollisionBounceRegularization.md)).

2. For $\ell \neq 0$ inside $\rho$, **no regularization exists**. The integrator must either:
   - Avoid initial conditions that lead to spiraling collisions
   - Stop integration before the collision (the trajectory is well-defined up to any finite time before collision)
   - Use partial radial regularization to extend the radial dynamics through collision, accepting that the angle is undefined at the collision point

3. In multi-body simulations, trajectories that develop $\ell \neq 0$ inside $\rho$ will inevitably fail. The timescale to collision can be estimated from the energy equation, providing advance warning.

## References

- Frauenfelder, U., Weber, J. "A mathematical description of the Weber nucleus as a classical and quantum mechanical system." *Anal. Math. Phys.* **14**:31 (2024). DOI: [10.1007/s13324-024-00891-5](https://doi.org/10.1007/s13324-024-00891-5).
- McGehee, R. "Triple collision in the collinear three-body problem." *Inventiones math.* **27**, 191--227 (1974).
- Sundman, K. F. "M\'{e}moire sur le probl\`{e}me des trois corps." *Acta Math.* **36**, 105--179 (1913).

See also: [CriticalRadiusAndLikeChargeAttraction.md](CriticalRadiusAndLikeChargeAttraction.md), [Regularization.md](../theory/Regularization.md), [WeberElectrodynamics.md](../theory/WeberElectrodynamics.md).
