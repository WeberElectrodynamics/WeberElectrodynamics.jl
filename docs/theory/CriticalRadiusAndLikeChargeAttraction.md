# Critical Radius and Like-Charge Attraction in Weber Electrodynamics

## Scope

This document describes the sub-critical regime of Weber's force law in which two like-charged particles attract rather than repel. The velocity- and acceleration-dependent terms in Weber's force define a critical separation $\rho$ that partitions the dynamics into two permanently distinct states: ordinary Coulomb-like repulsion at large distance, and bound oscillatory attraction at molecular distance. Weber developed this theory in his Sixth Memoir (1871) and refined it in the Seventh (1878) and Eighth (posthumous, 1894).

## Notation

We consider two charged particles with charges $e, e'$ and inertial masses $\varepsilon, \varepsilon'$. The reduced mass is

$$\mu = \frac{\varepsilon\varepsilon'}{\varepsilon + \varepsilon'}.$$

Separation $r$, radial velocity $\dot{r} = dr/dt$, radial acceleration $\ddot{r} = d^2r/dt^2$. The constant $c$ is Weber's electrodynamic constant, related to the speed of light by $c_W = \sqrt{2}\,c_{\text{light}}$. For full notation see [WeberElectrodynamics.md](WeberElectrodynamics.md).

## Weber Force and Effective Inertial Mass

The radial Weber force between $e$ and $e'$ is (see [WeberElectrodynamics.md](WeberElectrodynamics.md)):

$$F = \frac{ee'}{r^2}\left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right).$$

Newton's second law for the relative coordinate gives

$$\mu\,\ddot{r} = F + \text{(centrifugal terms)}.$$

Substituting the Weber force and isolating $\ddot{r}$:

$$\mu\,\ddot{r} = \frac{ee'}{r^2}\left(1 - \frac{\dot{r}^2}{2c^2}\right) + \frac{ee'}{rc^2}\,\ddot{r}.$$

Collecting the acceleration on the left:

$$\mu\left(1 - \frac{ee'}{\mu r c^2}\right)\ddot{r} = \frac{ee'}{r^2}\left(1 - \frac{\dot{r}^2}{2c^2}\right).$$

Define the critical radius

$$\rho = \frac{ee'}{\mu c^2} = \frac{(\varepsilon + \varepsilon')}{\varepsilon\varepsilon'}\cdot\frac{ee'}{c^2}$$

so that the equation of motion becomes

$$\mu\left(1 - \frac{\rho}{r}\right)\ddot{r} = \frac{ee'}{r^2}\left(1 - \frac{\dot{r}^2}{2c^2}\right).$$

The factor multiplying $\ddot{r}$ is the effective inertial mass:

$$\mu_{\text{eff}}(r) = \mu\left(1 - \frac{\rho}{r}\right).$$

For like charges ($ee' > 0$), $\rho > 0$ and $\mu_{\text{eff}}$ changes sign at $r = \rho$.

## The Critical Radius

### Definition

Weber defines (Sixth Memoir, §9.9):

$$\rho = 2\left(\frac{1}{\varepsilon} + \frac{1}{\varepsilon'}\right)\frac{ee'}{c^2} = \frac{2ee'}{\mu c^2}.$$

The factor of 2 compared to the effective-mass derivation above reflects Weber's convention for the radial force law; both forms yield the same critical distance.

### Sign dependence

- **Like charges** ($ee' > 0$): $\rho > 0$, a real positive distance. The dynamics bifurcate at $r = \rho$.
- **Unlike charges** ($ee' < 0$): $\rho < 0$, no positive critical distance exists. The force remains attractive at all separations and no sign flip occurs.

### Physical scale

Weber estimated (§9.14) that for oscillation frequencies matching visible light, $\rho$ must lie between $1/4000$ and $1/8000$ of a millimetre, i.e. on the order of $10^{-7}$ m. This places the critical radius firmly in the molecular domain — far below any macroscopically measurable distance.

### General force law

In the Seventh Memoir (§13.5), Weber derived the general expression for the repulsive force including external acceleration $f$:

$$-\frac{dV}{dr} = \frac{ee'}{r(r - \rho)}\left(1 - \frac{\dot{r}^2}{c^2} + \frac{2r}{c^2}f\right).$$

The denominator $r(r - \rho)$ reveals the structure: a singularity at $r = \rho$ where the force diverges, and a sign change across it.

## Two States of Aggregation

For two like particles ($\rho > 0$) with purely radial motion and $\dot{r} = 0$ at $r = r_0$, Weber's integrated equation of motion (§9.9) is

$$\frac{\dot{r}^2}{c^2} = \frac{r - r_0}{r - \rho}\cdot\frac{\rho}{r_0}.$$

The requirement $\dot{r}^2 \geq 0$ partitions the motion into two regimes depending on where $r_0$ lies relative to $\rho$:

### Distant state ($r_0 > \rho$)

The right-hand side is non-negative only for $r \geq r_0$. The particles approach from infinity, decelerate, reach closest approach $r_0 > \rho$, and recede. This is ordinary Coulomb-like scattering. The maximum velocity satisfies $\dot{r}^2 \leq (\rho/r_0)c^2 < c^2$.

### Molecular state ($r_0 < \rho$)

The right-hand side is non-negative only for $r \leq r_0$. The particles oscillate between $r = r_0$ and $r = 0$, never reaching $\rho$ from below. As $r \to 0$, the velocity approaches $\dot{r}^2 \to c^2$. Weber called this a "molecular movement."

### Impenetrability of the barrier

No continuous trajectory can cross $r = \rho$. At that point $\dot{r}^2$ would need to diverge or pass through imaginary values. Weber (§9.10):

> No transition from one state of aggregation to the other takes place so long as both particles move in consequence of their reciprocal action only.

The critical radius is an absolute dynamical barrier for an isolated pair.

## The Oscillatory Atomic Pair

Two like charges in the molecular state ($r_0 < \rho$) form what Weber called an "electrical atomic pair" (§9.14). The pair oscillates between $r_0$ and $r = 0$ with a well-defined period.

### Period of oscillation

For small amplitude ($r_0 \ll \rho$) and vanishing transverse velocity ($\alpha_0 \to 0$):

$$2\vartheta = \frac{4r_0}{c}.$$

The period is proportional to the amplitude $r_0$, with coefficient $4/c$ for small oscillations, decreasing to $2/c$ at the maximum amplitude $r_0 = \rho$.

### Transverse motion

When the particles also have a relative transverse velocity $\alpha_0$, the equation of motion generalises (§9.11) to

$$\frac{\dot{r}^2}{c^2} = \frac{r - r_0}{r - \rho}\cdot\frac{\rho}{r_0} + \frac{r + r_0}{r}\cdot\frac{\alpha_0^2}{c^2}.$$

The same two-state structure persists. However, no stable circular orbit exists for like charges below $\rho$ (§9.13): the attraction always exceeds the centripetal requirement, so any transverse motion still results in radial oscillation with inward spiralling tendency.

## Physical Interpretation: Negative Effective Mass

The sign flip at $r = \rho$ can be understood through the effective inertial mass:

| Regime | $r$ vs $\rho$ | $\mu_{\text{eff}}$ | Force character | Dynamics |
| --- | --- | --- | --- | --- |
| Distant | $r > \rho$ | positive | repulsive | scattering |
| Critical | $r = \rho$ | zero (singular) | $F \to \infty$ | barrier |
| Molecular | $r < \rho$ | negative | attractive | bound oscillation |

Below $\rho$, the Coulomb factor $ee'/r^2$ is still formally positive (repulsive). But the negative effective mass inverts the dynamical response: the particle accelerates *toward* the source of repulsion rather than away from it. This is not a change in the force but a change in how inertia responds to it.

Equivalently, the velocity-dependent Weber potential

$$U = \frac{ee'}{r}\left(1 - \frac{\dot{r}^2}{2c^2}\right)$$

creates, for $r < \rho$, a confining well in the combined position-velocity phase space. The kinetic coupling $\dot{r}^2/c^2$ grows as the particles approach, lowering the effective potential faster than the Coulomb barrier rises.

## Connection to Zöllner Gravitation and the Planetary Atom

Weber's sub-critical binding of like charges is the dynamical foundation for the program described in [ZollnerElectrogravitationalTheory.md](ZollnerElectrogravitationalTheory.md).

In the Eighth Memoir (§15.1), Weber proposed that all ponderable matter consists of bound pairs of positive and negative electrical molecules. The molecular state provides the mechanism: like charges can be permanently bound at sub-$\rho$ separations, forming a stable nucleus. Unlike charges then orbit this nucleus under the combined electrodynamic and gravitational residual forces.

Zöllner's mismatch parameter $\alpha$ and the critical radius $\rho$ play complementary roles:

- $\alpha$ produces the gravitational residual (a small net attraction between neutral dyads).
- $\rho$ enables the internal binding that holds the dyad together in the first place.

Without the sub-critical regime, there is no mechanism for like charges to cohere, and the entire electrogravitational program loses its atomic foundation.

## References

- Weber, W. E. "Electrodynamic Measurements, Sixth Memoir" (1871). In Assis (ed.), *Wilhelm Weber's Main Works on Electrodynamics*, Vol. IV, Ch. 9, §§9.8--9.17.
- Weber, W. E. "Electrodynamic Measurements, Seventh Memoir" (1878). In Assis (ed.), Vol. IV, Ch. 13, §§13.5--13.6.
- Weber, W. E. "Electrodynamic Measurements, Eighth Memoir" (1894, posthumous). In Assis (ed.), Vol. IV, Ch. 15, §§15.1--15.2.
- Assis, A. K. T. (ed.) *Wilhelm Weber's Main Works on Electrodynamics Translated into English, Volume IV: Conservation of Energy, Weber's Planetary Model of the Atom and the Unification of Electromagnetism and Gravitation*. Apeiron, Montreal, 2021. ISBN 978-1-987980-29-5.
- Assis, A. K. T., Tajmar, M. "Superluminal potentials, forces and accelerations in Weber electrodynamics." *Journal of Advanced Physics* 8 (2019).

See also: [WeberElectrodynamics.md](WeberElectrodynamics.md), [ZollnerElectrogravitationalTheory.md](ZollnerElectrogravitationalTheory.md), [Regularization.md](Regularization.md).
