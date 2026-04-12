# Canonical Two-Body Reference Problem

This document specifies the reference two-body Weber problem used as the
canonical hello-world for this package. It is implementation-independent:
any integrator of the Weber Hamiltonian should be able to reproduce the
expected observables listed below to the stated tolerances.

The companion artifacts are:

- [`two_body_reference.ipynb`](two_body_reference.ipynb) — annotated notebook tutorial.
- [`two_body_reference.jl`](two_body_reference.jl) — script that regenerates every figure to `output/`.

## Hamiltonian

Two charged point particles in a plane, evolved under the Weber Hamiltonian
(Gauss–Weber absolute units):

$$H = \frac{|\mathbf{p}_1|^2}{2 m_1} + \frac{|\mathbf{p}_2|^2}{2 m_2} + \frac{q_1 q_2}{r_{12}} \left(1 - \frac{\dot{r}_{12}^2}{2 c^2}\right)$$

with $r_{12} = |\mathbf{r}_1 - \mathbf{r}_2|$ and
$\dot{r}_{12} = (\mathbf{r}_1 - \mathbf{r}_2) \cdot (\dot{\mathbf{r}}_1 - \dot{\mathbf{r}}_2) / r_{12}$.
See [`theory/WeberElectrodynamics.md`](../theory/WeberElectrodynamics.md) for the full derivation.

## Parameters

All values are exact rationals; every parameter is a whole number.

| Symbol     | Value        | Meaning                           |
|------------|--------------|-----------------------------------|
| $m_1, m_2$ | $1, 1$       | Particle masses                   |
| $q_1, q_2$ | $+1, -1$     | Charges (unlike → attractive)     |
| $c$        | $4$          | Speed of light                    |
| $r_0$      | $2$          | Initial separation                |
| Dimension  | $2$          | Planar motion                     |

Derived from the above, with $k = q_1 q_2 = -1$ and reduced mass $\mu = m_1 m_2 / (m_1 + m_2) = 1/2$:

| Quantity             | Expression                       | Value           |
|----------------------|----------------------------------|-----------------|
| Circular relative $v$| $\sqrt{\lvert k\rvert / (\mu r_0)}$ | $1$             |
| Chosen relative $v$  | $v_{\text{rel,circ}} / 2$        | $1/2$           |
| Per-particle speed   | $v_{\text{rel}} \cdot m_j / M$   | $1/4$           |

The relative velocity is exactly half the circular value. This is deliberate:
it produces a large but closed Keplerian ellipse in the Coulomb limit so that
the Weber-induced precession is unambiguously visible against the background
motion.

## Initial conditions

The configuration is symmetric about the origin, with the centre of mass at
rest at the origin. Initial momenta are perpendicular to the separation
vector, so $\dot{r}_{12}(0) = 0$ exactly (see the [zero-radial-velocity
principle](../theory/InitialConditions.md#the-zero-radial-velocity-principle)).
Consequently the Weber potential reduces to the Coulomb potential at $t=0$
and the initial energy is analytic.

$$\mathbf{r}_1(0) = (-1, 0), \qquad \mathbf{r}_2(0) = (+1, 0)$$

$$\mathbf{p}_1(0) = (0, -1/4), \qquad \mathbf{p}_2(0) = (0, +1/4)$$

As flat state vectors:

```
q(0) = [-1,  0,  +1,  0]
p(0) = [ 0, -1/4, 0, +1/4]
```

## Derived analytical quantities

Because $\dot{r}_{12}(0) = 0$, the initial state lies on a Keplerian ellipse
of the same $k$ and $\mu$. The Coulomb two-body formulas apply at $t=0$:

| Quantity                   | Expression                                 | Value           |
|----------------------------|--------------------------------------------|-----------------|
| Total energy $E_0$         | $\tfrac{1}{2}\mu v_{\text{rel}}^2 + k/r_0$ | $-7/16 = -0.4375$ |
| Total angular momentum $L_z$ | $\sum_i (x_i p_{y,i} - y_i p_{x,i})$     | $+1/2$          |
| Total linear momentum $\mathbf{P}$ | $\sum_i \mathbf{p}_i$                | $\mathbf{0}$    |
| Semi-major axis $a$        | $-k / (2 E_0)$                             | $8/7 \approx 1.1429$ |
| Eccentricity $e$           | $\sqrt{1 + 2 E_0 L_z^2 / (\mu k^2)}$       | $3/4 = 0.75$    |
| Apoapsis $r_{\max}$        | $a(1 + e)$                                 | $2$ (equals $r_0$) |
| Periapsis $r_{\min}$       | $a(1 - e)$                                 | $2/7 \approx 0.2857$ |
| Kepler period $T$          | $2\pi \sqrt{a^3 \mu / |k|}$                | $\approx 5.428$ |

The starting configuration is at apoapsis. A canonical run covers
$5 T \approx 27.14$ time units — five Keplerian periods, long enough for the
Weber precession to produce a clearly visible rosette pattern and short
enough for the integration to finish in a fraction of a second.

## Expected qualitative behaviour

- Bounded orbital motion for all time.
- In the pure-Coulomb limit ($c \to \infty$) the orbit would close after one
  Kepler period. The finite $c = 4$ introduces the velocity-dependent term,
  which advances the periapsis by a small angle per orbit (a rosette).
- Separation oscillates between $r_{\min} \approx 0.2857$ and $r_{\max} = 2$.
- Radial velocity $\dot{r}_{12}$ is zero at both turning points and
  maximum in magnitude between them; its square is the Weber correction.

## Acceptance tolerances

These bounds characterise a "correct" integration of this specific problem
and are the ones satisfied by the reference implementation in this repository
with $dt = 10^{-3}$ over $5 T$:

| Quantity                         | Bound                |
|----------------------------------|----------------------|
| Maximum relative energy drift    | $< 10^{-5}\,\%$      |
| Maximum $|\mathbf{P}|$ drift     | machine precision    |
| Maximum $|L_z - L_z(0)|$         | $< 10^{-10}$         |
| Symbolic-Hamiltonian self-consistency | $< 10^{-15}$    |

Any symplectic integrator with comparable step size and regularization policy
should meet these. A non-symplectic integrator will typically violate the
energy bound long before it violates angular momentum.

## Cross-references

- [`theory/WeberElectrodynamics.md`](../theory/WeberElectrodynamics.md) — Weber Hamiltonian and equations of motion.
- [`theory/InitialConditions.md`](../theory/InitialConditions.md) — general construction of consistent initial conditions.
- [`theory/SemiExplicitIntegrator.md`](../theory/SemiExplicitIntegrator.md) — integrator theory.
