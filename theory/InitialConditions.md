# Initial Conditions for Weber Electrodynamic Systems

## Scope

This document presents systematic methods for constructing initial conditions for $n$-body Weber electrodynamic systems. Given a target total energy and a choice of geometry or symmetry class, the methods here yield explicit positions and momenta in closed form. The treatment covers two-body orbits, symmetric planar polygon configurations, three-dimensional symmetric polyhedra, and the Zöllner-modified case. All derivations use absolute (Gauss–Weber) units.

## Notation

We follow the conventions established for the Weber Hamiltonian. For an $n$-particle system:

- $m_i$, $q_i$ — mass and electric charge of particle $i$
- $\vec{r}_i = (x_i, y_i, z_i)$ — position vector
- $\vec{v}_i = \dot{\vec{r}}_i$ — velocity vector
- $\vec{p}_i = m_i \vec{v}_i$ — canonical momentum
- $M = \sum_i m_i$ — total mass
- $\kappa_{ij} \geq 1$ — pair coupling factor (equals 1 in standard Weber; see Section 7 for the Zöllner generalisation)

For each pair $(i,j)$ with $i < j$:

$$\vec{r}_{ij} = \vec{r}_i - \vec{r}_j, \qquad r_{ij} = |\vec{r}_{ij}|, \qquad \hat{r}_{ij} = \frac{\vec{r}_{ij}}{r_{ij}}$$

$$\dot{r}_{ij} = \frac{\vec{r}_{ij} \cdot \dot{\vec{r}}_{ij}}{r_{ij}}$$

The Weber Hamiltonian for $n$ particles is

$$H = \sum_{i=1}^n \frac{|\vec{p}_i|^2}{2m_i} + \sum_{i < j} U_{ij}$$

with pair potential

$$U_{ij} = \frac{\kappa_{ij} q_i q_j}{r_{ij}} \left(1 - \frac{\dot{r}_{ij}^2}{2c^2}\right)$$

where $c$ is the speed of light. The initial value problem is to choose $\vec{r}_i(0)$ and $\vec{p}_i(0)$ such that $H(0)$ takes a prescribed value $E$.

## The Zero-Radial-Velocity Principle

### Simplification

The Weber potential depends on the radial velocity $\dot{r}_{ij}$, making $H$ depend on momenta in a nonlinear way. This complicates the direct inversion of the energy equation. The central simplification is:

**If $\dot{r}_{ij}(0) = 0$ for every pair $(i,j)$, the Weber potential reduces to the Coulomb potential at $t = 0$.**

With this condition, the initial potential energy is

$$U_0 = \sum_{i < j} \frac{\kappa_{ij} q_i q_j}{r_{ij}(0)}$$

and the energy constraint becomes the explicit relation

$$T_0 = E - U_0, \qquad T_0 = \sum_{i=1}^n \frac{|\vec{p}_i(0)|^2}{2m_i}$$

The required kinetic energy $T_0$ must be positive, which imposes $E > U_0$.

### How to Satisfy $\dot{r}_{ij}(0) = 0$

The condition $\dot{r}_{ij} = 0$ means the rate of change of the separation $r_{ij}$ is zero: particles are neither approaching nor receding at $t = 0$. It is satisfied whenever the instantaneous velocity of every particle is perpendicular to all separation vectors connecting it to other particles. Two constructions cover all cases considered in this document:

**Tangential momenta (2D).** Place all particles in a plane. Assign each particle a momentum directed perpendicular to the position vector in that plane (i.e.,tangential to the circle on which it sits). For any pair, the relative velocity is then tangential to the line connecting the two particles, so $\vec{r}_{ij} \cdot \dot{\vec{r}}_{ij} = 0$.

**Rigid rotation (3D).** For any configuration rotating rigidly with angular velocity $\vec{\omega}$, set

$$\vec{p}_i(0) = m_i \left(\vec{\omega} \times \vec{r}_i(0)\right)$$

Because rigid rotation preserves all pairwise distances, $d\,r_{ij}/dt = 0$ exactly for every pair. This construction reduces to tangential momenta in 2D when $\vec{\omega} = \omega \hat{z}$.

## Conservation Laws and the Centre-of-Mass Frame

Three quantities are conserved by the Weber Hamiltonian for an isolated system:

$$E = H(\vec{r}(t), \vec{p}(t)) = \text{const}$$

$$\vec{P} = \sum_{i=1}^n \vec{p}_i = \text{const}$$

$$\vec{L} = \sum_{i=1}^n \vec{r}_i \times \vec{p}_i = \text{const}$$

Working in the centre-of-mass (COM) frame removes two free parameters. The COM frame conditions at $t = 0$ are:

$$\sum_{i=1}^n m_i \vec{r}_i(0) = \vec{0}, \qquad \sum_{i=1}^n \vec{p}_i(0) = \vec{0}$$

Both tangential-momentum and rigid-rotation constructions automatically satisfy $\sum \vec{p}_i = 0$ when the COM is at the origin, which follows from:

$$\sum_i m_i (\vec{\omega} \times \vec{r}_i) = \vec{\omega} \times \left(\sum_i m_i \vec{r}_i\right) = \vec{\omega} \times \vec{0} = \vec{0}$$

## Method A — Two-Body Bound Orbits

### Centre-of-Mass Positions

For two particles with masses $m_1$, $m_2$ and total mass $M = m_1 + m_2$, place them along the $x$-axis with separation $r_0 > 0$:

$$\vec{r}_1(0) = -\frac{m_2}{M} r_0 \,\hat{x}, \qquad \vec{r}_2(0) = +\frac{m_1}{M} r_0 \,\hat{x}$$

This satisfies $m_1 \vec{r}_1 + m_2 \vec{r}_2 = \vec{0}$.

Tangential momenta (perpendicular to $\hat{x}$, ensuring $\dot{r}_{12}(0) = 0$) take the form

$$\vec{p}_1(0) = +p_\perp \,\hat{y}, \qquad \vec{p}_2(0) = -p_\perp \,\hat{y}$$

for some scalar $p_\perp > 0$. This satisfies $\vec{p}_1 + \vec{p}_2 = \vec{0}$ and $\dot{r}_{12}(0) = 0$.

### Initial Energy

With $\dot{r}_{12}(0) = 0$ the initial Hamiltonian is

$$H(0) = \frac{p_\perp^2}{2m_1} + \frac{p_\perp^2}{2m_2} + \frac{\kappa_{12} q_1 q_2}{r_0} = \frac{p_\perp^2}{2\mu} + \frac{k}{r_0}$$

where $\mu = m_1 m_2 / M$ is the reduced mass and $k = \kappa_{12} q_1 q_2$. The energy constraint becomes

$$E = \frac{p_\perp^2}{2\mu} + \frac{k}{r_0}$$

### Circular Orbit Velocity

For an attractive pair ($k < 0$), a circular orbit at separation $r_0$ requires force balance between the Coulomb attraction and the centripetal acceleration. The radial equation in the reduced-mass frame is

$$F_r = \frac{\mu v^2}{r_0}$$

where $v = p_\perp / \mu$ is the speed in the relative frame.

The Weber force at circular orbit reduces exactly to the Coulomb force. This can be verified in both standard formulations. In the explicit time-parameterisation,

$$F_r^{\mathrm{Weber}} = \frac{|k|}{r_0^2}\left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right)$$

For a circular orbit $r = r_0 = \text{const}$, so $\dot{r} = 0$ and $\ddot{r} = 0$, giving $F_r^{\mathrm{Weber}} = |k|/r_0^2$. In the vector formulation,

$$F_r^{\mathrm{Weber}} = \frac{|k|}{r_0^2}\left(1 + \frac{\vec{v} \cdot \vec{v} + \vec{r} \cdot \vec{a} - \tfrac{3}{2}(\hat{r} \cdot \vec{v})^2}{c^2}\right)$$

For circular motion $\hat{r} \cdot \vec{v} = 0$ (motion is tangential), $\vec{v} \cdot \vec{v} = v^2$, and $\vec{r} \cdot \vec{a} = -v^2$ (centripetal), so the correction is $(v^2 - v^2 - 0)/c^2 = 0$. The Weber force equals the Coulomb force exactly for all circular orbits, regardless of speed.

Force balance then gives the circular orbit velocity:

$$v_{\mathrm{circ}} = \sqrt{\frac{|k|}{\mu r_0}}$$

This is identical to the Coulomb result. In terms of the canonical momentum magnitude,

$$p_{\mathrm{circ}} = \mu v_{\mathrm{circ}} = \sqrt{\mu |k| / r_0}$$

The circular orbit energy follows from $H(0)$ with $p_\perp = p_{\mathrm{circ}}$ and $\dot{r} = 0$:

$$E_{\mathrm{circ}} = \frac{p_{\mathrm{circ}}^2}{2\mu} + \frac{k}{r_0} = \frac{|k|}{2r_0} - \frac{|k|}{r_0} = -\frac{|k|}{2r_0}$$

This is the Weber virial relation at circular orbit, and it is identical to the Coulomb result because $\dot{r} = 0$ at all times for a circular orbit.

### Elliptical Orbits via Velocity Scale

Introduce a dimensionless velocity scale $\eta_v > 0$ and set

$$p_\perp = \eta_v \, p_{\mathrm{circ}} = \eta_v \sqrt{\mu |k| / r_0}$$

The initial energy as a function of $\eta_v$ is

$$E(\eta_v) = \frac{\eta_v^2 p_{\mathrm{circ}}^2}{2\mu} + \frac{k}{r_0} = \frac{\eta_v^2 |k|}{2r_0} - \frac{|k|}{r_0} = \frac{|k|}{r_0}\!\left(\frac{\eta_v^2}{2} - 1\right)$$

Key values:

| $\eta_v$ | Orbit type | Energy |
| --- | --- | --- |
| $0$ | radial free fall | $-|k|/r_0$ |
| $(0, 1)$ | elliptical, sub-circular | $E < E_{\mathrm{circ}}$ |
| $1$ | circular | $E_{\mathrm{circ}} = -|k|/(2r_0)$ |
| $(1, \sqrt{2})$ | elliptical, super-circular | $E_{\mathrm{circ}} < E < 0$ |
| $\sqrt{2}$ | parabolic escape | $E = 0$ |
| $> \sqrt{2}$ | hyperbolic (unbound) | $E > 0$ |

The orbit is bound if and only if $E < 0$, i.e.,$\eta_v < \sqrt{2}$.

**Inverse problem.** Given target energy $E < 0$ and separation $r_0$, the required velocity scale is

$$\eta_v = \sqrt{2\left(1 + \frac{E r_0}{|k|}\right)}$$

This requires $|E| \leq |k|/r_0$, i.e.,the target energy cannot be more negative than the purely radial-fall energy.

### Method from Energy and Angular Momentum

An alternative parameterisation uses the total energy $E$ and the angular momentum magnitude $L = |\vec{L}|$. At a turning point where $\dot{r}_{12} = 0$ (periapsis or apoapsis), the Weber potential reduces to Coulomb and the effective one-body energy equation reads

$$E = \frac{L^2}{2\mu r^2} + \frac{k}{r}$$

Substituting $u = 1/r$ this becomes the quadratic

$$\frac{L^2}{2\mu} u^2 + k\, u - E = 0$$

The two positive roots for an attractive bound orbit ($k < 0$, $E < 0$) are

$$u_{p,a} = \frac{-k \pm \sqrt{k^2 + 2EL^2/\mu}}{L^2/\mu}$$

giving periapsis $r_p = 1/u_p$ (smaller) and apoapsis $r_a = 1/u_a$ (larger). The discriminant is positive when

$$L^2 < \frac{\mu k^2}{2|E|}$$

The eccentricity of the Keplerian reference orbit is

$$e = \frac{r_a - r_p}{r_a + r_p} = \sqrt{1 + \frac{2E L^2}{\mu k^2}}$$

The tangential speed at apoapsis (where $\dot{r} = 0$ and the initial condition is most conveniently set) is $v_a = L / (\mu r_a)$, giving momentum magnitude $p_\perp = \mu v_a = L / r_a$.

**Note.** The orbit will not be a closed Keplerian ellipse under Weber dynamics; once radial motion develops the velocity-dependent terms activate and cause orbital precession. The Keplerian elements $E$, $L$, $r_p$, $r_a$ characterise the orbit only at the chosen initial instant.

### When Weber Corrections Activate

The result that circular orbits are unaffected by Weber corrections extends to any initial instant with $\dot{r}_{ij}(0) = 0$: the Hamiltonian value, and therefore $T_0 = E - U_0$, is computed entirely from the Coulomb potential. Weber velocity-dependent terms enter only once radial motion develops ($\dot{r} \neq 0$), typically within the first quarter of an orbital period. Their principal dynamical effects are:

- **Orbital precession.** The radial restoring force is modified by a velocity-dependent term proportional to $\dot{r}^2/c^2$, leading to apsidal advance or recession depending on orbit shape.
- **Energy redistribution.** Kinetic and potential energy oscillate differently from the Keplerian case; the total $H$ remains conserved.

These are dynamical consequences of the initial conditions; they do not alter the initial condition calculation itself.

## Method B — Symmetric N-Body Configurations (Planar)

### The Energy-Fraction Recipe

Consider $N$ particles of equal mass $m$ placed at the vertices of a regular $N$-gon of circumradius $R$ in the $xy$-plane. Assign charges $\pm Q$ in alternating order around the polygon. For even $N$, adjacent vertices always carry opposite charges, which maximises the number of unlike (attractive) pairs.

Place particle $i$ at

$$\vec{r}_i(0) = R\left(\cos\!\frac{2\pi(i-1)}{N},\;\sin\!\frac{2\pi(i-1)}{N},\; 0\right), \qquad i = 1, \ldots, N$$

and assign charges $q_i = (-1)^{i+1} Q$.

Assign momenta tangentially in the counter-clockwise sense:

$$\vec{p}_i(0) = m v\left(-\sin\!\frac{2\pi(i-1)}{N},\;\cos\!\frac{2\pi(i-1)}{N},\; 0\right)$$

for a common speed $v$ to be determined. By symmetry $\sum_i \vec{p}_i = \vec{0}$ and $\dot{r}_{ij}(0) = 0$ for all pairs.

The initial potential energy $U_0$ depends only on the geometry and charges. Define the energy fraction

$$\eta = \frac{T_0}{|U_0|}, \qquad T_0 = \frac{1}{2} N m v^2$$

For attractive $U_0 < 0$ (which holds for all alternating-charge polygons with $N \geq 2$), the total energy is

$$E = T_0 + U_0 = (\eta - 1)|U_0|$$

and the orbit is bound if and only if $\eta < 1$.

Given a target $\eta$, the common particle speed is

$$v = \sqrt{\frac{2\eta |U_0|}{Nm}}$$

Equivalently, given a target $E < 0$:

$$v = \sqrt{\frac{2(E - U_0)}{Nm}}$$

### Computing $U_0$ for a Regular Alternating-Charge Polygon

For even $N$, the vertex-to-vertex distance between particles $k$ steps apart around the polygon is

$$d_k = 2R \sin\!\frac{k\pi}{N}, \qquad k = 1, \ldots, N/2$$

The charge product for a pair $k$ steps apart is $(-1)^k Q^2$: negative (unlike, attractive) when $k$ is odd; positive (like, repulsive) when $k$ is even. There are $N$ pairs for each $k < N/2$ and $N/2$ pairs for $k = N/2$.

The initial potential energy is therefore

$$U_0 = \frac{NQ^2}{2R}\left[\sum_{k=1}^{N/2-1} \frac{(-1)^k}{\sin(k\pi/N)} + \frac{(-1)^{N/2}}{2}\right]$$

For all even $N$ with alternating charges, $U_0 < 0$ because the sum is dominated by the $k = 1$ adjacent unlike pairs.

### N = 2: Two-Body Dimer

With $N = 2$, $R = r_0/2$, and a single unlike pair at $r_0 = 2R$:

$$U_0 = \frac{-Q^2}{r_0}$$

Setting $q_1 q_2 = -Q^2$ and $\kappa_{12} = 1$, this reproduces $k/r_0$ with $k = -Q^2 < 0$. The polygon method gives each particle speed $v = \sqrt{\eta Q^2/(m r_0)}$; the two-body reduced-mass circular speed (relative frame) is $v_{\mathrm{circ}} = \sqrt{2Q^2/(m r_0)}$ for $\mu = m/2$. Since for equal masses the relative speed is twice the individual particle speed, the two parameterisations are related by $\eta_v = \sqrt{2\eta}$, or equivalently the circular orbit ($\eta_v = 1$) corresponds to $\eta = 1/2$.

### N = 4: Alternating Square

Place four particles at the cardinal points of circumradius $R$:

$$\vec{r}_1(0) = (R, 0), \quad \vec{r}_2(0) = (0, R), \quad \vec{r}_3(0) = (-R, 0), \quad \vec{r}_4(0) = (0, -R)$$

with charges $q_1 = +Q$, $q_2 = -Q$, $q_3 = +Q$, $q_4 = -Q$.

Pairs and distances:

| Pair | Charges | Distance | Type | Contribution to $U_0$ |
| --- | --- | --- | --- | --- |
| $(1,2)$ | $(+Q,{-Q})$ | $R\sqrt{2}$ | unlike | $-Q^2/(R\sqrt{2})$ |
| $(2,3)$ | $(-Q,{+Q})$ | $R\sqrt{2}$ | unlike | $-Q^2/(R\sqrt{2})$ |
| $(3,4)$ | $(+Q,{-Q})$ | $R\sqrt{2}$ | unlike | $-Q^2/(R\sqrt{2})$ |
| $(4,1)$ | $(-Q,{+Q})$ | $R\sqrt{2}$ | unlike | $-Q^2/(R\sqrt{2})$ |
| $(1,3)$ | $(+Q,{+Q})$ | $2R$ | like | $+Q^2/(2R)$ |
| $(2,4)$ | $(-Q,{-Q})$ | $2R$ | like | $+Q^2/(2R)$ |

Summing:

$$U_0 = -\frac{4Q^2}{R\sqrt{2}} + \frac{2Q^2}{2R} = \frac{Q^2}{R}\!\left(1 - 2\sqrt{2}\right)$$

Since $2\sqrt{2} \approx 2.828 > 1$, we have $U_0 < 0$ as expected.

Tangential momenta (counter-clockwise):

$$\vec{p}_1(0) = mv(0, 1), \quad \vec{p}_2(0) = mv(-1, 0), \quad \vec{p}_3(0) = mv(0, -1), \quad \vec{p}_4(0) = mv(1, 0)$$

Kinetic energy: $T_0 = 4 \cdot mv^2/2 = 2mv^2$. The speed at energy fraction $\eta$ is

$$v = \sqrt{\frac{\eta|U_0|}{2m}} = \sqrt{\frac{\eta Q^2 (2\sqrt{2}-1)}{2mR}}$$

### N = 6: Alternating Hexagon

Place six particles at the vertices of a regular hexagon of circumradius $R$:

$$\vec{r}_i(0) = R\!\left(\cos\!\tfrac{\pi(i-1)}{3},\;\sin\!\tfrac{\pi(i-1)}{3}\right), \qquad i = 1,\ldots,6$$

with charges $q_i = (-1)^{i+1} Q$. The vertex-to-vertex distances by ring are $d_1 = R$, $d_2 = R\sqrt{3}$, $d_3 = 2R$.

Pairs by ring:

| Ring $k$ | Count | Distance | Type | Subtotal |
| --- | --- | --- | --- | --- |
| $k=1$ (adjacent) | 6 | $R$ | unlike | $-6Q^2/R$ |
| $k=2$ (skip-one) | 6 | $R\sqrt{3}$ | like | $+6Q^2/(R\sqrt{3})$ |
| $k=3$ (opposite) | 3 | $2R$ | unlike | $-3Q^2/(2R)$ |

Initial potential energy:

$$U_0 = \frac{Q^2}{R}\!\left(-6 + \frac{6}{\sqrt{3}} - \frac{3}{2}\right) = \frac{Q^2}{R}\!\left(-\frac{15}{2} + 2\sqrt{3}\right)$$

Numerically, $2\sqrt{3} \approx 3.464$, so $U_0 \approx -4.036\,Q^2/R < 0$.

Kinetic energy: $T_0 = 6 \cdot mv^2/2 = 3mv^2$. Speed at energy fraction $\eta$:

$$v = \sqrt{\frac{\eta|U_0|}{3m}} = \sqrt{\frac{\eta Q^2\!\left(\tfrac{15}{2} - 2\sqrt{3}\right)}{3mR}}$$

### General Even N

The general formula

$$U_0 = \frac{NQ^2}{2R}\left[\sum_{k=1}^{N/2-1} \frac{(-1)^k}{\sin(k\pi/N)} + \frac{(-1)^{N/2}}{2}\right]$$

holds for any even $N$ with alternating $\pm Q$ charges on a regular $N$-gon of circumradius $R$. The speed for all particles at energy fraction $\eta$ is

$$v = \sqrt{\frac{2\eta|U_0|}{Nm}}$$

Note that as $N \to \infty$ the alternating-charge ring approaches a neutral charged ring, and $U_0$ converges to a finite limit; this is consistent with the fact that all rings with $N \geq 4$ have $U_0 < 0$.

## Method C — Three-Dimensional Symmetric Configurations

### The Rigid Rotation Construction

For a three-dimensional configuration of $n$ particles at positions $\{\vec{r}_i(0)\}$ with COM at the origin, choose a rotation axis $\hat{\omega}$ and angular speed $\omega > 0$, and set

$$\vec{p}_i(0) = m_i \left(\omega\hat{\omega} \times \vec{r}_i(0)\right)$$

Properties of this construction:

- $\sum_i \vec{p}_i = M \omega\hat{\omega} \times \vec{0} = \vec{0}$ (zero total momentum) ✓
- $\dot{r}_{ij}(0) = 0$ for all pairs (rigid rotation preserves all distances) ✓
- Angular momentum: $\vec{L} = I \omega \hat{\omega}$ where $I = \sum_i m_i r_{\perp,i}^2$ is the moment of inertia about the chosen axis and $r_{\perp,i} = |\vec{r}_i - (\hat{\omega} \cdot \vec{r}_i)\hat{\omega}|$ is the perpendicular distance of particle $i$ from the axis.
- Kinetic energy: $T_0 = I\omega^2/2$.

Given a target energy $E$ and computed $U_0$, the required angular speed is

$$\omega = \sqrt{\frac{2(E - U_0)}{I}}$$

This requires $E > U_0$. Different choices of axis $\hat{\omega}$ yield the same energy but different angular momenta; this freedom can be used to control $|\vec{L}|$.

### Regular Tetrahedron (N = 4)

Place four particles at the vertices of a regular tetrahedron of circumradius $R$. A convenient choice of coordinates with the $C_3$ symmetry axis along $\hat{z}$ is

$$\begin{aligned}
\vec{r}_1(0) &= \left(0,\; 0,\; R\right) \qquad \text{(apex)}\\
\vec{r}_2(0) &= \left(\tfrac{2\sqrt{2}}{3}R,\; 0,\; -\tfrac{R}{3}\right)\\
\vec{r}_3(0) &= \left(-\tfrac{\sqrt{2}}{3}R,\; \tfrac{\sqrt{6}}{3}R,\; -\tfrac{R}{3}\right)\\
\vec{r}_4(0) &= \left(-\tfrac{\sqrt{2}}{3}R,\; -\tfrac{\sqrt{6}}{3}R,\; -\tfrac{R}{3}\right)
\end{aligned}$$

All six edges have equal length $d = R \cdot 2\sqrt{6}/3 \approx 1.633\,R$.

For an alternating charge assignment with two positive and two negative particles — say $q_1 = q_3 = +Q$ and $q_2 = q_4 = -Q$ — there are 4 unlike pairs and 2 like pairs, all at the same edge length $d$:

$$U_0 = \frac{4(-Q^2)}{d} + \frac{2(+Q^2)}{d} = -\frac{2Q^2}{d} = -\frac{3Q^2}{\sqrt{6}\,R}$$

For rigid rotation about the $C_3$ axis ($\hat{\omega} = \hat{z}$), the apex particle lies on the axis ($r_{\perp,1} = 0$) and the three base particles each have perpendicular distance $r_\perp = 2\sqrt{2}\,R/3$:

$$I_{C_3} = m \cdot 0^2 + 3m \cdot \left(\frac{2\sqrt{2}}{3}R\right)^2 = \frac{8mR^2}{3}$$

Angular speed from target energy $E$:

$$\omega = \sqrt{\frac{2(E - U_0)}{I_{C_3}}} = \sqrt{\frac{3(E - U_0)}{4mR^2}}$$

Momenta: particle 1 on the axis has $\vec{p}_1(0) = \vec{0}$; the three base particles have momenta tangential to their circular paths around $\hat{z}$.

### Octahedron (N = 6)

Place six particles at the vertices of a regular octahedron of circumradius $R$:

$$\vec{r}_1(0) = (R,0,0),\; \vec{r}_2(0) = (0,R,0),\; \vec{r}_3(0) = (0,0,R)$$
$$\vec{r}_4(0) = (-R,0,0),\; \vec{r}_5(0) = (0,-R,0),\; \vec{r}_6(0) = (0,0,-R)$$

Assign charges so that opposite vertices carry unlike charges and same-axis positive and negative ends alternate: $q_1 = q_2 = q_3 = +Q$ and $q_4 = q_5 = q_6 = -Q$. Then each positive particle is directly opposite a negative particle.

Pair classification. The 15 pairs split as follows. Same-axis pairs (opposite ends):

| Pair | Charges | Distance | Type |
| --- | --- | --- | --- |
| $(1,4)$, $(2,5)$, $(3,6)$ | $(+Q,-Q)$ | $2R$ | unlike |

Cross-axis pairs (perpendicular axes):

| Charge combination | Count | Distance | Type |
| --- | --- | --- | --- |
| $(+Q,+Q)$ | 3 | $R\sqrt{2}$ | like |
| $(-Q,-Q)$ | 3 | $R\sqrt{2}$ | like |
| $(+Q,-Q)$ | 6 | $R\sqrt{2}$ | unlike |

The like and unlike cross-axis contributions cancel:

$$U_0^{\text{cross}} = 6 \cdot \frac{(-Q^2)}{R\sqrt{2}} + 6 \cdot \frac{(+Q^2)}{R\sqrt{2}} = 0$$

The initial potential energy is determined entirely by the three opposite-end unlike pairs:

$$U_0 = -\frac{3Q^2}{2R}$$

For rigid rotation about the $C_4$ axis ($\hat{\omega} = \hat{z}$), the two axial particles (vertices 3 and 6) lie on the axis and the four equatorial particles orbit at perpendicular distance $R$:

$$I_{C_4} = 0 + 0 + 4mR^2 = 4mR^2$$

$$\omega = \sqrt{\frac{E - U_0}{2mR^2}} = \sqrt{\frac{E + 3Q^2/(2R)}{2mR^2}}$$

### Axis Choice and Angular Momentum

For a fixed configuration and energy, different choices of rotation axis $\hat{\omega}$ give different moments of inertia $I$ and therefore different $\omega$, while preserving $T_0 = I\omega^2/2 = E - U_0$. The angular momentum magnitude $|\vec{L}| = I\omega$ then varies:

$$|\vec{L}| = I\omega = \sqrt{2 I (E - U_0)}$$

A large-$I$ axis (rotation in a plane that sweeps far from the axis) gives larger angular momentum at the same energy. A small-$I$ axis (rotation close to a symmetry axis) gives smaller angular momentum and larger $\omega$.

This degree of freedom is useful when simulating different dynamical regimes — e.g.,slowly tumbling vs.rapidly spinning — at the same total energy.

## Zöllner-Modified Initial Conditions

### Modified Pair Coupling

When the Zöllner electrogravitational extension is active with mismatch parameter $a > 0$, unlike-sign pairs receive an enhanced coupling

$$\kappa_{ij} = 1 + a \quad (q_i q_j < 0), \qquad \kappa_{ij} = 1 \quad (q_i q_j > 0)$$

The initial potential energy at $\dot{r}_{ij}(0) = 0$ becomes

$$U_0^{\mathrm{Z}} = \sum_{\substack{i < j \\ q_i q_j < 0}} \frac{(1+a) q_i q_j}{r_{ij}(0)} + \sum_{\substack{i < j \\ q_i q_j > 0}} \frac{q_i q_j}{r_{ij}(0)}$$

Since unlike pairs contribute negatively (attractive), increasing $a$ makes $U_0^{\mathrm{Z}}$ more negative than the standard $U_0$. Splitting:

$$U_0^{\mathrm{Z}} = U_0 + a \sum_{\substack{i < j \\ q_i q_j < 0}} \frac{q_i q_j}{r_{ij}(0)}$$

where $U_0$ is the standard ($\kappa = 1$) initial potential and the second term is always negative.

### Two-Body Circular Orbit with Zöllner

The circular orbit derivation in Section 4 carries through unchanged with $k$ replaced by $k^{\mathrm{Z}} = (1+a) q_1 q_2$ for an attractive pair:

$$v_{\mathrm{circ}}^{\mathrm{Z}} = \sqrt{\frac{(1+a)|q_1 q_2|}{\mu r_0}}$$

The circular orbit is faster than in standard Weber by the factor $\sqrt{1+a}$. Equivalently, for the same circular speed, the Zöllner radius is smaller by a factor of $1+a$.

### Energy Budget Comparison

For a fixed target energy $E < 0$ and the same geometric configuration:

- Standard Weber: $T_0 = E - U_0$
- Zöllner: $T_0^{\mathrm{Z}} = E - U_0^{\mathrm{Z}} > T_0$ (more kinetic energy required)

For a fixed energy fraction $\eta$ (same $\eta = T_0/|U_0|$):

- Standard Weber: $E = (\eta - 1)|U_0|$
- Zöllner: $E^{\mathrm{Z}} = (\eta - 1)|U_0^{\mathrm{Z}}|$, with $|E^{\mathrm{Z}}| > |E|$ for $\eta < 1$

The Zöllner system at the same $\eta$ is therefore more deeply bound (more negative energy) than the standard Weber system in the same geometry.

To compare Zöllner and standard Weber at the **same physical energy** $E$, use $U_0^{\mathrm{Z}}$ in the speed formula with $\kappa_{ij}$ values set appropriately:

$$v = \sqrt{\frac{2(E - U_0^{\mathrm{Z}})}{Nm}}$$

## Method D — Sub-Critical Like-Charge Pairs

Methods A–C construct initial conditions for attractive ($k = q_1 q_2 < 0$) or mixed-sign systems. For two like charges ($k = q_1 q_2 > 0$), Weber's velocity-dependent force creates a fundamentally different regime: below a critical separation $\rho$, the pair forms a bound oscillatory state with *positive* total energy. This section gives closed-form initial conditions for such systems. For the underlying physics see [CriticalRadiusAndLikeChargeAttraction.md](../research/investigations/CriticalRadiusAndLikeChargeAttraction.md).

### The Critical Radius

For a like-charge pair with $k = q_1 q_2 > 0$, define

$$\rho = \frac{k}{\mu c^2} = \frac{q_1 q_2}{\mu c^2}$$

where $\mu = m_1 m_2 / M$ is the reduced mass and $c$ is the speed of light. For like charges $\rho > 0$; for unlike charges $\rho < 0$ and the sub-critical regime does not exist.

The pair dynamics bifurcate at $r = \rho$:

- **Molecular state** ($r_0 < \rho$): bound radial oscillation between $r_0$ and $r = 0$.
- **Distant state** ($r_0 > \rho$): Coulomb-like repulsive scattering.

No continuous trajectory crosses $r = \rho$; it is an absolute dynamical barrier.

### Radial Oscillation from Turning Point

The simplest sub-critical initial condition places two like charges at the outer turning point of their oscillation.

**Positions** (COM frame, along $x$-axis):

$$\vec{r}_1(0) = -\frac{m_2}{M} r_0 \,\hat{x}, \qquad \vec{r}_2(0) = +\frac{m_1}{M} r_0 \,\hat{x}$$

with $r_0 < \rho$.

**Momenta**: both zero.

$$\vec{p}_1(0) = \vec{0}, \qquad \vec{p}_2(0) = \vec{0}$$

At the turning point $\dot{r}_{12} = 0$, the Weber potential reduces to Coulomb:

$$E = H(0) = \frac{k}{r_0} > 0$$

The energy is positive, yet the system is permanently bound. The pair oscillates between $r_0$ and $r = 0$, never escaping past $\rho$.

### Radial Oscillation from Arbitrary Phase

To start at separation $r < r_0$ with nonzero radial velocity (mid-oscillation), use the integrated energy equation. From Hamiltonian conservation with $E = k/r_0$:

$$\dot{r}^2 = \frac{2\rho c^2 (r_0 - r)}{r_0 (\rho - r)}$$

Both $r_0 - r > 0$ and $\rho - r > 0$ (since $r < r_0 < \rho$), so $\dot{r}^2 > 0$. As $r \to 0$, $\dot{r}^2 \to 2c^2$.

Choose the sign of $\dot{r}$: negative for approaching, positive for receding. The momenta are

$$\vec{p}_1(0) = +\mu\dot{r}\,\hat{r}, \qquad \vec{p}_2(0) = -\mu\dot{r}\,\hat{r}$$

where $\hat{r} = (\vec{r}_1 - \vec{r}_2)/r$. These satisfy $\vec{p}_1 + \vec{p}_2 = \vec{0}$.

**Verification.** Since $\dot{r} \neq 0$, the full velocity-dependent potential enters:

$$H(0) = \frac{\mu\dot{r}^2}{2} + \frac{k}{r}\left(1 - \frac{\dot{r}^2}{2c^2}\right) = \frac{k}{r_0}$$

Substituting $\dot{r}^2$ from the integrated equation confirms this identity.

### Oscillation with Transverse Velocity

To add angular momentum, assign tangential momenta at the turning point. Positions as in the radial case with $r_0 < \rho$. Momenta:

$$\vec{p}_1(0) = +p_\perp \,\hat{y}, \qquad \vec{p}_2(0) = -p_\perp \,\hat{y}$$

where $p_\perp = \mu \alpha_0$ and $\alpha_0$ is the relative transverse speed. The initial energy is

$$E = \frac{p_\perp^2}{2\mu} + \frac{k}{r_0} = \frac{\mu\alpha_0^2}{2} + \frac{k}{r_0}$$

and the angular momentum magnitude is $L = p_\perp r_0 = \mu \alpha_0 r_0$.

Unlike the attractive case, no stable circular orbit exists for like charges below $\rho$. The sub-critical attraction always exceeds the centripetal requirement, so the pair oscillates radially with an inward-spiralling tendency regardless of the transverse velocity.

### Prescribed-Energy Construction

**Inverse problem.** Given target energy $E > 0$, find the turning-point amplitude.

For purely radial oscillation ($\alpha_0 = 0$):

$$r_0 = \frac{k}{E}$$

Sub-critical feasibility requires $r_0 < \rho$, i.e.,$E > k/\rho = \mu c^2$. The minimum sub-critical energy (maximum amplitude $r_0 = \rho$) is

$$E_{\min} = \frac{k}{\rho} = \mu c^2$$

For $E < \mu c^2$, the turning point lies in the distant state ($r_0 > \rho$) and the pair scatters rather than oscillates.

With transverse velocity:

$$r_0 = \frac{k}{E - \mu\alpha_0^2/2}$$

This requires $E > \mu\alpha_0^2/2$.

| $E$ | $r_0$ vs $\rho$ | Regime |
| --- | --- | --- |
| $E > \mu c^2$ | $r_0 < \rho$ | molecular (bound oscillation) |
| $E = \mu c^2$ | $r_0 = \rho$ | critical (barrier) |
| $0 < E < \mu c^2$ | $r_0 > \rho$ | distant (scattering) |

### Period and Timescale

The oscillation period for small amplitude ($r_0 \ll \rho$) with $\alpha_0 = 0$ is

$$T = \frac{2\sqrt{2}\, r_0}{c}$$

Period is proportional to amplitude $r_0$. The coefficient decreases from $2\sqrt{2}/c$ at small amplitude to $\sqrt{2}/c$ at maximum amplitude $r_0 = \rho$.

Simulation parameters should be chosen accordingly:

- $\mathrm{tspan}$: several multiples of $T$ to capture multiple oscillation cycles.
- $\mathrm{dt}$: small enough to resolve the oscillation; at least 20 steps per period.

**Regularization.** As $r \to 0$ the radial velocity approaches $\sqrt{2}\,c$ and the Coulomb singularity $k/r$ diverges. Regularization (see [Regularization.md](Regularization.md)) is essential for accurate integration. Use `:adaptive_cartesian` for 3D or `:lifted_pair` for 2D.

## Distant-State Like-Charge Scattering

For completeness, this section treats the super-critical regime $r_0 > \rho$.

Two like charges approaching from infinity with kinetic energy $E > 0$ reach closest approach $r_0 = k/E > \rho$ and recede. The initial conditions at the turning point have the same form as the molecular case — COM positions with zero momenta — but with $r_0 > \rho$. After the turning point, the pair separates to infinity rather than oscillating.

The distant state requires $E < \mu c^2$ (equivalently $r_0 > \rho$). For $E = \mu c^2$ the turning point coincides with the barrier; for $E > \mu c^2$ the pair is in the molecular state.

To start the scattering pair at separation $r > r_0$ with inward radial velocity:

$$\dot{r}^2 = \frac{2\rho c^2 (r - r_0)}{r_0 (r - \rho)}$$

Both $r - r_0 > 0$ and $r - \rho > 0$ (since $r > r_0 > \rho$), so $\dot{r}^2 > 0$. The maximum radial speed is bounded: $\dot{r}^2 < 2(\rho/r_0) c^2 < 2c^2$.

## Mixed Systems — Sub-Critical Nucleus with Orbiter

Weber's "planetary model of the atom" combines a sub-critical like-charge pair with orbiting unlike charges (see [CriticalRadiusAndLikeChargeAttraction.md](../research/investigations/CriticalRadiusAndLikeChargeAttraction.md)). The simplest realisation is a three-body system.

### Three-Body Configuration

Consider two like charges $q_1 = q_2 = +q$ with mass $m$ forming a sub-critical nucleus at separation $r_{\text{nuc}} < \rho$, plus one unlike charge $q_3 = -q'$ with mass $m_3$ orbiting at distance $R \gg r_{\text{nuc}}$.

**Nucleus.** Use the turning-point construction for particles 1 and 2. Place them symmetrically about the nucleus centre with $r_{\text{nuc}} < \rho$ and zero radial momenta.

**Orbiter.** At distance $R \gg r_{\text{nuc}}$, the orbiter sees an effective charge $Q_{\text{eff}} \approx 2q$ at the nucleus centre. The orbital speed for an approximate circular orbit is

$$v_{\text{orb}} \approx \sqrt{\frac{2q \cdot q'}{\mu_{\text{orb}}\, R}}$$

where $\mu_{\text{orb}} = (2m) m_3 / (2m + m_3)$ is the reduced mass of the orbiter relative to the nucleus.

**COM adjustment.** All three positions must satisfy $\sum_i m_i \vec{r}_i(0) = \vec{0}$. After placing the nucleus pair and the orbiter, shift all positions by $-\vec{R}_{\text{COM}}$ and adjust momenta to ensure $\sum_i \vec{p}_i(0) = \vec{0}$.

Scale the orbital speed by a fraction $\eta_{\text{orb}} \in (0, 1)$ to produce sub-circular (elliptical) orbits:

$$v = \eta_{\text{orb}} \cdot v_{\text{orb}}$$

Bound orbits require $\eta_{\text{orb}} \lesssim 0.9$; at $\eta_{\text{orb}} \geq 1$ the orbiter escapes. The exact bound/unbound threshold depends on the nucleus size and Zöllner coupling.

**Collision bounce.** The nucleus pair oscillation passes through $r = 0$ (in the $\ell = 0$ regularisable case). Use a collision bounce radius $r_{\text{bounce}} > 0$ to reflect the relative coordinate at small separation (see [../exploratory/CollisionBounceRegularization.md](../exploratory/CollisionBounceRegularization.md)).

**Zöllner enhancement.** With mismatch parameter $a > 0$, the unlike-pair coupling strengthens by factor $(1+a)$. This tightens the orbiter's orbit and can circularise an otherwise eccentric trajectory. For $\eta_{\text{orb}} = 0.8$ and $a = 0.5$, the orbit becomes nearly circular. The nucleus (like-charge pair) is unaffected by $a$ since $\kappa_{ij} = 1$ for like signs.

This construction is approximate: the nucleus oscillation and orbital motion couple dynamically. The scale separation $R \gg r_{\text{nuc}}$ controls the quality of the approximation. Validated configurations and parameter sweeps are documented in [../exploratory/ThreeBodyBoundStates.md](../exploratory/ThreeBodyBoundStates.md).

## Reference Tables

### Bound Orbit Criteria

| Configuration | Bound iff | Critical energy |
| --- | --- | --- |
| Any attractive configuration ($k < 0$) | $E < 0$ | $E = 0$ (escape) |
| Two-body attractive, fixed $r_0$ | $\eta_v < \sqrt{2}$ | $E = 0$, $\eta_v = \sqrt{2}$ |
| Two-body attractive, fixed $r_0$, fix $E$ | $r_0 < |k|/|E|$ | $r_0 = |k|/|E|$ (zero velocity) |
| $N$-body polygon, $\eta$ | $\eta < 1$ | $\eta = 1$, $E = 0$ |
| $N$-body polygon, fix $E < 0$ | $T_0 < |U_0|$ | $T_0 = |U_0|$, $E = 0$ |
| Like-charge pair, sub-critical | $r_0 < \rho$ | $E = \mu c^2$ (barrier) |
| Like-charge pair, distant | $r_0 > \rho$ | always unbound (scattering) |

### Common Configurations

| System | $U_0$ | $v$ or $\omega$ |
| --- | --- | --- |
| 2-body attractive, sep.$r_0$ | $k/r_0$ | $\eta_v\sqrt{|k|/(\mu r_0)}$ |
| 4-body square, radius $R$ | $(Q^2/R)(1-2\sqrt{2})$ | $\sqrt{\eta|U_0|/(2m)}$ |
| 6-body hexagon, radius $R$ | $(Q^2/R)(-\tfrac{15}{2}+2\sqrt{3})$ | $\sqrt{\eta|U_0|/(3m)}$ |
| Tetrahedron, radius $R$ | $-3Q^2/(\sqrt{6}\,R)$ | $\omega = \sqrt{3(E-U_0)/(4mR^2)}$ |
| Octahedron $C_4$, radius $R$ | $-3Q^2/(2R)$ | $\omega = \sqrt{(E-U_0)/(2mR^2)}$ |
| Like-charge pair, turning pt $r_0 < \rho$ | $k/r_0$ | $\vec{p}_i = \vec{0}$ |
| Like-charge pair, transverse $\alpha_0$ | $k/r_0$ | $p_\perp = \mu\alpha_0$ tangential |
| Planetary atom $(+\!+\!-)$, $R \gg r_{\text{nuc}}$ | $\approx 2q^2/r_{\text{nuc}} - 2q^2/R$ | $\eta_{\text{orb}}\sqrt{2q^2/(\mu_{\text{orb}} R)}$ tangential |

## Step-by-Step Recipe

Given: masses $m_i$, charges $q_i$, speed of light $c$, coupling factors $\kappa_{ij}$, and a target energy $E$.

1. **Choose geometry.** Select a symmetry class (two-body, $N$-gon, polyhedron) and geometry parameters ($r_0$ or $R$). For like-charge pairs in the sub-critical regime ($k > 0$, $r_0 < \rho$), see Method D. Note that sub-critical pairs have $E > 0$ (bound with positive energy), unlike the attractive case where bound means $E < 0$.

2. **Place particles.** Assign positions satisfying the COM condition $\sum_i m_i \vec{r}_i(0) = \vec{0}$.

3. **Compute $U_0$.** Evaluate the Coulomb sum (with $\kappa_{ij}$):

$$U_0 = \sum_{i < j} \frac{\kappa_{ij} q_i q_j}{r_{ij}(0)}$$

4. **Check feasibility.** Require $E > U_0$ so that $T_0 = E - U_0 > 0$. If $E \leq U_0$, increase $r_0$ (or $R$) or choose a less negative target energy.

5. **Determine momenta.** For 2D configurations assign tangential momenta; for 3D use the rigid rotation construction with a chosen axis $\hat{\omega}$ and

$$\omega = \sqrt{2(E - U_0)/I}, \qquad \vec{p}_i(0) = m_i \omega\hat{\omega} \times \vec{r}_i(0)$$

6. **Verify.** Compute $H(0) = \sum_i |\vec{p}_i|^2/(2m_i) + U_0$ and confirm it equals $E$. Check $\sum_i \vec{p}_i = \vec{0}$.
