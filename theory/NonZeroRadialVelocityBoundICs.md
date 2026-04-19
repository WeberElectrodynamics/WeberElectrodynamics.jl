# Non-Zero-Radial-Velocity Initial Conditions for Bound Weber Orbits

## Scope

This document gives closed-form initial conditions for $n$-body Weber Hamiltonian systems in which **at least one pair has non-zero relative radial velocity at $t = 0$** (i.e. $\dot{r}_{ij}(0) \neq 0$) while the total system is **attractively bound** (i.e. $E = T_0 + U_0 < 0$ with $|U_0| > T_0$).

It complements [`InitialConditions.md`](InitialConditions.md), whose Methods A (two-body), B (planar $N$-gons), and C (polyhedra) all enforce the Zero-Radial-Velocity Principle $\dot{r}_{ij}(0) = 0$ to simplify the Hamiltonian constraint, and whose Method D admits $\dot{r} \neq 0$ only for like-charge **sub-critical** pairs with $E > 0$ (repulsive barrier).

This document fills the remaining regime — $E < 0$, $\dot{r}_{ij}(0) \neq 0$ — which is relevant to:

- Mid-flight snapshots of attractive two-body orbits (any eccentric anomaly, not just apsides).
- Infalling / expanding binaries inside multi-body configurations.
- Breathing regular $N$-gons.

All conventions (absolute units, Zöllner couplings $\kappa_{ij}$, notation for relative coordinates) follow [`WeberElectrodynamics.md`](WeberElectrodynamics.md) and [`InitialConditions.md`](InitialConditions.md) unchanged.

---

## 1. Canonical vs. kinetic momentum

The Weber Lagrangian depends on velocities through both $T$ and $\dot{r}_{ij}$, so the canonical momentum $\vec{p}_i \equiv \partial L/\partial \vec{v}_i$ differs from the kinetic momentum $\vec{p}_i^{\,\mathrm{kin}} \equiv m_i \vec{v}_i$ whenever any pair has $\dot{r}_{ij} \neq 0$. All Methods A/B/C/D-turning-point constructions impose $\dot{r}_{ij}(0) = 0$, so for those constructions $\vec{p}_i = m_i \vec{v}_i$ and no conversion is needed. For the constructions in this document the conversion is non-trivial — it is derived once in Section 6 and used throughout.

The short statement: the user supplies physical velocities $\{\vec{v}_i\}$, and the integrator initial state is $\{\vec{r}_i, \vec{p}_i\}$ with $\vec{p}_i$ from the forward map $(F)$ in Section 6. Setting $\vec{p}_i = m_i \vec{v}_i$ is **wrong** in this regime and silently corrupts the energy invariant.

---

## 2. The two-body radial-phase equation

### 2.1 Derivation

For a two-body Weber system with reduced mass $\mu = m_1 m_2/M$ (total mass $M = m_1 + m_2$), coupling $k = \kappa_{12} q_1 q_2$, and conserved angular momentum $L = \mu r^2 \dot{\theta}$, the Hamiltonian (in the COM frame) is

$$H \;=\; T + U \;=\; \tfrac{1}{2}\mu\!\left(\dot{r}^2 + \frac{L^2}{\mu^2 r^2}\right) \;+\; \frac{k}{r}\!\left(1 - \frac{\dot{r}^2}{2c^2}\right) \;=\; E.$$

Collecting terms quadratic in $\dot{r}$ and solving:

$$\boxed{\;\dot{r}^2 \;=\; \frac{2\bigl(E \;-\; k/r \;-\; L^2/(2\mu r^2)\bigr)}{\mu \;-\; k/(c^2 r)}\;} \qquad (\star)$$

This is the **Weber radial energy equation**. The form $(\star)$ matches the equation given by Clemente & Assis (*Int. J. Theor. Phys.* **30**, 537 (1991), DOI [10.1007/BF00672899](https://doi.org/10.1007/BF00672899)) and is equivalent to eq. (2.3) of Frauenfelder & Weber (*Anal. Math. Phys.* **14**:31 (2024), DOI [10.1007/s13324-024-00891-5](https://doi.org/10.1007/s13324-024-00891-5)). Sanity check: in the Kepler limit $c\to\infty$ the denominator $\to\mu$ and $(\star)$ reduces to the textbook $\dot{r}^2 = (2/\mu)[E - k/r - L^2/(2\mu r^2)]$. The specialization $L=0$ appears in §2.5.

### 2.2 Denominator sign

For attractive $k < 0$, the denominator $\mu - k/(c^2 r) > \mu$ for every $r > 0$; $(\star)$ is globally regular on the bound radial interval $[r_p, r_a]$.

For repulsive $k > 0$ the denominator vanishes at $r_\star = k/(\mu c^2) = \rho$, where $\rho = k/(\mu c^2)$ is the Method-D critical radius. **The pole of $(\star)$ coincides exactly with $\rho$.** This matches the physical intuition that like-charge sub-critical oscillation (Method D, $r < \rho$) lies in the regime where the energy-to-radial-velocity inversion is not analytic, and never arises in the attractive bound case considered here.

### 2.3 Turning points

Setting $\dot{r} = 0$ in $(\star)$ reduces the numerator to the Kepler quadratic

$$\frac{L^2}{2\mu}\,u^2 \;+\; k\,u \;-\; E \;=\; 0, \qquad u = 1/r,$$

with roots $u_{p,a}$ giving periapsis $r_p = 1/u_p$ and apoapsis $r_a = 1/u_a$. The Weber correction term vanishes at turning points, so the Weber turning-point radii are identical to the Coulomb/Kepler ones. A bound two-body Weber orbit exists iff

$$E < 0 \qquad\text{and}\qquad L^2 < \frac{\mu k^2}{2|E|}.$$

### 2.4 Recipe from $(E, L, r, s)$

Given target energy $E < 0$, angular momentum $L \geq 0$ with $L^2 < \mu k^2/(2|E|)$, a radius $r \in [r_p, r_a]$, and a sign $s \in \{-1, +1\}$ (inbound / outbound):

1. **Radial velocity.**
$$\dot{r} \;=\; s\,\sqrt{\frac{2\bigl(E \;-\; k/r \;-\; L^2/(2\mu r^2)\bigr)}{\mu \;-\; k/(c^2 r)}}.$$

2. **Tangential velocity.** $v_\perp = L/(\mu r)$.

3. **Relative physical velocity.** $\vec{v}_{\mathrm{rel}} = \dot{r}\,\hat{r} + v_\perp\,\hat{\theta}$. COM-frame particle velocities:
$$\vec{v}_1 = \frac{m_2}{M}\,\vec{v}_{\mathrm{rel}}, \qquad \vec{v}_2 = -\frac{m_1}{M}\,\vec{v}_{\mathrm{rel}}.$$

4. **Canonical momenta via $(F)$.** Placing the pair on the $x$-axis so that $\hat{r} = \hat{x}$, $\hat{\theta} = \hat{y}$:
$$\vec{p}_1 \;=\; \left(\mu\dot{r}\!\left(1 \;-\; \frac{k}{\mu c^2 r}\right),\; \mu v_\perp,\; 0\right), \qquad \vec{p}_2 \;=\; -\vec{p}_1.$$

The factor $1 - k/(\mu c^2 r) = 1 - \rho/r$ is the Weber correction on the radial component; for attractive $k < 0$ it is strictly greater than $1$.

### 2.5 Limits

- $r \to r_p$ or $r \to r_a$ ⇒ $\dot{r} \to 0$ and the recipe collapses continuously to Method A's turning-point IC.
- $L = 0$ (radial fall): $v_\perp = 0$, and $(\star)$ becomes $\dot{r}^2 = 2(E - k/r)/(\mu - k/(c^2 r))$, which matches Newtonian radial fall up to an $\mathcal{O}(1/c^2)$ Weber correction.

### 2.6 Phase parameterisation

An IC may be labelled by the eccentric anomaly along the Keplerian reference ellipse with semi-major axis $a = -k/(2E)$ and eccentricity $e = \sqrt{1 + 2EL^2/(\mu k^2)}$:

$$r(E_{\mathrm{anom}}) \;=\; a\,(1 - e\cos E_{\mathrm{anom}}), \qquad E_{\mathrm{anom}} \in [0, 2\pi),$$

with $s = +1$ on $(0, \pi)$ (outbound from periapsis) and $s = -1$ on $(\pi, 2\pi)$ (inbound). The reference ellipse is bookkeeping only — Weber orbits precess; the subsequent integrated trajectory is not Keplerian.

---

## 3. Hot binary + cold outer orbiter (3-body, asymptotic)

### 3.1 Scope

Particles 1 and 2 form a "hot binary" with $\dot{r}_{12}(0) \neq 0$; particle 3 is a "cold orbiter" at large distance. Cross-pair radial rates $\dot{r}_{13}, \dot{r}_{23}$ are not exactly zero — they are $\mathcal{O}(r_0/R)$ where $r_0$ is the binary separation and $R$ is the orbiter's distance from the binary COM. For $R/r_0 \gtrsim 10$ the cross-pair Weber corrections are a $\sim 1\%$ perturbation and the construction is practically exact. Strictly closed-form configurations in which every pair has $\dot{r}_{ij} \neq 0$ are covered in Section 4.

### 3.2 Geometry

Orient the binary internal axis along $\hat{z}$ and place the orbiter in the $xy$-plane:

$$\vec{r}_1 \;=\; \vec{R}_B \;-\; \frac{m_2}{m_1+m_2}\,r_0\,\hat{z}, \qquad \vec{r}_2 \;=\; \vec{R}_B \;+\; \frac{m_1}{m_1+m_2}\,r_0\,\hat{z}, \qquad \vec{r}_3 \;=\; \vec{R}_B + R\,\hat{x},$$

with the system COM at the origin fixing the binary COM:

$$\vec{R}_B \;=\; -\frac{m_3}{M}\,R\,\hat{x}, \qquad M = m_1 + m_2 + m_3.$$

### 3.3 Binary internal motion

Pick target binary energy $E_B < 0$, angular momentum $L_B$ (can be zero), radius $r_0 \in [r_p^B, r_a^B]$, sign $s$. Apply the Section 2 recipe with the internal axis along $\hat{z}$ to obtain $\dot{r}_{12}$, $v_{\perp, B}$, and the internal relative velocity $\vec{v}_{\mathrm{rel}}^B$. Split into binary-internal particle velocities:

$$\vec{v}_1^{\,B} \;=\; \frac{m_2}{m_1+m_2}\,\vec{v}_{\mathrm{rel}}^{\,B}, \qquad \vec{v}_2^{\,B} \;=\; -\frac{m_1}{m_1+m_2}\,\vec{v}_{\mathrm{rel}}^{\,B}.$$

### 3.4 Orbiter motion

The orbiter sees an effective attractor of charge $q_1 + q_2$ at the binary COM. Apply Method A of [`InitialConditions.md`](InitialConditions.md) (circular or elliptical) with reduced mass $\mu_{\mathrm{out}} = (m_1 + m_2)m_3/M$ and coupling $k_{\mathrm{out}} = \kappa_{\mathrm{out}}(q_1 + q_2)q_3$, obtaining the orbital speed $V_{\mathrm{rel}}$ perpendicular to $\vec{R}$. Place $\vec{V}_B$ and $\vec{v}_3$ along $\hat{y}$ with $\sum_i m_i \vec{v}_i = \vec{0}$.

### 3.5 Cross-pair radial rates

A direct computation of $\dot{r}_{13}(0) = (\vec{r}_1 - \vec{r}_3)\cdot(\vec{v}_1 - \vec{v}_3)/r_{13}$ yields a leading term $\mathcal{O}(r_0/R) \cdot v_{\perp, B}$. When the binary is radial ($L_B = 0$, $v_{\perp, B} = 0$), the leading term vanishes identically and the residual is $\mathcal{O}((r_0/R)^2) \cdot \dot{r}_{12}$.

### 3.6 Energy budget

To leading order in $r_0/R$:

$$U_0 \;=\; \frac{k_{12}}{r_0}\!\left(1 \;-\; \frac{\dot{r}_{12}^2}{2c^2}\right) \;+\; \frac{k_{13}}{R} \;+\; \frac{k_{23}}{R},$$

$$T_0 \;=\; \tfrac{1}{2}\mu_{12}\dot{r}_{12}^2 \;+\; \frac{L_B^2}{2\mu_{12}r_0^2} \;+\; \tfrac{1}{2}\mu_{\mathrm{out}}V_{\mathrm{rel}}^2.$$

The system is bound iff the binary and outer orbits are individually bound ($E_B < 0$ and $E_{\mathrm{out}} < 0$) — a standard hierarchical sufficient condition at large scale separation.

### 3.7 Canonical momenta

Only the $(1,2)$ pair contributes a Weber correction to leading order:

$$\vec{p}_1 \;=\; m_1 \vec{v}_1 \;-\; \frac{k_{12}\dot{r}_{12}}{c^2 r_{12}^2}(\vec{r}_1 - \vec{r}_2), \qquad \vec{p}_2 \;=\; m_2 \vec{v}_2 \;+\; \frac{k_{12}\dot{r}_{12}}{c^2 r_{12}^2}(\vec{r}_1 - \vec{r}_2), \qquad \vec{p}_3 \;=\; m_3 \vec{v}_3.$$

### 3.8 Worked sketch

Set $m_1 = m_2 = m_3 = 1$, $q_1 = +Q$, $q_2 = -Q$, $q_3 = +Q$, $r_0 = 1$, $R = 20$. Choose $L_B = 0$ and target binary energy $E_B = -0.3\,|k_{12}|/r_0$ at radius $r_0$ with $s = -1$ (infalling). Section 2 gives $\dot{r}_{12}(0)$. The binary's net charge is zero, so to give the orbiter a net attraction one must use a charge asymmetry (e.g. give particle 3 a different magnitude) or enable Zöllner enhancement so that $\kappa_{13} \neq \kappa_{23}$ produces a net attractive $k_{13} + k_{23}$.

---

## 4. Breathing regular $N$-gon (exact, $N$-body)

### 4.1 Geometry

$N$ equal-mass particles at angular positions $\theta_i = 2\pi(i-1)/N$ on a circle of radius $R$ in the $xy$-plane, with alternating charges $q_i = (-1)^{i+1}Q$ (so $N$ must be even):

$$\vec{r}_i(0) \;=\; R\,\hat{r}_i, \qquad \hat{r}_i = (\cos\theta_i, \sin\theta_i, 0).$$

### 4.2 Velocities

Superimpose a common radial breathing rate $\dot{R} \neq 0$ and a rigid rotation $\omega$:

$$\vec{v}_i(0) \;=\; \dot{R}\,\hat{r}_i \;+\; R\omega\,\hat{\theta}_i, \qquad \hat{\theta}_i = (-\sin\theta_i, \cos\theta_i, 0).$$

Both terms individually satisfy $\sum_i m\vec{v}_i = \vec{0}$ by the centred-polygon identity $\sum_i \hat{r}_i = \vec{0}$.

### 4.3 Key lemma — common radial-rate ratio

Pair separation is $r_{ij} = 2R\sin(k\pi/N) \equiv R\,d_k$ with $k = \min(|i-j|, N-|i-j|)$ and $d_k = 2\sin(k\pi/N)$. Since the angular separation $\theta_i - \theta_j$ is time-independent under common rotation, $r_{ij}(t) = R(t)\,d_k$ and therefore

$$\boxed{\;\dot{r}_{ij}(0) \;=\; d_k\,\dot{R}\;}$$

is non-zero for every pair whenever $\dot{R} \neq 0$. All pairs in ring $k$ share the same $\dot{r}_{ij}$, and **every pair shares the same ratio** $\dot{r}_{ij}/r_{ij} = \dot{R}/R$. This common-ratio property is what reduces the canonical-momentum formula $(F)$ to a single scalar.

### 4.4 Initial potential energy

Split Coulomb and Weber contributions. Using the ring counting of [`InitialConditions.md`](InitialConditions.md) (§"Computing $U_0$ for a Regular Alternating-Charge Polygon") — with $n_k = N$ for $k < N/2$ and $n_{N/2} = N/2$, signs $s_k = (-1)^k$:

$$U_0 \;=\; U_0^{\mathrm{Coul}} \;-\; \frac{\dot{R}^2\,Q^2}{2c^2 R}\,\Sigma_N, \qquad \Sigma_N \;=\; \sum_{k=1}^{N/2} n_k\,s_k\,d_k.$$

Explicit values:

- **$N = 4$** (square): $d_1 = \sqrt{2}$, $d_2 = 2$.  
  $\Sigma_4 \;=\; 4(-1)\sqrt{2} \;+\; 2(+1)(2) \;=\; 4 - 4\sqrt{2} \;\approx\; -1.657$.

- **$N = 6$** (hexagon): $d_1 = 1$, $d_2 = \sqrt{3}$, $d_3 = 2$.  
  $\Sigma_6 \;=\; 6(-1)(1) + 6(+1)\sqrt{3} + 3(-1)(2) \;=\; -12 + 6\sqrt{3} \;\approx\; -1.608$.

In both cases $\Sigma_N < 0$, so the Weber correction $-\dot{R}^2 Q^2 \Sigma_N/(2c^2 R)$ is **positive**: breathing motion reduces $|U_0|$ relative to the static polygon. The construction therefore requires $|\dot{R}|$ small enough that $|U_0^{\mathrm{Coul}}|$ still dominates the kinetic energy.

### 4.5 Kinetic energy

$$T_0 \;=\; \tfrac{1}{2}\,N\,m\,\bigl(\dot{R}^2 + R^2\omega^2\bigr).$$

### 4.6 Bound criterion

$|U_0| > T_0$, i.e.

$$|U_0^{\mathrm{Coul}}| \;-\; \frac{\dot{R}^2 Q^2\,|\Sigma_N|}{2c^2 R} \;>\; \tfrac{1}{2}\,N\,m\,\bigl(\dot{R}^2 + R^2\omega^2\bigr).$$

In the $\dot{R} \to 0$ limit this recovers the Method-B criterion exactly.

### 4.7 Canonical momenta — closed form

Apply $(F)$ to particle $i$ with $\dot{r}_{ij} = d_k\dot{R}$, $r_{ij}^2 = R^2 d_k^2$, $\vec{r}_i - \vec{r}_j = R(\hat{r}_i - \hat{r}_j)$:

$$\vec{p}_i \;=\; m\vec{v}_i \;-\; \frac{\dot{R}}{c^2 R}\sum_{j\ne i}\frac{\kappa_{ij}q_i q_j}{d_k}\,(\hat{r}_i - \hat{r}_j).$$

Using $(\hat{r}_i - \hat{r}_j)\cdot \hat{r}_i = 1 - \cos(2\pi k/N)$ and the identity $(1 - \cos 2\theta)/(2\sin\theta) = \sin\theta$, so $(1 - \cos(2\pi k/N))/d_k = d_k/2$. The projection onto $\hat{r}_i$ collapses by symmetry: the transverse components from ring-$k$ partners $j = i \pm k$ cancel. Counting per-particle multiplicities $m_k = 2$ for $k < N/2$ (two partners $j = i\pm k$) and $m_{N/2} = 1$, together with $n_k = N\,m_k/2$ (each pair belongs to one ring but is counted at both endpoints), the per-particle sum evaluates to

$$\sum_{j\ne i}\frac{\kappa_{ij}q_i q_j}{d_k}\,(\hat{r}_i - \hat{r}_j)\cdot\hat{r}_i \;=\; \frac{Q^2}{2}\sum_{k=1}^{N/2} m_k\,s_k\,d_k \;=\; \frac{Q^2\,\Sigma_N}{N}.$$

Therefore:

$$\boxed{\;\vec{p}_i \;=\; m\,\vec{v}_i \;-\; \frac{\dot{R}\,Q^2}{c^2 R}\,\frac{\Sigma_N}{N}\,\hat{r}_i.\;}$$

Since $\sum_i \hat{r}_i = \vec{0}$, the canonical total automatically satisfies $\sum_i \vec{p}_i = \sum_i m\vec{v}_i = \vec{0}$.

### 4.8 Worked example ($N = 4$)

$m = Q = R = 1$, $\dot{R} = -0.1\,c$ (contracting), alternating charges. Method B gives $U_0^{\mathrm{Coul}} = (1 - 2\sqrt{2})Q^2/R \approx -1.828\,Q^2/R$. The Weber correction:

$$-\frac{\dot{R}^2 Q^2 \Sigma_4}{2c^2 R} \;=\; -\frac{(0.01 c^2)(4 - 4\sqrt{2})}{2c^2}\cdot\frac{Q^2}{R} \;\approx\; +0.00828\,Q^2/R.$$

So $U_0 \approx -1.820\,Q^2/R$. Bound iff $T_0 < 1.820$, i.e. $2m\dot{R}^2 + 2mR^2\omega^2 < 1.820$, i.e. $0.02 + 2\omega^2 < 1.820$, giving $\omega^2 < 0.9$. Canonical momentum correction:

$$-\frac{\dot{R}\,Q^2}{c^2 R}\,\frac{\Sigma_4}{4}\,\hat{r}_i \;\approx\; +0.0414\,c^{-1}\,\hat{r}_i$$

— small in absolute units.

### 4.9 Periodicity warning

Bound total energy ($E < 0$) does **not** imply a periodic trajectory. Periodic breathing-square orbits exist but are generally highly unstable. This document supplies **bound initial conditions only**; long-time dynamics requires numerical integration and is generically aperiodic.

---

## 5. Out of scope

The following are physically interesting but do not reduce to a clean closed-form recipe of the kind above:

- **Pulsating regular polyhedra beyond the tetrahedron.** Only the tetrahedron has a single inter-vertex distance. In the octahedron and larger polyhedra, equatorial and axial chords scale differently with the radius, so the common-ratio property of §4.3 fails. A two-scale breathing ansatz is possible but loses the single-scalar inversion.
- **Inspiraling binary + orbiter at comparable scales.** When $R/r_0 \sim \mathcal{O}(1)$, cross-pair radial rates are not small and the Section 3 asymptotic expansion breaks; cross-pair Weber corrections enter at the same order as the binary's.
- **General many-body mid-flight configurations.** For an arbitrary geometry and velocity field satisfying only $\sum_i m_i \vec{v}_i = \vec{0}$, there is no closed-form map from target energy to canonical momenta. Users should either integrate backward from a Method A/B/C IC or run numerical shooting with a target Hamiltonian.

---

## 6. Canonical momentum map

### 6.1 Derivation from the Weber Lagrangian

From [`WeberElectrodynamics.md`](WeberElectrodynamics.md), the Weber Lagrangian is $L = T - S$ with

$$T \;=\; \sum_i \tfrac{1}{2}m_i|\vec{v}_i|^2, \qquad S \;=\; \sum_{j<k}\frac{\kappa_{jk}q_j q_k}{r_{jk}}\!\left(1 \;+\; \frac{\dot{r}_{jk}^2}{2c^2}\right).$$

The canonical momentum is

$$\vec{p}_i \;=\; \frac{\partial L}{\partial \vec{v}_i} \;=\; m_i \vec{v}_i \;-\; \frac{\partial S}{\partial \vec{v}_i}.$$

Only pair terms containing particle $i$ contribute to $\partial S/\partial \vec{v}_i$. Using $\dot{r}_{jk} = (\vec{r}_j - \vec{r}_k)\cdot(\vec{v}_j - \vec{v}_k)/r_{jk}$ gives $\partial \dot{r}_{jk}/\partial \vec{v}_i = (\delta_{ij} - \delta_{ik})(\vec{r}_j - \vec{r}_k)/r_{jk}$, and chain-ruling yields

$$\boxed{\;\vec{p}_i \;=\; m_i \vec{v}_i \;-\; \sum_{j \neq i}\frac{\kappa_{ij}q_i q_j}{c^2}\,\frac{\dot{r}_{ij}}{r_{ij}^2}\,(\vec{r}_i - \vec{r}_j)\;}\qquad(F)$$

For the two-body case $(F)$ inverts the $\dot{x}_1 = \partial H/\partial p_{x_1}$ equation of [`WeberElectrodynamics.md`](WeberElectrodynamics.md) exactly.

### 6.2 Newton's third-law symmetry

The pairwise correction in $(F)$ is antisymmetric under $i \leftrightarrow j$, so summing over $i$ the pair contributions cancel:

$$\sum_i \vec{p}_i \;=\; \sum_i m_i \vec{v}_i.$$

In the COM frame, canonical total momentum vanishes iff kinetic total momentum does — the velocity-dependent correction cannot produce a spurious COM drift.

### 6.3 Inverse map

The inverse $\vec{p} \mapsto \vec{v}$ is implicit. Hamilton's equations give

$$\vec{v}_i \;=\; \frac{\vec{p}_i}{m_i} \;+\; \frac{1}{m_i}\sum_{j\ne i}\frac{\kappa_{ij}q_i q_j}{c^2}\,\frac{\dot{r}_{ij}}{r_{ij}^2}\,(\vec{r}_i - \vec{r}_j),$$

with $\dot{r}_{ij}$ on the right-hand side depending on $\vec{v}_i - \vec{v}_j$. Collecting the $n(n-1)/2$ unknowns $\dot{r}_{ij}$ into a vector and the sources $\hat{r}_{ij}\cdot(\vec{p}_i/m_i - \vec{p}_j/m_j)$ into another, the system takes the form $(\mathbb{I} - A)\boldsymbol{\varrho} = \boldsymbol{\pi}$, with $A$'s entries $\mathcal{O}(\kappa q^2/(m c^2 r))$. For sub-relativistic velocities $A$ is a contraction and the solution is unique.

**Users do not need the inverse map for IC construction** — the forward map $(F)$ is explicit. The integrator handles the inverse internally.

### 6.4 Canonical and physical Hamiltonians

The "kinetic-like" term $\sum_i|\vec{p}_i|^2/(2m_i)$ is **not** the physical kinetic energy when $\vec{p}$ includes the Weber correction from $(F)$. Writing $\vec{A}_i \equiv \vec{p}_i - m_i\vec{v}_i$ (the pairwise correction in $(F)$), the double-sum
$\sum_i \vec{v}_i\cdot\vec{A}_i = \sum_{i<j}(\kappa_{ij}q_iq_j/c^2)(\dot{r}_{ij}^2/r_{ij})$ follows from the pair-antisymmetry of $\vec{A}_i$ and from $\dot{r}_{ij} = (\vec{r}_i-\vec{r}_j)\cdot(\vec{v}_i-\vec{v}_j)/r_{ij}$. Expanding $|\vec{p}_i|^2 = m_i^2|\vec{v}_i|^2 - 2m_i\vec{v}_i\cdot\vec{A}_i + |\vec{A}_i|^2$ and dropping the $\mathcal{O}(c^{-4})$ tail $\sum_i |\vec{A}_i|^2/(2m_i)$:

$$\sum_i \frac{|\vec{p}_i|^2}{2m_i} \;=\; T_{\mathrm{phys}} \;-\; \sum_{i<j}\frac{\kappa_{ij}q_i q_j\,\dot{r}_{ij}^2}{c^2\,r_{ij}} \;+\; \mathcal{O}(c^{-4}).\qquad(H_1)$$

The Legendre-transform identity of [`WeberElectrodynamics.md`](WeberElectrodynamics.md) is $H = T_{\mathrm{phys}} + U_{\mathrm{Weber}}$, expressed as a function of $(\vec{r},\vec{v})$. Re-expressed in canonical coordinates $(\vec{q},\vec{p})$ via the implicit inverse map of §6.3, $H$ is **not** equal to $\sum|\vec{p}|^2/(2m) + U_{\mathrm{Weber}}$. Substituting $(H_1)$ into $H = T_{\mathrm{phys}} + U_{\mathrm{Weber}}$ gives, to $\mathcal{O}(c^{-2})$,

$$\boxed{\;H(\vec{q},\vec{p}) \;=\; \sum_i \frac{|\vec{p}_i|^2}{2m_i} \;+\; \sum_{i<j}\frac{\kappa_{ij}q_i q_j}{r_{ij}}\!\left(1 + \frac{\dot{r}_{ij}^2}{2c^2}\right) \;+\; \mathcal{O}(c^{-4}),\;}$$

with $\dot{r}_{ij}$ understood as $\dot{r}_{ij}(\vec{q},\vec{p})$ through §6.3. Observe the **sign flip** in the velocity correction relative to the physical potential: the canonical Hamiltonian carries $+\dot{r}^2/(2c^2)$ where $U_{\mathrm{Weber}}$ carries $-\dot{r}^2/(2c^2)$. The flip absorbs the $\sum_i\vec{v}_i\cdot\vec{A}_i$ cross-term produced by the Legendre transform.

The primary numerical cross-check in Section 7 step 5 is $(H_1)$ itself: given $(\vec{r},\vec{v})$ and $\vec{p}$ from $(F)$, the residual $\sum_i|\vec{p}_i|^2/(2m_i) - T_{\mathrm{phys}} + \sum_{i<j}(\kappa_{ij}q_iq_j/c^2)(\dot{r}_{ij}^2/r_{ij})$ must vanish to the $\mathcal{O}(c^{-4})$ tail $\sum_i|\vec{A}_i|^2/(2m_i)$ — typically $\lesssim 10^{-10}$ in absolute units at sub-relativistic velocities.

### 6.5 Gotcha — stationary particles can carry canonical momentum

A particle with $\vec{v}_i = \vec{0}$ in the physical frame does **not** automatically have $\vec{p}_i = \vec{0}$ in the canonical frame. If any pair containing particle $i$ has $\dot{r}_{ij} \neq 0$, then

$$\vec{p}_i \;=\; \vec{0} \;-\; \sum_{j\ne i}\frac{\kappa_{ij}q_i q_j}{c^2}\,\frac{\dot{r}_{ij}}{r_{ij}^2}\,(\vec{r}_i - \vec{r}_j) \;\ne\; \vec{0}.$$

Always compute $\vec{p}_i$ from physical $\vec{v}_i$ via $(F)$; never shortcut with "$\vec{v}_i = \vec{0} \implies \vec{p}_i = \vec{0}$".

---

## 7. Verification checklist

Every IC built with the recipes in this document should pass this six-step check:

1. **Positions and physical velocities.** Compute $\{\vec{r}_i(0), \vec{v}_i(0)\}$ from the chosen recipe. Verify $\sum_i m_i \vec{r}_i(0) = \vec{0}$ and $\sum_i m_i \vec{v}_i(0) = \vec{0}$.

2. **Pair radial rates.** For every pair, $\dot{r}_{ij}(0) = (\vec{r}_i - \vec{r}_j)\cdot(\vec{v}_i - \vec{v}_j)/r_{ij}$, computed kinematically.

3. **Physical Hamiltonian.** Evaluate
$$T_{\mathrm{phys}} \;=\; \tfrac{1}{2}\sum_i m_i|\vec{v}_i|^2, \qquad U_{\mathrm{Weber}} \;=\; \sum_{i<j}\frac{\kappa_{ij}q_i q_j}{r_{ij}}\!\left(1 - \frac{\dot{r}_{ij}^2}{2c^2}\right).$$
Confirm $E = T_{\mathrm{phys}} + U_{\mathrm{Weber}} < 0$ (bound).

4. **Canonical momenta.** Apply $(F)$ to get $\{\vec{p}_i(0)\}$. Confirm $\sum_i \vec{p}_i(0) = \vec{0}$.

5. **Canonical-momentum consistency.** Compute the residual $X \;\equiv\; \sum_i\frac{|\vec{p}_i|^2}{2m_i} \;-\; T_{\mathrm{phys}} \;+\; \sum_{i<j}\frac{\kappa_{ij}q_i q_j\,\dot{r}_{ij}^2}{c^2\,r_{ij}}$. By $(H_1)$ of §6.4 this is the $\mathcal{O}(c^{-4})$ tail $\sum_i|\vec{A}_i|^2/(2m_i)$ with $\vec{A}_i = \vec{p}_i - m_i\vec{v}_i$; confirm $|X|/|E| \lesssim 10^{-10}$ at sub-relativistic velocities. Note the common pitfall of comparing $\sum|\vec{p}|^2/(2m) + U_{\mathrm{Weber}}$ against $T_{\mathrm{phys}} + U_{\mathrm{Weber}}$: these differ by $\sum(\kappa qq/c^2)(\dot r^2/r)$ whenever any pair has $\dot{r}_{ij}\neq 0$ and should **not** agree.

6. **Forward integration.** Integrate for several natural timescales. Confirm $|H(t) - E|/|E|$ and $|\vec{L}(t) - \vec{L}(0)|/|\vec{L}(0)|$ are within the integrator's tolerance.

---

## 8. References

- **Clemente, R. A. and Assis, A. K. T.** (1991). *Two-body problem for Weber-like interactions.* Int. J. Theor. Phys. **30**(4), 537–545. DOI [10.1007/BF00672899](https://doi.org/10.1007/BF00672899). — Source of the two-body radial energy equation $(\star)$ in Section 2.
- **Frauenfelder, U. and Weber, J.** (2024). *A mathematical description of the Weber nucleus as a classical and quantum mechanical system.* Anal. Math. Phys. **14**:31. DOI [10.1007/s13324-024-00891-5](https://doi.org/10.1007/s13324-024-00891-5). — Eq. (2.3) is equivalent to $(\star)$ in the like-charge regime.
- **Wesley, J. P.** (1990). *Weber electrodynamics, Part I. General theory, steady current effects.* Found. Phys. Lett. **3**(5), 443–469. DOI [10.1007/BF00665929](https://doi.org/10.1007/BF00665929). — Background on the Weber Lagrangian and Hamiltonian used in Section 6.
- **Assis, A. K. T.** (1994). *Weber's Electrodynamics.* Springer. DOI [10.1007/978-94-017-3670-1](https://doi.org/10.1007/978-94-017-3670-1). — Standard reference for the Weber formalism and two-body dynamics.
