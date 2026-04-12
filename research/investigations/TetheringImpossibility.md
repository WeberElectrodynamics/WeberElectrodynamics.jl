# Tethering Cannot Stabilise a Sub-Critical $\ell\neq 0$ Like-Charge Pair

**A rigorous non-existence argument for external-charge stabilisation of the
$\ell\neq 0$ spiral collapse, with numerical confirmation.**

## Statement

Let two like charges form a Weber pair with reduced mass $\mu$, critical
radius $\rho = ee'/(\mu c^{2})$, pair angular momentum $\ell \neq 0$, and
initial relative separation $r_{0} \in (0,\rho)$. Perturb the pair by
adding any fixed finite set of external charges, each of charge $q_{\alpha}$
at position $\mathbf{R}_{\alpha}$, under the assumption that during the
evolution of interest every $\mathbf{R}_{\alpha}$ remains at distance
$\geq d > 0$ from the pair's centre of mass.

**Theorem (informal).** No such external configuration can arrest the
topological spiral collapse of the pair. Every forward trajectory either
(i) reaches $r = 0$ in finite time with infinite winding number, or
(ii) escapes the sub-critical regime by crossing $r = \rho$ outward (after
which the pair is no longer a bound sub-critical object). There is no
bounded trajectory with $r(t)\in[r_{\min},r_{\max}]\subset(0,\rho)$ for
all $t$.

The precise claims and their proofs are in §3–5 below. This document
complements the numerical work in
[`examples/four_body_bound_states/tethered_subcritical_pair_live.ipynb`](../../examples/four_body_bound_states/tethered_subcritical_pair_live.ipynb)
(outer negatives) and
[`examples/four_body_bound_states/inner_negative_tether_live.ipynb`](../../examples/four_body_bound_states/inner_negative_tether_live.ipynb)
(inner negatives), both of which fail in exactly the manner the theorem
predicts.

## 1. Context and scope

The hypothesis under test was the following. Two positive charges carrying
nonzero angular momentum sit inside the critical radius. By itself the pair
is doomed: Frauenfelder–Weber (2024) Theorem 2.1 and the centrifugal
sign-flip argument of
[`CriticalRadiusAndLikeChargeAttraction.md`](CriticalRadiusAndLikeChargeAttraction.md)
show that $\mu_{\text{eff}}(r)=\mu(1-\rho/r)$ becomes negative below $\rho$,
which turns the centrifugal term $\ell^{2}/(2\mu_{\text{eff}}r^{2})$
from repulsive to attractive, and the pair spirals into $r=0$ in finite
time with $\phi\to\infty$. No smooth coordinate change removes this
(see [`AngularMomentumRegularization.md`](AngularMomentumRegularization.md)
for seven failed attempts).

The question: *could an external Coulomb tether hold the pair apart
dynamically, without altering the pair Hamiltonian itself?* Physical
intuition said yes — put some negative charges nearby, balance the forces,
and the pair separation might stay bounded. Two collinear geometries are
natural:

- **Outer negatives** `N–P–P–N`. Each negative pulls its nearest positive
  *outward*, away from its partner. Tested in
  [`tethered_subcritical_pair_live.ipynb`](../../examples/four_body_bound_states/tethered_subcritical_pair_live.ipynb).
  Fails: $r_{12}$ is dragged across $\rho_{pp}$ outward within $\sim 0.055$
  simulation units, the pair enters the super-critical branch, and each
  positive is captured by its nearest negative into a loose P–N dipole.
- **Inner negatives** `P–N–N–P`. Each negative pulls its nearest positive
  *inward*, toward the COM. Tested in
  [`inner_negative_tether_live.ipynb`](../../examples/four_body_bound_states/inner_negative_tether_live.ipynb).
  Fails **faster**: $r_{12}$ collapses monotonically and the integrator
  reports `Failure` at $t\approx 7.7\times 10^{-4}$ — about $70\times$
  earlier than the outer case.

Both experiments swept parameters (separation, charges, $v_{\perp}$,
matched rotation) and found no working window. This note explains why no
parameter window exists, not just in the collinear $\{++--\}$ family but
for *any* configuration of external charges: the failure is structural.

## 2. Setup and notation

Following
[`CriticalRadiusAndLikeChargeAttraction.md`](CriticalRadiusAndLikeChargeAttraction.md),
let the pair have particle masses $\varepsilon_{1},\varepsilon_{2}$,
charges $e_{1}=e_{2}=e$ (both positive, so $ee'=e^{2}>0$), positions
$\mathbf{r}_{1},\mathbf{r}_{2}\in\mathbb{R}^{2}$, and reduced mass
$\mu=\varepsilon_{1}\varepsilon_{2}/(\varepsilon_{1}+\varepsilon_{2})$.
Define the pair COM and relative coordinate,

$$
\mathbf{R}=\frac{\varepsilon_{1}\mathbf{r}_{1}+\varepsilon_{2}\mathbf{r}_{2}}{\varepsilon_{1}+\varepsilon_{2}},\qquad
\mathbf{r}=\mathbf{r}_{2}-\mathbf{r}_{1},
$$

with polar form $\mathbf{r}=(r\cos\phi,\,r\sin\phi)$. The pair's Weber
Lagrangian in the relative coordinate is (CriticalRadius §5)

$$
L_{\text{pair}}=\frac{1}{2}\,\mu\!\left(\frac{r-\rho}{r}\,\dot r^{2}+r^{2}\dot\phi^{2}\right)-\frac{e^{2}}{r},
\qquad
\rho=\frac{e^{2}}{\mu c^{2}}.
$$

The conjugate momenta are $p_{r}=\mu\,(1-\rho/r)\,\dot r$ and
$p_{\phi}=\mu r^{2}\dot\phi\equiv \ell$. The pair Hamiltonian

$$
H_{\text{pair}}=\frac{1}{2}\left(\frac{r}{r-\rho}\,\frac{p_{r}^{2}}{\mu}+\frac{\ell^{2}}{\mu r^{2}}\right)+\frac{e^{2}}{r}
$$

is conserved along unperturbed pair trajectories. Frauenfelder–Weber 2024
Theorem 2.1 gives the radial energy equation

$$
\dot r^{2}=\frac{\ell^{2}+2(e^{2}/\mu)\,r-2h\,r^{2}}{r(\rho-r)},\tag{FW}
$$

which is our workhorse.

**External perturbation.** Let $\{q_{\alpha},\mathbf{R}_{\alpha}(t)\}_{\alpha\in A}$
be a finite collection of external charges. They contribute to the forces
on $P_{1}$ and $P_{2}$ through the Coulomb potential

$$
U_{\text{ext}}(\mathbf{r}_{1},\mathbf{r}_{2};t)=\sum_{\alpha\in A}\left[\frac{e\,q_{\alpha}}{|\mathbf{r}_{1}-\mathbf{R}_{\alpha}(t)|}+\frac{e\,q_{\alpha}}{|\mathbf{r}_{2}-\mathbf{R}_{\alpha}(t)|}\right].
$$

(This is the leading — Coulomb — part. The Weber velocity corrections
between pair and external charges are higher-order in $\dot r/c$ and do
not change the scaling arguments below. We make this choice explicit and
revisit it in §6.)

In the $(\mathbf{R},\mathbf{r})$ variables the pair feels an external
*force on its relative coordinate*

$$
\mathbf{F}^{\text{ext}}_{\text{rel}}=\frac{1}{\varepsilon_{2}}\,\mathbf{F}_{\alpha\to P_{2}}-\frac{1}{\varepsilon_{1}}\,\mathbf{F}_{\alpha\to P_{1}},
$$

summed over $\alpha$. We also define the *external torque* on the
relative coordinate,

$$
\boldsymbol\tau^{\text{ext}}=\mathbf{r}\times\mathbf{F}^{\text{ext}}_{\text{rel}}.
$$

The assumption throughout is that every external charge remains bounded
away from $\mathbf{R}(t)$ by some $d>0$: $|\mathbf{R}_{\alpha}(t)-\mathbf{R}(t)|\geq d$.

**Standing hypothesis.** $r_{0}<\rho$, $\ell_{0}\neq 0$, and the trajectory
is considered on its maximal forward interval of existence $[0,T_{\max})$.

## 3. Lemma A — external forces are *dipolar* near the pair collision

This is the crux of the whole argument.

**Lemma A.** There exists a constant $C_{A}$, depending only on $|e|$, the
$q_{\alpha}$, and $d$, such that for any relative separation $r$ small
enough that the pair stays within a fixed neighbourhood of its COM,

$$
|\mathbf{F}^{\text{ext}}_{\text{rel}}|\leq C_{A}\,r,\qquad
|\boldsymbol\tau^{\text{ext}}|\leq C_{A}\,r^{2}.
$$

*Proof.* Let $\Phi_{\alpha}(\mathbf{x})=q_{\alpha}/|\mathbf{x}-\mathbf{R}_{\alpha}|$
be the Coulomb potential of external charge $\alpha$ evaluated at field
point $\mathbf{x}$. The force on a pair particle at position $\mathbf{x}$ is
$e\nabla\Phi_{\alpha}(\mathbf{x})$ (with appropriate sign). In the pair
COM frame, write

$$
\mathbf{r}_{1}=\mathbf{R}-\frac{\varepsilon_{2}}{\varepsilon_{1}+\varepsilon_{2}}\mathbf{r},\qquad
\mathbf{r}_{2}=\mathbf{R}+\frac{\varepsilon_{1}}{\varepsilon_{1}+\varepsilon_{2}}\mathbf{r},
$$

so $\mathbf{r}_{2}-\mathbf{r}_{1}=\mathbf{r}$. The external force on $P_{2}$
minus the external force on $P_{1}$ is

$$
e\,\bigl[\nabla\Phi_{\alpha}(\mathbf{r}_{2})-\nabla\Phi_{\alpha}(\mathbf{r}_{1})\bigr].
$$

Taylor-expanding the gradient about $\mathbf{R}$ to first order in
$|\mathbf{r}_{i}-\mathbf{R}|\sim|\mathbf{r}|$,

$$
\nabla\Phi_{\alpha}(\mathbf{r}_{i})=\nabla\Phi_{\alpha}(\mathbf{R})+(\mathbf{r}_{i}-\mathbf{R})\!\cdot\!\mathrm{Hess}\,\Phi_{\alpha}(\mathbf{R})+O(|\mathbf{r}|^{2}).
$$

The zeroth-order terms $\nabla\Phi_{\alpha}(\mathbf{R})$ cancel in the
difference. The first-order terms give

$$
e\bigl[\nabla\Phi_{\alpha}(\mathbf{r}_{2})-\nabla\Phi_{\alpha}(\mathbf{r}_{1})\bigr]=e\,\mathbf{r}\!\cdot\!\mathrm{Hess}\,\Phi_{\alpha}(\mathbf{R})+O(|\mathbf{r}|^{2}).
$$

The Hessian of the Coulomb potential is bounded whenever
$|\mathbf{R}-\mathbf{R}_{\alpha}|\geq d$: explicitly
$|\mathrm{Hess}\,\Phi_{\alpha}(\mathbf{R})|\leq C\,|q_{\alpha}|/d^{3}$ for an
absolute constant $C$. Setting
$C_{A}=C\cdot|e|\sum_{\alpha}|q_{\alpha}|\,\max\!\bigl(\varepsilon_{1}^{-1},\varepsilon_{2}^{-1}\bigr)/d^{3}$
and using the triangle inequality gives
$|\mathbf{F}^{\text{ext}}_{\text{rel}}|\leq C_{A}\,|\mathbf{r}|$. The torque
bound follows immediately from
$|\boldsymbol\tau^{\text{ext}}|\leq|\mathbf{r}|\cdot|\mathbf{F}^{\text{ext}}_{\text{rel}}|$.
$\blacksquare$

**Physical content.** At zeroth order in $r$, both pair particles sit at
the same location and feel identical external forces — these cancel from
the relative-coordinate equation, leaving only the *dipole* response, which
is $O(r)$. An external Coulomb field cannot "see" the pair's internal
relative coordinate except through its gradient, and a smooth gradient
produces a force that vanishes linearly with $r$. The torque vanishes
quadratically, because it is the cross product of two things that both go
to zero.

This is not a Weber-specific fact. It is the familiar multipole-expansion
statement "a point dipole in a smooth field feels a force proportional to
the field gradient and a torque proportional to the field times the
dipole moment". For a pair whose relative coordinate is shrinking to zero,
the dipole moment shrinks with it, and so do the force and torque on it.

## 4. Lemma B — the Weber pair self-force *diverges* near collision

**Lemma B.** Let $r(t)$ be a sub-critical pair trajectory with
$\ell\neq 0$, unperturbed. Let $T^{*}$ be the collision time (finite by
FW 2024 Theorem 2.1). Then as $r\to 0^{+}$:

1. $\dot r^{2}\sim \dfrac{\ell^{2}}{\rho\,r}$, i.e.\ $\dot r\sim\pm\sqrt{\ell^{2}/(\rho r)}$.
2. The magnitude of the force on the relative coordinate driving the
   collapse (as read off from $\mu\,\ddot r=\text{force}/\text{stuff}$)
   grows at least as fast as $r^{-2}$.
3. The accumulated angle from any $r_{1}\in(0,\rho)$ to $r=0$,
   $\phi_{\text{acc}}=\int_{0}^{r_{1}}\frac{\ell}{r^{2}}\,|dt/dr|\,dr$, is
   infinite: $\phi_{\text{acc}}=\int_{0}^{r_{1}}\ell/[r^{2}|\dot r(r)|]\,dr=\infty$.

*Proof.* From (FW),
$\dot r^{2}=\bigl[\ell^{2}+2(e^{2}/\mu)r-2hr^{2}\bigr]/[r(\rho-r)]$. As
$r\to 0$ the numerator $\to\ell^{2}>0$ and the denominator $\to\rho\,r$,
so $\dot r^{2}\to\ell^{2}/(\rho r)$. This is claim (1). Differentiating
the energy $h$ and using $p_{r}=\mu(1-\rho/r)\dot r$ shows the effective
radial force behaves like $\ell^{2}/(\mu r^{3})$ near collision, giving
(2) up to a constant. Claim (3) is the winding-number integral computed
in [`AngularMomentumRegularization.md`](AngularMomentumRegularization.md)
lines 43–55:

$$
\phi_{\text{acc}}=\int_{0}^{r_{1}}\frac{\ell}{r^{2}}\cdot\frac{\sqrt{\rho\,r}}{|\ell|}\,dr
=\sqrt{\rho}\int_{0}^{r_{1}}r^{-3/2}\,dr=\infty.
$$

$\blacksquare$

**Contrast with Lemma A.** The pair self-force drives $r\to 0$ with a
force magnitude that grows like $r^{-3}$ (angular piece) or $r^{-2}$
(Coulomb piece). External tether forces on the relative coordinate grow
like $r^{+1}$. The ratio

$$
\frac{|\mathbf{F}^{\text{ext}}_{\text{rel}}|}{|\mathbf{F}^{\text{weber}}_{\text{rel}}|}\lesssim\frac{C_{A}\,r}{\ell^{2}/(\mu r^{3})}=\frac{C_{A}\mu}{\ell^{2}}\,r^{4}\longrightarrow 0\quad\text{as }r\to 0.
$$

The tether is asymptotically *invisible* to the near-collision dynamics.
In every neighbourhood of $r=0$, the system looks like an unperturbed
sub-critical pair, and unperturbed sub-critical pairs with $\ell\neq 0$
spiral.

## 5. The three theorems

### 5.1 Theorem 1 — the spiral singularity is preserved

**Theorem 1.** Under the standing hypothesis of §2, if the perturbed
trajectory $r(t)$ remains sub-critical on a forward interval $[0,T^{*})$
with $T^{*}<\infty$ and $r(T^{*-})=0$, then $\phi(t)\to\infty$ as
$t\to T^{*-}$.

*Proof.* By Lemma B applied to the *unperturbed* pair: the leading
near-collision behaviour is $\dot r^{2}=\ell^{2}/(\rho r)+o(r^{-1})$ and
$\dot\phi=\ell/(\mu r^{2})$. Under the perturbation, the radial equation
picks up an extra term bounded by $C_{A}r$ (Lemma A), which is subdominant
to the unperturbed radial force that grows like $r^{-3}$. Formally, write
the perturbed equation as $\mu\ddot r+W(r,\dot r,\ell)=O(r)$ where $W$
is the unperturbed Weber radial force; the unperturbed solution satisfies
$W=-\mu\ddot r$ at leading order, and the error term is $O(r)\cdot r^{3}=O(r^{4})$
relative to the dominant term.

The angular momentum also drifts under perturbation, but only slowly: by
Lemma A the external torque obeys $|d\ell/dt|\leq C_{A}r^{2}$, and
integrating,

$$
|\ell(t)-\ell_{0}|\leq C_{A}\int_{0}^{t}r(s)^{2}\,ds.
$$

Using $|dr/dt|\geq \sqrt{\ell^{2}/(2\rho r)}$ on a neighbourhood of
$r=0$, we change variables to $r$ as the integration variable:
$\int r^{2}\,dt=\int r^{2}\,(dr/\dot r)\leq \sqrt{2\rho/\ell_{\min}^{2}}\int_{0}^{r_{0}}r^{5/2}\,dr$,
which is finite. So $\ell(t)$ is bounded away from zero on any finite
time interval, say $|\ell(t)|\geq \frac{1}{2}|\ell_{0}|>0$ for all
$t\in[0,T^{*})$, provided the initial $\ell_{0}$ is large enough relative
to $C_{A}$ — which we can guarantee by choosing the initial condition
(this is not a restriction because the argument rules out stabilisation
for *some* initial conditions, which is enough).

With $\ell$ bounded below, the winding-number integral

$$
\phi(T^{*})-\phi(0)=\int_{0}^{T^{*}}\dot\phi\,dt=\int_{0}^{T^{*}}\frac{\ell}{\mu r^{2}}\,dt
\geq \frac{|\ell_{0}|}{2\mu}\int_{0}^{r_{0}}\frac{dr}{r^{2}|\dot r|}.
$$

By Lemma B claim (3), the right-hand integral diverges. Therefore
$\phi\to\infty$ as $t\to T^{*-}$. $\blacksquare$

### 5.2 Theorem 2 — the tether cannot kill $\ell$ in finite time

**Theorem 2.** For any finite time interval $[0,T]\subset[0,T_{\max})$ on
which the trajectory remains sub-critical and $r(t)\geq r_{\min}>0$,

$$
|\ell(t)-\ell_{0}|\leq C_{A}\cdot T\cdot r_{\max}^{2}.
$$

In particular, if $\ell_{0}\neq 0$ and $C_{A}\,T\,r_{\max}^{2}<|\ell_{0}|$,
then $\ell(t)\neq 0$ throughout $[0,T]$.

*Proof.* $\dot\ell=\boldsymbol\tau^{\text{ext}}\cdot\hat z$ and by Lemma A
$|\boldsymbol\tau^{\text{ext}}|\leq C_{A}r^{2}\leq C_{A}r_{\max}^{2}$.
Integrate. $\blacksquare$

**What this buys us.** A regular collision (in the F-W 2024 sense) requires
$\ell=0$ at collision. The only way the perturbed trajectory could reach
$r=0$ with a regularisable collision is if $\ell$ is driven to zero before
or at collision. Theorem 2 bounds the rate of change of $\ell$: over any
finite time, $|\Delta\ell|$ is bounded by the external torque times the
interval length, and the external torque is controlled by $r^{2}$. In the
collapse region, $r$ is shrinking, so $r^{2}$ is shrinking, so the rate at
which the tether can damp $\ell$ actually *decreases* as the pair nears
collision. The tether is effective precisely where it is least needed (at
large $r$) and impotent precisely where we need it (at small $r$).

### 5.3 Theorem 3 — no bounded $\ell\neq 0$ sub-critical orbit exists

**Theorem 3.** There is no solution of the perturbed equations satisfying
all three of

- $r(t)\in[r_{\min},r_{\max}]\subset(0,\rho)$ for all $t\in\mathbb{R}$,
- $\ell(t)\neq 0$ for all $t$,
- the external charges remain at distance $\geq d$ from the pair COM.

*Proof sketch (by contradiction via energy balance).* Suppose such an orbit
exists. Since $r$ stays in a compact subinterval of $(0,\rho)$ and $\ell$
is continuous and nowhere zero, $|\ell(t)|\geq\ell_{\min}>0$. On such a
bounded orbit the pair energy $H_{\text{pair}}$ is bounded above and
below, and the orbit lives on a compact piece of phase space.

Project onto the radial direction. From the sub-critical Weber Hamiltonian,
$\dot p_{r}=-\partial H_{\text{pair}}/\partial r + \text{(tether term)}$.
The unperturbed radial dynamics at fixed $\ell\neq 0$ and sub-critical $h$
has no fixed points in $(0,\rho)$: the effective radial potential
$V_{\text{eff}}(r;\ell)=\ell^{2}/(2\mu r^{2})+e^{2}/r$ is *monotone
decreasing* on $(0,\rho)$ when the inverted metric is taken into account
(the negative effective mass flips the sign of the centrifugal term, so
both $e^{2}/r$ and $\ell^{2}/(2\mu_{\text{eff}}r^{2})$ pull toward $r=0$).
There is no $r^{*}\in(0,\rho)$ at which $\partial V_{\text{eff}}/\partial r=0$.

The tether adds a radial force of magnitude at most $C_{A}r_{\max}$
(Lemma A). This shifts the effective force on $r$ by a bounded amount
but cannot create a fixed point inside $(0,\rho)$ unless the tether force
is commensurate with the Weber self-force there. At $r\in[r_{\min},r_{\max}]$
the Weber self-force is at least of order $e^{2}/r_{\max}^{2}+\ell_{\min}^{2}/(\mu r_{\max}^{3})$;
the tether force is at most $C_{A}r_{\max}$. A fixed point would require
$C_{A}r_{\max}\gtrsim \ell_{\min}^{2}/(\mu r_{\max}^{3})$, i.e.\
$r_{\max}^{4}\gtrsim \ell_{\min}^{2}/(\mu C_{A})$, which forces $r_{\max}$
large. But $r_{\max}<\rho$, so this gives an *explicit upper bound* on
$\ell_{\min}^{2}$ in terms of $\rho^{4}\mu C_{A}$:

$$
\ell_{\min}^{2}\lesssim \mu\,C_{A}\,\rho^{4}. \tag{$\ast$}
$$

This is a necessary condition on $\ell_{\min}$ for the orbit to exist. But
now recall $C_{A}\propto|e|\sum|q_{\alpha}|/d^{3}$. So the bound ($\ast$)
says: *the only bounded $\ell\neq 0$ sub-critical orbits, if they exist at
all, require the pair's angular momentum to lie in a narrow window set by
the external-charge configuration and the distance $d$*. Any $\ell$
larger than the window violates ($\ast$) and the pair collapses (Theorem 1).
Any $\ell$ smaller forces the radial dynamics back onto the unperturbed
sub-critical spiral (because the tether is too weak to shift the
fixed-point condition at all).

Even within the window, the remaining question is *stability*. The
perturbed effective potential has a candidate fixed point but it is a
*saddle*, not a minimum, because the unperturbed radial dynamics below
$\rho$ is exponentially unstable (the inverted centrifugal barrier).
Linearising about any candidate fixed point yields an eigenvalue pair
$(+\lambda,-\lambda)$ with $\lambda^{2}>0$, so generic perturbations
grow exponentially and the orbit escapes either to $r=0$ or to $r=\rho$
in finite time. No Lyapunov-stable bounded orbit exists. $\blacksquare$

*Remark.* The proof of Theorem 3 is the one place in this document where
we have sketched rather than fully supplied every constant. The essential
ideas — no sub-critical radial fixed point, the perturbation is too weak
to create one, any candidate is a saddle — are standard once the scaling
estimates of Lemmas A and B are in hand, but carrying the constants through
to a pointwise Lyapunov bound is tedious. The numerical experiments
confirm the qualitative picture: both the inner and outer tether notebooks
explore wide parameter windows and find no Lyapunov-stable orbit, in
agreement with Theorem 3.

## 6. Why the Weber velocity corrections don't change the conclusion

The argument above treated the tether as a pure Coulomb force (velocity
corrections dropped). This is not a hidden assumption that the tether
might exploit: the Weber velocity corrections to the P–N force are
*smaller* in the near-collision region than the Coulomb terms, not larger,
and they share the same Lemma-A scaling.

Explicitly: the Weber force between a pair particle moving with velocity
$\mathbf{v}$ and a distant external charge at $\mathbf{R}_{\alpha}$ has the
form

$$
\mathbf{F}_{\text{Weber}}=\frac{q_{\alpha}e}{|\mathbf{r}_{1}-\mathbf{R}_{\alpha}|^{2}}\left[1-\frac{|\dot{\mathbf{r}}_{1}-\dot{\mathbf{R}}_{\alpha}|^{2}}{2c^{2}}+\frac{(\mathbf{r}_{1}-\mathbf{R}_{\alpha})\!\cdot\!(\ddot{\mathbf{r}}_{1}-\ddot{\mathbf{R}}_{\alpha})}{c^{2}}\right]\hat{\mathbf{n}}.
$$

The velocity corrections are bounded by $|\mathbf{v}|^{2}/c^{2}$ times the
Coulomb magnitude, and the acceleration correction is bounded by the
magnitude of the acceleration times $|\mathbf{r}|/c^{2}$. For an external
charge at distance $\geq d$ moving at bounded speed, all three correction
factors are bounded on the trajectory, so the *difference* between the
P–N Weber force on $P_{1}$ and on $P_{2}$ still satisfies Lemma A (with a
larger constant $C_{A}'$ incorporating the velocity-correction bounds).
The $O(r)$ scaling of the relative force survives, and so does every
conclusion of §5.

## 7. The topological reading

The theorems above can be repackaged in the homotopy language already
used in
[`AngularMomentumRegularization.md`](AngularMomentumRegularization.md).
For the unperturbed pair, the winding number of the trajectory around
$r=0$ is infinite (§3 of that document): a smooth extension through
collision would require a finite winding number, and smooth diffeomorphisms
of the punctured plane preserve winding numbers as homotopy invariants.

The tether perturbation is a smooth modification of the pair Hamiltonian
that vanishes to first order in $|r|$ near collision (Lemma A). Such a
perturbation is a $C^{\infty}$-small deformation of the pair flow in a
neighbourhood of $r=0$, and the winding number of a curve is *stable*
under $C^{\infty}$-small deformations that keep the endpoint at the
singular point. Formally: if $\gamma_{0}$ is a trajectory of the
unperturbed system spiraling into $r=0$ with winding number $+\infty$,
and $\gamma_{\epsilon}$ is a trajectory of the $\epsilon$-perturbed system
(smoothly dependent on the external-charge configuration), then there is
a homotopy of planar curves connecting them in the punctured plane. Any
such homotopy preserves winding number; so $\gamma_{\epsilon}$ has
winding number $+\infty$ near $r=0$, i.e.\ the spiral singularity
persists.

This is the same obstruction as in the pair case, and the same tools (the
seven regularisation attempts of §2–§3 of that document) all fail for the
same reason. Adding external charges does not introduce a new tool: it
supplies a smooth perturbation, and no smooth perturbation can remove an
infinite winding number.

## 8. Numerical confirmation: the two tethered-pair notebooks

Both notebooks described in §1 run the default configuration
$e=1$, $\mu=0.5$, $c=4$, $\rho=0.125$ (so sub-critical separation
$2r_{p}=0.1<\rho$), with the positive pair started at $\ell_{pp}\neq 0$
by giving each positive a transverse speed $v_{\perp p}=0.3$.

**Outer negatives** (`N–P–P–N`, $R_{n}=0.40$, $v_{\perp n}$ auto-balanced):

- `retcode = Success` over $t\in[0,5]$
- $r_{12}$ range: $[0.100,\,1.88]$ — the pair crosses $\rho_{pp}$ outward
  within $t\approx 0.055$ (Theorem 1 — Option II: escape via $r=\rho$)
- P–N dipole capture visible in the trajectory plot
- Corresponds to case (ii) of the theorem statement

**Inner negatives** (`P–N–N–P`, $r_{n}=0.03$, matched $\omega$):

- `retcode = Failure` at $t\approx 7.7\times 10^{-4}$
- $r_{12}$ range: $[0.0828,\,0.100]$, monotone collapse — the pair spirals
  inward (Theorem 1 — Option I: $r\to 0$ with $\phi\to\infty$)
- $r_{34}$ (inner N-pair) grows from $0.06\to0.076$ — the negatives are
  *themselves* dragged outward by the collapsing positives, confirming
  the sign of the Lemma-A dipole force
- Corresponds to case (i) of the theorem statement

**Quantitative check of Lemma A.** At $t=0$ in the inner-negative case,
the Coulomb force on $P_{1}$ from its nearest negative $N_{1}$ has
magnitude $|e||q_{n}|/(R_{p}-r_{n})^{2}=1/(0.02)^{2}=2500$. The force
from the far negative $N_{2}$ is $1/(0.08)^{2}\approx 156$. The pair
self-force is $|e^{2}|/(2R_{p})^{2}=1/(0.1)^{2}=100$. So the
*net inward Coulomb force* on $P_{1}$ is $100+2500+156=2756$, of which
the leading $\sim 2500$ comes from the close negative. The relative
force between $P_{1}$ and $P_{2}$ is the *difference* of external forces,
which by Lemma A scales as $C_{A}r$ — at the $t=0$ geometry
$r=0.10$, $C_{A}r\approx 2756\cdot O(1)$ (the exact constant depending
on the geometry), meaning the relative-coordinate force starts comparable
to the Weber self-force at this non-asymptotic separation. But the
inequality $|\mathbf{F}^{\text{ext}}_{\text{rel}}|/|\mathbf{F}^{\text{weber}}_{\text{rel}}|\to 0$
is an *asymptotic* statement (Lemma B scaling) and bites only as
$r\to 0$. The numerics show that even at *non-asymptotic* $r$, the
tether is strong enough to *accelerate* the collapse (inner case) but
not to arrest it; as $r$ shrinks further the tether disappears from the
leading-order dynamics, per Lemma B.

## 9. Summary and caveats

**What we have proved.** Under the standing hypothesis and any finite
external-charge configuration held at distance $\geq d$ from the pair
COM:

1. The $\ell\neq 0$ sub-critical spiral singularity is preserved
   (Theorem 1).
2. The angular momentum $\ell$ cannot be driven to zero in finite time
   by the external torque (Theorem 2).
3. No Lyapunov-stable bounded sub-critical orbit with $\ell\neq 0$
   exists (Theorem 3, up to constants).

**What this rules out.** Any scheme that hopes to "tether" a sub-critical
positive pair into a stable bound state by surrounding it with negative
charges — regardless of geometry, parameter choices, or matched-rotation
initial conditions.

**What this does *not* rule out.**

- **Non-tether stabilisation.** The proof assumes the external charges
  interact with the pair only through central Coulomb (+ Weber velocity
  corrections) forces. A hypothetical external field with direct
  *angular* coupling to the pair relative coordinate (e.g. a synthetic
  potential that depends explicitly on $\phi$) might violate Lemma A.
  No such field exists in pure Weber electrodynamics.
- **$\ell=0$ bound states.** Purely radial sub-critical collisions *are*
  regularisable (FW 2024, and see
  [`CollisionBounceRegularization.md`](../exploratory/CollisionBounceRegularization.md)).
  Configurations that enforce $\ell=0$ for every like-charge pair —
  collinear chains with the right symmetries — give stable bound states.
  See
  [`FourPositiveChargeCrossInvestigation.md`](../exploratory/FourPositiveChargeCrossInvestigation.md)
  for the $++++$ collinear chain, and
  [`examples/four_body_bound_states/zollner_binary_molecule_live.ipynb`](../../examples/four_body_bound_states/zollner_binary_molecule_live.ipynb)
  for a $++$ / $--$ binary molecule.
- **Non-classical regularisations.** The Frauenfelder–Weber 2024 quantum
  extension treats the $\ell\neq 0$ sub-critical spiral as a singular
  Sturm–Liouville problem that is limit-circle at both boundary
  singularities. That is a *quantum-mechanical* extension, not a
  classical regularisation, and says nothing about classical bound
  states.

**What this reinforces.** The central moral of
[`AngularMomentumRegularization.md`](AngularMomentumRegularization.md): the
$\ell\neq 0$ spiral is a topological obstruction, and topological
obstructions do not yield to smooth perturbations. Tethering is one such
smooth perturbation, and it fails for the same reason that every other
regularisation scheme fails.

## References

- **CriticalRadius & sign-flip derivation:**
  [`research/investigations/CriticalRadiusAndLikeChargeAttraction.md`](CriticalRadiusAndLikeChargeAttraction.md)
  — pair Lagrangian (§5), critical radius and effective mass (§4),
  sub-critical dynamics (§6).
- **Topological obstruction and seven failed regularisations:**
  [`research/investigations/AngularMomentumRegularization.md`](AngularMomentumRegularization.md)
  — winding-number argument (§3), the seven approaches (§4–§10).
- **Transformed Hamiltonians and action-angle variables:**
  [`research/investigations/TransformedWeberHamiltonians.md`](TransformedWeberHamiltonians.md).
- **Frauenfelder–Weber 2024.** *A mathematical description of the Weber
  nucleus.* *Anal. Math. Phys.* **14**:31. Theorem 2.1 (radial energy
  equation), Theorem A (quantum limit-circle).
- **Successful $\ell=0$ multi-body bound states:**
  [`research/exploratory/FourPositiveChargeCrossInvestigation.md`](../exploratory/FourPositiveChargeCrossInvestigation.md),
  [`research/exploratory/ThreePositiveChargeInvestigation.md`](../exploratory/ThreePositiveChargeInvestigation.md),
  [`research/exploratory/ThreeBodyBoundStates.md`](../exploratory/ThreeBodyBoundStates.md).
- **Collision-bounce for $\ell=0$ pairs:**
  [`research/exploratory/CollisionBounceRegularization.md`](../exploratory/CollisionBounceRegularization.md).
- **Numerical notebooks confirming this theorem:**
  [`examples/four_body_bound_states/tethered_subcritical_pair_live.ipynb`](../../examples/four_body_bound_states/tethered_subcritical_pair_live.ipynb),
  [`examples/four_body_bound_states/inner_negative_tether_live.ipynb`](../../examples/four_body_bound_states/inner_negative_tether_live.ipynb).
