# Weber Hamiltonian Correctness Finding and Remediation Plan

> **Status: confirmed physics issue; correction not yet implemented.**
>
> This document records where the companion paper and the Julia package depart
> from the Weber Lagrangian, which surrounding formulas remain correct, and the
> work required for a future correction. It is a remediation plan, not a
> description of the behavior currently implemented by the package.
>
> The `_research/` and `python-frontend/` directories are outside the scope of
> this review.

## Executive summary

The companion paper correctly writes the Weber Lagrangian and correctly derives
the conserved energy as a function of positions and **physical velocities**.
The error occurs when the paper converts that velocity-space energy into a
canonical Hamiltonian: it replaces every velocity component by $p/m$.

That replacement is invalid because the Weber Lagrangian depends on velocity.
Its canonical momentum is

$$
\vec p_i = \frac{\partial L}{\partial \vec v_i}
\ne m_i\vec v_i
$$

in general when pair radial velocities are nonzero. In an $n$-body system,
different pair corrections can cancel for a particular particle, so nonzero
$\dot r_{ij}$ is sufficient to activate the correction but does not imply that
every individual $\vec p_i-m_i\vec v_i$ is nonzero. The Julia implementation
follows the paper's substitution, so it consistently integrates the Hamiltonian
derived in the paper, but that Hamiltonian is not the Legendre transform of the
Weber Lagrangian.

This is not a one-sign typo. A correct repair must:

1. use the canonical momentum obtained from the Lagrangian;
2. solve the coupled linear relation between canonical momenta and velocities;
3. replace the canonical Hamiltonian;
4. correct both Hamilton equations, including both $1/c^2$ signs in
   $\dot{\vec p}$;
5. update code that currently interprets $\vec p/m$ as physical velocity; and
6. update internal-consistency tests for the corrected Hamiltonian and
   supplement them with independent Weber-force and Legendre-transform checks.

The current and corrected systems agree in the Coulomb limit $c\to\infty$ and
at special instants where every pair has $\dot r_{ij}=0$. They generally differ
at finite $c$ once radial motion is present. Away from a singular or
ill-conditioned momentum map, the discrepancy begins at the same $O(1/c^2)$
order as the Weber correction itself.

## Background: kinetic and canonical momentum are different

For the formulas in this section, let $k_{ij}$ denote the coefficient that
multiplies the **entire** Weber pair term:

$$
k_{ij}=q_iq_j
$$

for the companion paper and for unmodified Weber electrodynamics. The package
currently scales the entire pair by $\kappa_{ij}$, which would instead make
$k_{ij}=\kappa_{ij}q_iq_j$ in every formula below.

That package convention is an unresolved, separate model question.
[`ZollnerElectrogravitationalTheory.md`](ZollnerElectrogravitationalTheory.md)
describes the Zöllner mismatch as an additional **static** residual
interaction, while the implementation and user documentation multiply both
the static and velocity-dependent Weber terms by $\kappa_{ij}$. Under the
static-residual interpretation, $M$ below must be built with $q_iq_j$, and a
separate configuration-only potential must be added to $H$; $\kappa_{ij}$ does
not enter $M$. Under the current full-pair-scaling interpretation, it does.
The future correction must choose and document one interpretation before
changing the Zöllner code. The canonical-momentum error exists under either
interpretation.

For

$$
\vec R_{ij}=\vec r_i-\vec r_j,\qquad
r_{ij}=\lVert\vec R_{ij}\rVert,\qquad
\hat r_{ij}=\frac{\vec R_{ij}}{r_{ij}},
$$

the physical relative radial velocity is

$$
\dot r_{ij}
=\hat r_{ij}\mathbin{\cdot}(\vec v_i-\vec v_j).
$$

The Weber Lagrangian used in the paper is

$$
L(\vec q,\vec v)
=
\frac12\sum_i m_i\lVert\vec v_i\rVert^2
-
\sum_{i<j}\frac{k_{ij}}{r_{ij}}
\left(1+\frac{\dot r_{ij}^{\,2}}{2c^2}\right).
$$

Differentiating with respect to a particle velocity gives

$$
\boxed{
\vec p_i
=m_i\vec v_i-\sum_{j\ne i}\vec\alpha_{ij}
}
$$

with

$$
\boxed{
\vec\alpha_{ij}
=
\frac{k_{ij}}{c^2}
\frac{\dot r_{ij}\,\vec R_{ij}}{r_{ij}^2}
=
\frac{k_{ij}}{r_{ij}c^2}\,
\dot r_{ij}\hat r_{ij}.
}
$$

Thus:

- $m_i\vec v_i$ is the **kinetic momentum**;
- $\vec p_i=\partial L/\partial\vec v_i$ is the **canonical momentum**;
- Hamilton's equations use the canonical momentum; and
- $\vec p_i=m_i\vec v_i$ only when the summed Weber correction vanishes.

For two particles, for example,

$$
\vec p_1=m_1\vec v_1-\vec\alpha_{12},\qquad
\vec p_2=m_2\vec v_2+\vec\alpha_{12}.
$$

Solving these relations for the velocities gives

$$
\vec v_1=\frac{\vec p_1+\vec\alpha_{12}}{m_1},\qquad
\vec v_2=\frac{\vec p_2-\vec\alpha_{12}}{m_2}.
$$

These equations are **implicit** because $\vec\alpha_{12}$ contains
$\dot r_{12}$, which itself contains the physical velocities. It is incorrect
to evaluate $\vec\alpha_{12}$ by first setting $\vec v_i=\vec p_i/m_i$.

## Why an effective-mass matrix appears

The paper does not introduce a matrix $M$. That omission is not evidence that
no matrix is needed; it is a consequence of the incorrect step
$\vec v_i\mapsto\vec p_i/m_i$.

The velocity-dependent part of the Weber Lagrangian is quadratic in the
velocities. Therefore the canonical-momentum equations form a linear system.
Stack all particle velocities and momenta into vectors

$$
\vec v=(\vec v_1,\ldots,\vec v_n),\qquad
\vec p=(\vec p_1,\ldots,\vec p_n).
$$

Then the momentum relation can be written compactly as

$$
\boxed{\vec p=M(\vec q)\vec v.}
$$

The $d\times d$ particle blocks of $M$ are

$$
M_{ii}
=m_iI-\sum_{j\ne i}\gamma_{ij}\,
\hat r_{ij}\hat r_{ij}^{\mathsf T},
$$

$$
M_{ij}
=\gamma_{ij}\,
\hat r_{ij}\hat r_{ij}^{\mathsf T}
\qquad(i\ne j),
$$

where

$$
\gamma_{ij}=\frac{k_{ij}}{r_{ij}c^2}.
$$

Equivalently,

$$
M(\vec q)=D-K(\vec q),\qquad
D=\operatorname{diag}(m_1I,\ldots,m_nI),
$$

where each pair contributes only along its radial relative direction.

This matrix is only bookkeeping for the coupled momentum equations already
contained in the Lagrangian. It introduces no new physical assumption. In the
Coulomb limit, $K\to0$, so $M\to D$ and the familiar relation
$\vec v_i=\vec p_i/m_i$ returns.

Where $M(\vec q)$ is invertible,

$$
\boxed{\vec v=M(\vec q)^{-1}\vec p.}
$$

For two particles, only the relative radial component is modified; tangential
components still satisfy kinetic momentum = canonical momentum. For a general
$n$-body configuration, all pairwise radial corrections participate in one
coupled $nd\times nd$ solve.

The invertibility qualification is essential. For a two-body like-charge pair
with $k_{12}>0$, the radial Legendre map becomes singular at

$$
r_{\mathrm{critical}}
=\frac{k_{12}}{\mu c^2},
\qquad
\mu=\frac{m_1m_2}{m_1+m_2}.
$$

The future implementation must define and test its behavior near this critical
surface rather than treating the inverse as unconditionally available.

## Correct canonical Hamiltonian and equations

The Legendre transform is

$$
H(\vec q,\vec p)
=\vec p^{\mathsf T}\vec v-L(\vec q,\vec v),
\qquad
\vec v=M(\vec q)^{-1}\vec p.
$$

It simplifies to

$$
\boxed{
H_{\mathrm{Weber}}(\vec q,\vec p)
=
\frac12\vec p^{\mathsf T}M(\vec q)^{-1}\vec p
+
\sum_{i<j}\frac{k_{ij}}{r_{ij}}.
}
$$

The Coulomb term is the configuration-only potential. The Weber velocity
dependence is carried by the configuration-dependent inverse mass matrix. This
Hamiltonian is still non-separable because its quadratic momentum term also
depends on $\vec q$.

This displayed Hamiltonian covers a coefficient $k_{ij}$ that scales the whole
Weber pair. If the Zöllner correction is retained as a static residual instead,
the pure-Weber coefficient $q_iq_j$ belongs in $M$ and the chosen residual
$V_{\mathrm Z}(\vec q)$ is added separately to the configuration-only
potential.

The first canonical equation is

$$
\boxed{\dot{\vec q}=M(\vec q)^{-1}\vec p=\vec v.}
$$

For one pair, the contribution to the second canonical equation is

$$
\boxed{
\dot{\vec p}_i^{\,(ij)}
=
\frac{k_{ij}}{r_{ij}^2}
\left[
\hat r_{ij}
\left(1+\frac{3\dot r_{ij}^{\,2}}{2c^2}\right)
-
\frac{\dot r_{ij}}{c^2}(\vec v_i-\vec v_j)
\right],
}
$$

with

$$
\dot{\vec p}_j^{\,(ij)}=-\dot{\vec p}_i^{\,(ij)}
$$

and $\vec v=M^{-1}\vec p$.

The signs of both $1/c^2$ terms differ from the canonical momentum equation in
the current paper and code.

### Canonical momentum rate is not the mechanical force

For a velocity-dependent Lagrangian,

$$
\dot{\vec p}_i=\frac{\partial L}{\partial\vec r_i}
$$

is a canonical equation, while the mechanical force is
$m_i\dot{\vec v}_i$. Since

$$
\vec p_i=m_i\vec v_i-\sum_{j\ne i}\vec\alpha_{ij},
$$

the two are related by

$$
m_i\dot{\vec v}_i
=
\dot{\vec p}_i
+
\frac{d}{dt}\sum_{j\ne i}\vec\alpha_{ij}.
$$

Combining the corrected canonical equation with this identity reproduces the
Weber force law. Comparing the paper's $\dot{\vec p}$ directly with
$m_i\dot{\vec v}_i$ would conflate canonical and kinetic momentum a second time.

## Companion paper: correct and incorrect inventory

The relevant file is
[`papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex`](../papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex).

### Unaffected by this Hamiltonian defect within the audited scope

The following parts do not cause this defect. This is not a blanket
re-certification of every claim in those sections.

1. **Weber potential and Weber force.** Equations labelled `potential` and
   `force` state the intended physical interaction.
2. **The auxiliary interaction $S$ and Lagrangian $L=T-S$.** Their signs produce
   the stated Weber force through the Euler–Lagrange equation.
3. **The Euler–Lagrange calculation.** Equation `euler_lagrange` is a statement
   about the mechanical Weber force and remains valid.
4. **The Legendre-transform identity in velocity variables.** Equations `H_def`,
   `H`, and `h_box` correctly give
   $H(\vec q,\vec v)=T_{\mathrm{phys}}+U_{\mathrm{Weber}}$ before the velocities
   are eliminated.
5. **Conservation and symmetry statements.** The corrected canonical
   Hamiltonian remains autonomous, translation-invariant, and
   rotation-invariant, so energy, canonical linear momentum, and canonical
   angular momentum remain conserved.
6. **The generic symplectic-integration formulas.** Their applicability to a
   general non-separable Hamiltonian is not invalidated by the incorrect Weber
   formula supplied to the method.
7. **Notation, radial-acceleration identities, units, and the qualitative
   critical-radius discussion.** These are not created by the faulty
   $v=p/m$ substitution.

### Correct only with an explicit qualification

1. **$H=T+U$.** This is correct as an energy expressed in physical velocities.
   It is not yet a canonical function $H(\vec q,\vec p)$.
2. **The definition of $\alpha_x$.** Its algebraic form is correct only when
   $\dot r$ is computed from the physical velocities obtained from
   $M\vec v=\vec p$. It must not use $\vec p/m$.
3. **Non-separability.** The conclusion is correct, but the explanation should
   refer to $\tfrac12\vec p^{\mathsf T}M(\vec q)^{-1}\vec p$, not to the current
   pair potential after a naive substitution.
4. **Computational complexity.** Pair geometry and assembly remain
   $O(n^2)$, but evaluating $\vec v=M^{-1}\vec p$ adds a coupled linear solve.
   A generic dense direct solve is $O(n^3)$ for fixed spatial dimension, though
   the pair structure may permit better implementations.

### Incorrect and requiring correction

1. **The conversion from velocities to momenta.** The sentence following
   equation `H` says each velocity is replaced by $p/m$. This is the primary
   derivation error.
2. **The two-particle velocity equations.** Equations `xdot_two` and
   `xdot_two_p2` have the signs appropriate to the current naive Hamiltonian.
   The canonical relations require

   $$
   \dot x_1=\frac{p_{x_1}+\alpha_x}{m_1},
   \qquad
   \dot x_2=\frac{p_{x_2}-\alpha_x}{m_2},
   $$

   with $\alpha_x$ evaluated implicitly through the physical velocity.
3. **The two-particle canonical momentum equation.** Equation `pdot_two` has
   both Weber-correction signs reversed relative to the exact Legendre
   transform. It must use

   $$
   \dot p_{x_1}
   =
   \frac{q_1q_2}{r_{12}^2}
   \left[
   \frac{x_1-x_2}{r_{12}}
   \left(1+\frac{3\dot r_{12}^{\,2}}{2c^2}\right)
   -
   \frac{\dot r_{12}(\dot x_1-\dot x_2)}{c^2}
   \right].
   $$
4. **The canonical $n$-body Hamiltonian.** Equation `hamiltonian` is correct as
   velocity-space energy, but it is presented as the Hamiltonian used with
   canonical momenta. Its canonical replacement is the matrix-inverse
   Hamiltonian above.
5. **The expanded $n$-body equations.** Equations `xdot_expanded` and
   `pdot_expanded` inherit the same errors as the two-particle equations.
6. **The appendix momentum definition.** The appendix defines
   $p_{x_i}=m_i\dot x_i$ and analogously for $y,z$. These are kinetic momenta,
   not the canonical variables used by Hamilton's equations.

### Minimal future paper edit

The eventual paper correction should be brief:

1. retain the correct potential, force, Lagrangian, Euler–Lagrange calculation,
   and velocity-space equation $H=T+U$;
2. identify that equation explicitly as energy temporarily expressed in
   physical velocities;
3. replace the $v=p/m$ sentence with $\vec p=M(\vec q)\vec v$;
4. state the canonical matrix Hamiltonian, define the blocks of $M$, and
   qualify $M^{-1}$ by its invertibility while referencing the existing
   critical-radius discussion;
5. correct the two-body and $n$-body $\dot q$ and $\dot p$ equations;
6. replace the appendix's kinetic-momentum definitions with the canonical
   momentum relation;
7. make the small non-separability and complexity qualifications above; and
8. bump the paper version from `1.2` to `1.3` under the repository's paper
   versioning convention.

No unrelated exposition or integrator derivation needs to be rewritten.

## Formula verifier: correct and incorrect inventory

The relevant file is
[`papers/Computational-Weber-Electrodynamics/verify_formulas.py`](../papers/Computational-Weber-Electrodynamics/verify_formulas.py).

### Checks that remain valid

- Group A: velocity-space potential identities.
- Group B: Euler–Lagrange derivation of the Weber force.
- Group C.1: canonical momentum
  $p_{x_1}=m_1\dot x_1-\alpha_x$.
- Group C.2: the Legendre transform evaluated in velocity variables gives
  $T+U$.
- Group F: radial-velocity and radial-acceleration identities.

### Checks that currently verify the wrong Hamiltonian

- The global definitions `rdot12_p`, `H_qp`, `H3`, and the corresponding
  three-particle radial rates set $\vec v=\vec p/m$.
- Group C.3 verifies that the paper's substitution was transcribed
  consistently; it does not verify that the substitution is a Legendre
  transform.
- Groups D and E differentiate the naive two-body Hamiltonian.
- Group H verifies the naive three-body equations.
- Groups G and I prove conservation properties of the naive Hamiltonian.
  Those symmetry arguments are structurally valid, but they must be rerun
  against the corrected Hamiltonian.

### Required verifier correction

The future verifier must check:

1. $\vec p=\partial L/\partial\vec v=M\vec v$;
2. $H=\vec p^{\mathsf T}\vec v-L
   =\tfrac12\vec p^{\mathsf T}M^{-1}\vec p+V_{\mathrm{Coulomb}}$;
3. $\partial H/\partial\vec p=M^{-1}\vec p$;
4. the corrected $\dot{\vec p}$ formula;
5. translation and rotation invariance of the corrected $H$; and
6. at least one multi-particle momentum-map check so that pair coupling is not
   tested only in the two-body special case.

For two particles, a closed-form inverse can keep the symbolic expressions
manageable. For three particles, it may be better to verify
$\vec p=\partial L/\partial\vec v$ and
$\dot{\vec p}=\partial L/\partial\vec q$ in velocity variables, supplemented by
numerical checks of the matrix inverse, rather than expanding a large symbolic
inverse.

The verifier is currently a manual script and SymPy is not provisioned by the
repository's paper workflow. The future fix must run it in an environment that
explicitly supplies SymPy and record a fully passing result.

## Julia package: correct and incorrect inventory

### Incorrect and requiring correction

1. **Core Weber builder**
   ([`src/hamiltonian/builders/weber.jl`](../src/hamiltonian/builders/weber.jl)).
   It constructs $\sum p^2/(2m)$ and evaluates every radial velocity with
   $\vec p/m$. It therefore builds the paper's naive Hamiltonian.
2. **Zöllner builder**
   ([`src/hamiltonian/builders/zollner.jl`](../src/hamiltonian/builders/zollner.jl)).
   It applies $\kappa-1$ to the same naive full Weber pair. This disagrees with
   the static-residual interpretation in the Zöllner theory note. If full-pair
   scaling is retained, $M^{-1}$ depends nonlinearly on all
   $\kappa_{ij}$ and the current additive canonical correction cannot simply be
   retained. If the static-residual interpretation is chosen, the pure-Weber
   $M$ is unchanged and the extra Zöllner term remains configuration-only.
   This model choice must be resolved before implementation.
3. **Default named-term decomposition**
   ([`src/hamiltonian_system.jl`](../src/hamiltonian_system.jl)).
   Its `:weber + :zollner` construction and pair closures assume the current
   full-pair scaling and naive additive Hamiltonian. Any changed term semantics
   must also be reflected in
   [`docs/src/api/system.md`](../docs/src/api/system.md).
4. **Energy statistics**
   ([`src/statistics/energy.jl`](../src/statistics/energy.jl)).
   They call $\sum p^2/(2m)$ physical kinetic energy, compute $\dot r$ with
   $\vec p/m$, and reconstruct the compiled Hamiltonian from those quantities.
   They must instead obtain $\vec v=M^{-1}\vec p$, compute
   $T_{\mathrm{phys}}=\tfrac12\sum m_i v_i^2$, and use the physical pair energy.
5. **Force statistics**
   ([`src/statistics/forces.jl`](../src/statistics/forces.jl)).
   They infer velocities and finite-difference accelerations from $\vec p/m$.
   They must obtain physical velocities from the corrected equation
   $\dot{\vec q}=\vec v$.
6. **Live Makie phase data**
   ([`ext/WeberElectrodynamicsMakieExt.jl`](../ext/WeberElectrodynamicsMakieExt.jl)).
   Its pair radial velocity also uses $\vec p/m$.
7. **Nonzero-radial two-body initial conditions**
   ([`src/initial_conditions.jl`](../src/initial_conditions.jl)).
   The `radial_velocity` keyword denotes a physical velocity, but the helper
   currently assigns radial canonical momentum as $\mu\dot r$. For two
   particles it must use

   $$
   p_r
   =
   \left(\mu-\frac{k_{12}^{(W)}}{r_{12}c^2}\right)\dot r.
   $$

   Here
   $p_r=\hat r_{12}\mathbin{\cdot}
   \mu(\vec p_1/m_1-\vec p_2/m_2)$ is the canonical momentum conjugate to the
   relative radial coordinate, and $k_{12}^{(W)}$ is the coefficient of the
   velocity-dependent Weber pair: $q_1q_2$ under the static-residual Zöllner
   interpretation, or $\kappa_{12}q_1q_2$ if full-pair scaling is deliberately
   retained.

   Zero-radial and rigid-rotation initial conditions remain valid at their
   initial instant because the Weber momentum correction then vanishes.

### Correct in principle and not the source of the defect

1. **Generic Hamiltonian constructor.** Given a correct symbolic $H$, it
   correctly forms Hamilton's equations by differentiation.
2. **Base symmetric-projection integrator.** It is a numerical method for a
   general non-separable Hamiltonian. It integrates the current wrong
   Hamiltonian faithfully; the physics error is in the Hamiltonian supplied to
   it.
3. **Pure kinetic and Coulomb custom builders.** Their standard canonical
   formulas do not contain the Weber velocity coupling.
4. **Canonical linear and angular momentum statistics.** Summing canonical
   $\vec p_i$ and $\vec r_i\times\vec p_i$ remains the appropriate Noether
   diagnostic.
5. **Trajectory extraction and generic persistence.** Copying stored
   coordinates/momenta does not itself assume $\vec p=m\vec v$.
6. **LC/KS coordinate maps.** Expressions resembling
   $\mu(p_i/m_i-p_j/m_j)$ in those maps are canonical center-of-mass coordinate
   transformations, not claims that $p_i/m_i$ is a physical velocity.

### Requires revalidation before being declared unaffected

1. **Regularized pair splitting.** The coordinate maps are canonical, but the
   corrected inverse-matrix Hamiltonian is not pair-additive in the same way as
   the current Hamiltonian. The existing full-system-minus-isolated-pair
   derivative design appears able to consume corrected compiled equations, so
   source changes are not currently identified; pair isolation and resulting
   fixtures still require revalidation.
2. **Collision callbacks.** Their canonical transformations may remain usable,
   and no source defect has been identified in this audit, but collision
   behavior and fixtures must be rechecked against the corrected singular
   dynamics.
3. **Plots extension.** Most plotting code consumes statistics without
   re-deriving physics. Its source can remain unchanged if the public
   `EnergyData` and `PairForceData` field meanings remain truthful; labels and
   plotted decompositions must still be revalidated. The Makie extension is
   separately listed above because it directly computes $\dot r$ from
   $\vec p/m$.
4. **Archive compatibility.** The serialization mechanism is generic, but old
   saved trajectories represent the old dynamical system and must not be
   presented as corrected Weber results.

## Theory and documentation inventory

### Incorrect or internally inconsistent

1. [`theory/WeberElectrodynamics.md`](WeberElectrodynamics.md) has the same
   problem. Its opening definitions $p=m\dot q$ are kinetic momenta, not the
   canonical momenta of its later Weber Lagrangian. The Lagrangian and the
   velocity-space energy $H=T+U$ are valid. Its later $+\alpha$ and $-\alpha$
   velocity equations have the correct signs for the **implicit** inverse
   momentum relation, but they contradict the opening definitions and must
   obtain $\dot r$ from the physical velocity solving $M\vec v=\vec p$ rather
   than from $\vec p/m$. Its $\dot p$ equation has both $1/c^2$ correction
   signs wrong. The Hamiltonian section and canonical equations therefore need
   the same matrix-Hamiltonian correction as the paper.
2. [`theory/InitialConditions.md`](InitialConditions.md) calls $m\vec v$
   canonical momentum and states the naive canonical Hamiltonian. Most of its
   zero-radial-velocity constructions remain numerically valid at the initial
   instant, but the general definitions must be corrected.
3. [`theory/NonZeroRadialVelocityBoundICs.md`](NonZeroRadialVelocityBoundICs.md)
   already records the important forward map
   $\vec p=\partial L/\partial\vec v$ and its implicit inverse. Its exact
   Hamiltonian discussion should be aligned with the $M^{-1}$ form, and the
   prose defining its auxiliary correction vector must be checked for a sign
   inconsistency.
4. [`theory/ZollnerElectrogravitationalTheory.md`](ZollnerElectrogravitationalTheory.md)
   correctly treats the Zöllner mismatch as a static residual, but its
   theory-level expression
   $H=\sum_i\lVert p_i\rVert^2/(2m_i)+U_W+U_g$ repeats the canonical-momentum
   error when $U_W$ is the velocity-dependent Weber interaction. Its corrected
   form should use the pure-Weber inverse-mass Hamiltonian plus the
   configuration-only $U_g$. Separately,
   [`docs/src/zollner.md`](../docs/src/zollner.md) and
   [`docs/src/quickstart.md`](../docs/src/quickstart.md) describe
   $\kappa_{ij}$ as scaling the whole velocity-dependent Weber pair. These
   sources must be reconciled before deciding whether $\kappa_{ij}$ enters
   $M$.
5. [`docs/src/api/system.md`](../docs/src/api/system.md) documents the current
   per-pair `:weber` and `:zollner` named-term closures. It must be updated if
   their fields or semantics change.
6. Other documentation describing energy as
   `kinetic + pair potential` or radial velocity as a function of $p/m$ must be
   updated when the statistics API is corrected.

### Examples and generated artifacts

1. The input momenta in examples that begin with every $\dot r_{ij}=0$ can
   often remain unchanged.
2. Subsequent finite-$c$ **Weber** trajectories were generated by the current
   Hamiltonian and must be regenerated. Unaffected custom-Coulomb sections need
   not be regenerated solely because of this finding.
3. Markdown cells that display the current Hamiltonian or describe
   $\vec v=\vec p/m$ must be corrected.
4. Stored outputs for affected Weber runs and the three committed
   `examples/figures/*.png` artifacts must be regenerated rather than retained
   as evidence for the corrected dynamics.

## Tests and regression data

The existing tests show that the implementation is internally consistent with
the current Hamiltonian. They do not establish that this Hamiltonian is the
Legendre transform of the Weber Lagrangian.

The future correction must update at least:

- `test/test_utils.jl`, whose manual Weber energy uses $\vec v=\vec p/m$;
- `test/test_hamiltonian_system.jl`, which checks current compiled formulas;
- `test/test_builders.jl`, which checks the current Weber/Zöllner additive
  decomposition;
- `test/test_named_term.jl` and `test/test_statistics.jl`, which check the
  current pair-energy decomposition;
- `test/test_initial_conditions.jl`, which checks conversion from physical
  radial velocity to canonical momentum;
- `test/test_zollner.jl`, whose expected behavior depends on resolving the
  static-residual versus full-pair-scaling model choice;
- `test/test_physics.jl`, whose energy test conserves the same Hamiltonian that
  the implementation integrates; and
- the regression fixtures under `test/regression/`, which record trajectories
  from the old system.

The corrected suite needs independent tests that do not derive expected values
from the implementation under test. Regularization and collision callback tests
also require revalidation and fixture updates, even though this audit has not
identified a necessary source edit in those subsystems.

## Future remediation plan

No partial sign-only patch should be merged. The work should proceed in the
following order.

### Phase 1: establish the corrected mathematics

1. Add a concise canonical-momentum derivation to the paper.
2. Define $M(\vec q)$ and state the exact canonical Hamiltonian.
3. Derive both Hamilton equations from that Hamiltonian.
4. State the invertibility domain and the critical-radius singularity.
5. Update `verify_formulas.py` to verify the corrected formulas independently.
6. Compile the paper and bump its version from `1.2` to `1.3`.

### Phase 2: correct the core Julia Hamiltonian

1. Implement construction of $M(\vec q)$ and the quadratic form
   $\tfrac12\vec p^{\mathsf T}M^{-1}\vec p$.
2. Route the default `HamiltonianSystem(n, dims)` through the corrected builder.
3. Decide and document how singular or ill-conditioned $M$ is handled.
4. Resolve the Zöllner model mismatch first: either retain full-pair scaling
   and derive the resulting $\kappa_{ij}$-dependent $M$, or follow the theory
   note and add a static residual to a pure-Weber $M$. Then make theory, code,
   documentation, initial conditions, and tests use that same choice.
5. Redesign named-term and pair decompositions so that they do not assume an
   additive canonical velocity correction; and
6. Provide one shared, tested way to obtain physical velocity from
   $(\vec q,\vec p)$.

### Phase 3: correct dependent APIs

1. Use physical velocity in energy, force, and live phase-space diagnostics.
2. Correct nonzero-radial-velocity initial-condition conversion.
3. Audit regularization and collision behavior against the corrected
   Hamiltonian;
4. Preserve public data structures where their existing meanings can remain
   truthful, and explicitly migrate meanings that cannot; and
5. Update theory and API documentation.

### Phase 4: replace validation artifacts

1. Add symbolic checks of the Legendre transform and both Hamilton equations.
2. Add numerical finite-difference checks of gradients of the corrected
   Hamiltonian.
3. Compare package trajectories with an independent integration of the
   implicit Weber force law and demonstrate convergence under timestep
   refinement.
4. Test the Coulomb limit and zero-radial-velocity limit.
5. Test two- and three-particle cases in 1D, 2D, and 3D, including nontrivial
   $\kappa_{ij}$.
6. Regenerate regression fixtures, notebooks, and figures only after the
   corrected tests pass.
7. Run the full Julia suite, formula verifier, paper build, and documentation
   build.

## Acceptance criteria for closing this finding

This finding is resolved only when all of the following hold:

1. The paper no longer substitutes $\vec v=\vec p/m$ in the Weber interaction.
2. The paper, verifier, theory docs, and Julia builder use the same canonical
   momentum and Hamiltonian.
3. Zöllner theory, builders, user documentation, named terms, and tests agree
   on whether $\kappa_{ij}$ scales only a static residual or the whole Weber
   pair.
4. Compiled $\dot q$ equals $M^{-1}p$.
5. Compiled $\dot p$ matches the corrected canonical equation.
6. Combining the canonical equations reproduces the stated mechanical Weber
   force.
7. Diagnostics obtain physical velocity from $\dot q$, not from $\vec p/m$.
8. Independent trajectory comparisons converge to the Weber force-law
   reference as the timestep is refined.
9. The behavior near singular $M$ is tested and documented.
10. All formula, paper, Julia, documentation, and regression checks pass.
11. Affected Weber notebook outputs and generated figures have been regenerated
    from the corrected system.

Only after these criteria are satisfied should this planning note be removed.
