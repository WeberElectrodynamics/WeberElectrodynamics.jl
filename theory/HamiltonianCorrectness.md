# Weber Hamiltonian Correctness Finding and Minimum Repair

> **Status: confirmed physics error; correction not yet implemented.**
>
> This note describes the behavior that the theory, paper, and package must
> eventually implement. It does not describe the equations currently integrated
> by the package.
>
> The `_research/` directory is outside the scope of this review.

## Result

The Weber Lagrangian in the theory notes and companion paper is correct. The
conserved energy written in terms of positions and physical velocities is also
correct. The error occurs when canonical momenta are introduced.

The current derivation treats

$$
\vec p_i=m_i\vec v_i
$$

as the canonical momentum of particle $i$, where $m_i$ is its mass and
$\vec v_i$ is its physical velocity. This equality gives kinetic momentum. It
does not give canonical momentum for the velocity-dependent Weber Lagrangian.

The correct definition is

$$
\vec p_i=\frac{\partial L}{\partial\vec v_i},
$$

where $L$ is the Weber Lagrangian, $\vec p_i$ is canonical momentum, and
$\vec v_i$ is the velocity with respect to which the derivative is taken.

Consequently:

1. physical velocity is generally not $\vec p_i/m_i$;
2. the velocity-space energy cannot be turned into a canonical Hamiltonian by
   replacing each velocity with momentum divided by mass;
3. the current paper and package use the wrong canonical Hamiltonian;
4. both Weber correction signs in the current canonical momentum-rate equation
   are wrong; and
5. diagnostics and initial-condition helpers that interpret $\vec p_i/m_i$ as
   physical velocity are wrong whenever a nonzero radial Weber correction is
   present.

The current and corrected dynamics agree in the Coulomb limit $c\to\infty$ and
at instants where every pair radial velocity is zero. They differ in general at
finite $c$ once radial motion develops.

This note uses the same particle, component, and pair notation as
[`WeberElectrodynamics.md`](WeberElectrodynamics.md). It uses only particle
and pair equations.

## Correct starting point

For particles $i$ and $j$, define

$$
\vec r_{ij}=\vec r_i-\vec r_j,\qquad
r_{ij}=\lVert\vec r_{ij}\rVert,\qquad
\hat r_{ij}=\frac{\vec r_{ij}}{r_{ij}}.
$$

Here $\vec r_i$ and $\vec r_j$ are the particle positions, $\vec r_{ij}$ is
their relative position, $r_{ij}$ is their separation, and $\hat r_{ij}$ is the
unit vector pointing from particle $j$ to particle $i$.

The physical relative radial velocity is

$$
\dot r_{ij}
=
\hat r_{ij}\mathbin{\cdot}(\vec v_i-\vec v_j)
=
\frac{\vec r_{ij}\mathbin{\cdot}(\vec v_i-\vec v_j)}{r_{ij}}.
$$

Here $\vec v_i$ and $\vec v_j$ are physical particle velocities and
$\dot r_{ij}$ is the time derivative of the pair separation.

In the absolute units used by the repository, the $n$-particle Weber
Lagrangian is

$$
L(\vec r,\vec v)
=
\frac12\sum_i m_i\lVert\vec v_i\rVert^2
-
\sum_{i<j}\frac{q_iq_j}{r_{ij}}
\left(1+\frac{\dot r_{ij}^{\,2}}{2c^2}\right).
$$

Here $\vec r$ and $\vec v$ denote all particle positions and velocities,
$m_i$ is the mass of particle $i$, $q_i$ and $q_j$ are electric charges, and
$c$ is the speed of light. The coefficient of a Weber pair is directly
$q_iq_j$. There is no additional per-pair coupling parameter in the current
package.

## Correct canonical momentum

Differentiating the Lagrangian with respect to a particle velocity gives

$$
\boxed{
\vec p_i
=
m_i\vec v_i
-
\sum_{j\ne i}
\frac{q_iq_j}{c^2}
\frac{\dot r_{ij}}{r_{ij}^2}
(\vec r_i-\vec r_j).
}
$$

Here $\vec p_i$ is the canonical momentum of particle $i$. Every term in the
sum is the Weber correction contributed by a pair containing particle $i$.
The correction is radial and depends on the physical radial velocity of that
pair.

This equation is linear in all particle velocities because each radial
velocity satisfies

$$
\dot r_{ij}
=
\frac{(\vec r_i-\vec r_j)\mathbin{\cdot}
(\vec v_i-\vec v_j)}{r_{ij}}.
$$

Here all symbols have the pair meanings defined above. Although the equations
are linear in the velocities, different pairs share particle velocities, so
all velocity components must generally be solved together.

### Two-particle form

For two particles, define the Weber momentum correction

$$
\vec\alpha
=
\frac{q_1q_2}{c^2}
\frac{\dot r(\vec r_1-\vec r_2)}{r^2}
=
\frac{q_1q_2}{rc^2}\dot r\,\hat r.
$$

Here $r=\lVert\vec r_1-\vec r_2\rVert$, $\hat r=(\vec r_1-\vec r_2)/r$,
$\dot r=\hat r\mathbin{\cdot}(\vec v_1-\vec v_2)$, and $\vec\alpha$ is the
pair correction.

The two canonical momenta are

$$
\boxed{
\vec p_1=m_1\vec v_1-\vec\alpha,\qquad
\vec p_2=m_2\vec v_2+\vec\alpha.
}
$$

Here $\vec p_1$ and $\vec p_2$ are canonical momenta and $m_1\vec v_1$ and
$m_2\vec v_2$ are kinetic momenta. The opposite correction signs preserve
total momentum.

Solving these equations for the physical velocities gives

$$
\boxed{
\vec v_1=\frac{\vec p_1+\vec\alpha}{m_1},\qquad
\vec v_2=\frac{\vec p_2-\vec\alpha}{m_2}.
}
$$

Here $\vec\alpha$ must be evaluated with the physical $\dot r$. These equations
are implicit because $\vec\alpha$ itself contains $\dot r$.

The implicit radial relation can be reduced to one scalar equation. Define

$$
\mu=\frac{m_1m_2}{m_1+m_2},\qquad
p_r=
\mu\hat r\mathbin{\cdot}
\left(\frac{\vec p_1}{m_1}-\frac{\vec p_2}{m_2}\right).
$$

Here $\mu$ is the reduced mass and $p_r$ is the canonical momentum conjugate to
the relative radial coordinate.

Then

$$
\boxed{
p_r=
\left(\mu-\frac{q_1q_2}{rc^2}\right)\dot r,
\qquad
\dot r=
\frac{p_r}{\mu-q_1q_2/(rc^2)}.
}
$$

Here the second expression is valid when its denominator is nonzero. It gives
the physical radial velocity directly from the two canonical momenta.

For like charges, the denominator vanishes at

$$
\rho=\frac{q_1q_2}{\mu c^2}.
$$

Here $\rho$ is Weber's critical radius for a like-charge pair. The corrected
implementation must define and test its behavior near this singular relation.

## Correct energy and canonical Hamiltonian

The Legendre transform first gives the conserved energy in physical
velocity variables:

$$
E(\vec r,\vec v)
=
\frac12\sum_i m_i\lVert\vec v_i\rVert^2
+
\sum_{i<j}\frac{q_iq_j}{r_{ij}}
\left(1-\frac{\dot r_{ij}^{\,2}}{2c^2}\right).
$$

Here $E$ is the conserved energy, the first sum is physical kinetic energy,
and the second sum is the Weber pair energy. This formula is correct in the
current theory and paper.

To obtain a canonical Hamiltonian, the physical velocities must first be
obtained from the canonical-momentum equations. The exact result can be
written without expanding that simultaneous solve:

$$
\boxed{
H(\vec r,\vec p)
=
E\bigl(\vec r,\vec v(\vec r,\vec p)\bigr)
=
\frac12\sum_i\vec p_i\mathbin{\cdot}\vec v_i(\vec r,\vec p)
+
\sum_{i<j}\frac{q_iq_j}{r_{ij}}.
}
$$

Here $H$ is the canonical Hamiltonian, $\vec p$ denotes all canonical
momenta, and $\vec v(\vec r,\vec p)$ denotes the physical velocities obtained
by solving the canonical-momentum equations. The last sum is the Coulomb
part of the interaction.

The equality between the two forms follows from the exact Legendre transform.
It must not be approximated by substituting
$\vec v_i=\vec p_i/m_i$.

### Explicit two-particle Hamiltonian

For two particles, define

$$
m_{\mathrm{tot}}=m_1+m_2,\qquad
\vec P=\vec p_1+\vec p_2,\qquad
\vec p_{\mathrm{rel}}
=
\mu\left(\frac{\vec p_1}{m_1}-\frac{\vec p_2}{m_2}\right).
$$

Here $m_{\mathrm{tot}}$ is the total mass, $\vec P$ is total canonical
momentum, $\mu$ is the reduced mass, and $\vec p_{\mathrm{rel}}$ is canonical
relative momentum.

Separate the radial and transverse parts by defining

$$
p_r=\hat r\mathbin{\cdot}\vec p_{\mathrm{rel}},\qquad
\vec p_\perp=\vec p_{\mathrm{rel}}-p_r\hat r.
$$

Here $p_r$ is relative radial canonical momentum and $\vec p_\perp$ is the
relative momentum perpendicular to the separation direction.

The exact two-particle canonical Hamiltonian is

$$
\boxed{
H
=
\frac{\lVert\vec P\rVert^2}{2m_{\mathrm{tot}}}
+
\frac{\lVert\vec p_\perp\rVert^2}{2\mu}
+
\frac{p_r^2}{2\left(\mu-q_1q_2/(rc^2)\right)}
+
\frac{q_1q_2}{r}.
}
$$

Here all quantities are the two-particle quantities defined immediately
above. Only relative radial motion receives the Weber modification; the
center-of-mass and relative transverse terms retain their ordinary forms.

## Correct canonical equations

The first canonical equation says that the coordinate rate is the physical
velocity:

$$
\boxed{
\dot{\vec r}_i=\vec v_i(\vec r,\vec p).
}
$$

Here $\dot{\vec r}_i$ is the time derivative of the position of particle $i$,
and $\vec v_i(\vec r,\vec p)$ is obtained from the simultaneous
canonical-momentum equations.

The pair contributions to the second canonical equation give

$$
\boxed{
\dot{\vec p}_i
=
\sum_{j\ne i}
\frac{q_iq_j}{r_{ij}^2}
\left[
\hat r_{ij}
\left(1+\frac{3\dot r_{ij}^{\,2}}{2c^2}\right)
-
\frac{\dot r_{ij}}{c^2}(\vec v_i-\vec v_j)
\right].
}
$$

Here $\dot{\vec p}_i$ is the canonical momentum rate. Every $\dot r_{ij}$ and
every $\vec v_i-\vec v_j$ on the right-hand side uses the physical velocities
obtained from the canonical momenta.

For the $x$ component of a two-particle system, this is

$$
\boxed{
\dot p_{x_1}
=
\frac{q_1q_2}{r^2}
\left[
\frac{x_1-x_2}{r}
\left(1+\frac{3\dot r^{\,2}}{2c^2}\right)
-
\frac{\dot r(\dot x_1-\dot x_2)}{c^2}
\right],
\qquad
\dot p_{x_2}=-\dot p_{x_1}.
}
$$

Here $p_{x_1}$ and $p_{x_2}$ are the canonical $x$ momenta, and
$\dot x_1$ and $\dot x_2$ are physical $x$ velocities. The analogous
$y$ and $z$ equations have the same structure.

Both $1/c^2$ signs differ from the current canonical momentum-rate equation.

Canonical momentum rate is not itself the mechanical Weber force. Because
canonical momentum contains the velocity-dependent pair corrections,
$m_i\dot{\vec v}_i$ also contains the time derivatives of those corrections.
Combining that identity with the corrected canonical equations reproduces the
Weber force law stated in the theory and paper.

## How the error is confirmed in the repository

### Theory

[`theory/WeberElectrodynamics.md`](WeberElectrodynamics.md) contains the
contradiction directly:

1. its opening momentum section defines each component as
   $p_{x_i}=m_i\dot x_i$;
2. its later Lagrangian contains the velocity-dependent Weber term;
3. its two-particle velocity equations have the signs obtained from the
   correct implicit momentum relation, contradicting the opening definitions;
4. its $H=T+U$ expression is correct only while written in physical
   velocities; and
5. its canonical momentum-rate equation has both Weber correction signs
   reversed.

### Companion paper

[`papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex`](../papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex)
correctly derives the velocity-space energy, then explicitly says that every
velocity is replaced by momentum divided by mass. That sentence is the
primary paper error. All later expanded Weber Hamilton equations inherit it.

### Julia package

[`src/hamiltonian/builders/weber.jl`](../src/hamiltonian/builders/weber.jl)
constructs the current Hamiltonian as ordinary
$\sum_i\lVert\vec p_i\rVert^2/(2m_i)$ plus pair terms whose radial velocities
are computed from $\vec p_i/m_i-\vec p_j/m_j$. Symbolics.jl differentiates
that expression correctly, but the expression being differentiated is not the
Legendre transform of the Weber Lagrangian.

The existing energy-conservation tests therefore show that the integrator
conserves the Hamiltonian it was given. They do not independently establish
that this Hamiltonian represents the stated Weber Lagrangian or Weber force.

## Minimum equation changes in the theory notes

### `theory/WeberElectrodynamics.md`

At minimum, eight existing display blocks require mathematical changes:

| Existing location | Minimum correction |
| --- | --- |
| Six component equations under “Momenta for Particle 1 and Particle 2” | Replace the six statements $p=m\dot x$, $p=m\dot y$, and $p=m\dot z$ with the canonical two-particle relation containing the opposite Weber corrections. The compact $\vec p_1$, $\vec p_2$, and $\vec\alpha$ equations above may replace all six blocks. |
| The $n$-particle equation under “Hamiltonian formulation” | Stop presenting $\sum_i T_i+\sum_{i<j}U_{ij}$ as a canonical function of positions and momenta. Retain it as velocity-space energy and add the exact canonical form $\tfrac12\sum_i\vec p_i\cdot\vec v_i+\sum_{i<j}q_iq_j/r_{ij}$, with the velocities obtained from the canonical-momentum equations. |
| The equation for $\dot p_{x_1}$ under “Equations of motion” | Change $1-3\dot r^2/(2c^2)$ to $1+3\dot r^2/(2c^2)$ and change the final plus sign to a minus sign. |

The following formulas can remain:

- the position, velocity, acceleration, and radial-derivative definitions;
- the Weber pair energy and Weber force;
- the Lagrangian;
- the Euler–Lagrange equation;
- the Legendre-transform identity $H=T+U$ when explicitly described as
  velocity-space energy;
- the existing two-particle $\dot x_1$ and $\dot x_2$ formulas, whose signs are
  already correct, provided their $\dot r$ is obtained from the implicit
  canonical-momentum relation; and
- $\dot p_{x_2}=-\dot p_{x_1}$.

### Other theory files

The minimum authoritative theory follow-up is:

| File | Required correction |
| --- | --- |
| [`theory/InitialConditions.md`](InitialConditions.md) | Replace the general statement that canonical momentum is $m_i\vec v_i$ and replace the naive canonical Hamiltonian. Correct every nonzero-radial construction that assigns $p_r=\mu\dot r$; it must use $p_r=(\mu-q_1q_2/(rc^2))\dot r$. Zero-radial and rigid-rotation constructions remain valid at their initial instant. |
| [`theory/NonZeroRadialVelocityBoundICs.md`](NonZeroRadialVelocityBoundICs.md) | Retain its correct forward canonical-momentum formula, but replace the approximate canonical-Hamiltonian discussion with the exact implicit form used in this note. Remove the claim that the current integrator already performs the correct inverse conversion. |

No file under `_research/` is part of the repair.

## Minimum equation changes in the companion paper

The paper has eight existing displayed equation blocks whose mathematical
content must change:

| Paper equation or location | Minimum correction |
| --- | --- |
| `xdot_two` | Replace the minus correction with $\dot x_1=(p_{x_1}+\alpha_x)/m_1$. |
| `xdot_two_p2` | Replace the plus correction with $\dot x_2=(p_{x_2}-\alpha_x)/m_2$. |
| `pdot_two` | Use the corrected two-particle canonical momentum-rate equation stated above. |
| `hamiltonian` | Retain its right-hand side as velocity-space energy $E_n(\vec r,\vec v)$ rather than a canonical Hamiltonian, and add the exact canonical form $\tfrac12\sum_i\vec p_i\cdot\vec v_i(\vec r,\vec p)+\sum_{i<j}q_iq_j/r_{ij}$. |
| `xdot_expanded` | Change the summed Weber correction from minus to plus and state that all pair radial velocities are obtained from the simultaneous canonical-momentum equations. |
| `pdot_expanded` | Change both $1/c^2$ signs to the corrected signs shown in the general canonical equation above. |
| The displayed velocity correction in the non-separability section | Replace the naive isolated pair term after momentum substitution with the canonical kinetic contribution $\tfrac12\sum_i\vec p_i\cdot\vec v_i(\vec r,\vec p)$. |
| The appendix momentum `align*` block | Replace all six kinetic-momentum definitions with the canonical two-particle momentum relation. |

At least two derivational relations must also be inserted:

1. the canonical momentum obtained from
   $\partial L/\partial\vec v_i$; and
2. the exact canonical Hamiltonian after the physical velocities have been
   obtained from the canonical momenta.

The paper equation `alpha_x` may remain algebraically unchanged, but
$\dot r$ must mean physical radial velocity. Equations `H_def`, `H`, and
`h_box` may remain if their surrounding prose makes clear that they still use
physical velocities. The abstract conservation claims, generic Hamilton
equations `xdot` and `pdot`, phase-space definitions, and generic
symplectic-integrator derivation do not require mathematical changes.

The paper's complexity section must also stop claiming that pair evaluation
alone completes an equation evaluation. Pair geometry still costs
$O(n^2)$, but obtaining all physical velocities requires solving a coupled
linear system. The eventual paper correction is a minor paper revision, so
[`papers/Computational-Weber-Electrodynamics/VERSION`](../papers/Computational-Weber-Electrodynamics/VERSION)
must change from `1.2` to `1.3`, and the tracked PDF must be rebuilt.

## Formula verifier

[`papers/Computational-Weber-Electrodynamics/verify_formulas.py`](../papers/Computational-Weber-Electrodynamics/verify_formulas.py)
currently defines radial velocities in canonical variables by using
$\vec p_i/m_i$. Its checks then prove that the paper and verifier contain the
same naive Hamiltonian.

At minimum, the verifier must replace:

- the global two- and three-particle momentum-space radial-rate definitions;
- Group C.3, which verifies the invalid velocity substitution;
- Groups D and E, which verify the current two-particle and expanded
  Hamilton equations;
- Group H, which verifies the current three-particle equations; and
- the Hamiltonian expressions used by Groups G and I for conservation checks.

The corrected verifier must independently check:

1. the canonical momentum obtained from the Lagrangian;
2. the exact Legendre-transform identity;
3. the two-particle scalar inverse for $\dot r$;
4. both corrected canonical equations;
5. at least one three-particle simultaneous velocity solve; and
6. translation and rotation invariance of the corrected canonical
   Hamiltonian.

SymPy is not provisioned by the repository's current Julia environment or
paper workflow. The repair must supply it explicitly and record a passing
verifier run.

## Minimum Julia source changes

The directly incorrect production paths are:

| File | Minimum source change |
| --- | --- |
| [`src/hamiltonian/builders/weber.jl`](../src/hamiltonian/builders/weber.jl) | Replace the naive symbolic Hamiltonian and its $\vec p/m$ radial rates with the exact canonical Hamiltonian. Its pair decomposition must use physical velocities obtained from the canonical momenta. |
| [`src/hamiltonian_system.jl`](../src/hamiltonian_system.jl) | Route the default Weber system through corrected Hamiltonian and equation evaluators. Preserve the generic custom-Hamiltonian constructor. If exact symbolic elimination is not practical for general $n$, the default Weber constructor needs a specialized numerical path while retaining the public compiled-function signatures. |
| [`src/statistics/energy.jl`](../src/statistics/energy.jl) | Compute physical kinetic energy from physical velocities, compute pair radial velocities from those velocities, and compare their decomposition with the corrected compiled Hamiltonian. |
| [`src/statistics/forces.jl`](../src/statistics/forces.jl) | Stop constructing velocities and accelerations from $\vec p/m$; obtain physical velocities from the corrected coordinate-rate equation. |
| [`src/initial_conditions.jl`](../src/initial_conditions.jl) | Convert the `radial_velocity` keyword to canonical radial momentum using $p_r=(\mu-q_1q_2/(rc^2))\dot r$. This may require adding `c` to the helper input. |
| [`ext/WeberElectrodynamicsMakieExt.jl`](../ext/WeberElectrodynamicsMakieExt.jl) | Replace the live pair radial-rate calculation based on $\vec p/m$ with the corrected physical velocities. |
| [`ext/WeberElectrodynamicsJLD2Ext.jl`](../ext/WeberElectrodynamicsJLD2Ext.jl) | Bump the solution archive format or otherwise prevent old default-Weber trajectories from being silently reconstructed as corrected dynamics. |

A shared internal routine should perform the canonical-momentum-to-velocity
solve. Repeating separate implementations in the Hamiltonian, statistics, and
Makie paths would risk another inconsistency. Performance work may also require
workspace buffers in `src/types.jl` and a richer named-term hook in
`src/hamiltonian/terms.jl`, but those are implementation choices rather than
additional physics corrections.

The following source areas do not presently contain the primary error:

- the generic Hamiltonian constructor and generic Symbolics differentiation;
- `kinetic_term` and `coulomb_term` for custom non-Weber systems;
- the symmetric-projection integrator;
- canonical linear and angular momentum statistics;
- trajectory extraction;
- the canonical center-of-mass and relative-momentum transformations in
  `src/regularization.jl` and `src/solve.jl`; and
- the generic collision callback transformations.

Regularization and collision behavior must nevertheless be revalidated because
the corrected Weber Hamiltonian changes the dynamics they receive.

## Minimum tests and generated artifacts

At minimum, the following tests contain expectations tied directly to the old
Hamiltonian or to $\vec p/m$ as physical velocity:

- [`test/test_utils.jl`](../test/test_utils.jl);
- [`test/test_hamiltonian_system.jl`](../test/test_hamiltonian_system.jl);
- [`test/test_builders.jl`](../test/test_builders.jl);
- [`test/test_custom_hamiltonian.jl`](../test/test_custom_hamiltonian.jl);
- [`test/test_named_term.jl`](../test/test_named_term.jl);
- [`test/test_statistics.jl`](../test/test_statistics.jl);
- [`test/test_initial_conditions.jl`](../test/test_initial_conditions.jl); and
- [`test/test_physics.jl`](../test/test_physics.jl).

The corrected suite must add independent checks that do not calculate expected
values from the implementation under test:

1. compare canonical momentum with a direct derivative of the Lagrangian;
2. compare the corrected Hamiltonian with the Legendre transform;
3. finite-difference the corrected Hamiltonian with respect to positions and
   momenta;
4. recover the mechanical Weber force from the canonical equations;
5. check the Coulomb and zero-radial-velocity limits;
6. test nonzero radial motion in two- and three-particle systems; and
7. test behavior near the critical denominator.

`test/test_integration.jl`, `test/test_regularization.jl`, and
`test/test_callbacks.jl` must be rerun and updated where trajectory behavior
changes. All seven committed files under `test/regression/fixtures/` represent
the old default Weber dynamics and must be regenerated after independent
validation passes.

The default-Weber portions of these notebooks must be rerun:

- [`examples/two_body_reference.ipynb`](../examples/two_body_reference.ipynb);
- [`examples/api_showcase.ipynb`](../examples/api_showcase.ipynb); and
- [`examples/unlike_collision_bounce.ipynb`](../examples/unlike_collision_bounce.ipynb).

The three tracked PNGs under `examples/figures/` must also be regenerated.
Pure custom-Coulomb notebook sections are unaffected.

## Minimum documentation changes

After the source API is settled, update:

- [`docs/src/api/system.md`](../docs/src/api/system.md) for the corrected
  default Weber builder and named-term behavior;
- [`docs/src/api/statistics.md`](../docs/src/api/statistics.md) for physical
  kinetic energy and radial-velocity semantics;
- [`docs/src/internals.md`](../docs/src/internals.md) if the default Weber
  system no longer follows the current all-symbolic construction path; and
- [`docs/src/quickstart.md`](../docs/src/quickstart.md) if the
  initial-condition helper signature or displayed momentum conventions change.

Generated `docs/build/` output is not an authoritative input and should be
rebuilt through Documenter rather than edited by hand.

## Work estimate

This is a medium-to-large correction, not a sign-only patch. A realistic
minimum estimate for one contributor already familiar with the package is:

| Work area | Estimate |
| --- | --- |
| Final derivation, paper source, and formula verifier | 1–2 focused days |
| Correct default Hamiltonian and shared velocity solve | 3–5 focused days |
| Statistics, initial conditions, Makie, and archive handling | 2–3 focused days |
| Independent tests and regularization/collision revalidation | 2–4 focused days |
| Documentation, notebooks, fixtures, figures, and final full validation | 1–3 focused days |

The total is approximately 9–17 focused working days. The largest uncertainty
is whether the current symbolic representation remains practical for general
particle counts or whether the default Weber system needs specialized runtime
equation evaluators.

## Minimum implementation order

1. Correct `WeberElectrodynamics.md`, the paper derivation, and the independent
   formula verifier.
2. Implement one tested canonical-momentum-to-velocity solve.
3. Implement the corrected canonical Hamiltonian and both canonical equations.
4. Update energy, force, initial-condition, and Makie consumers.
5. Decide archive compatibility for trajectories generated by the old
   dynamics.
6. Add independent physics tests before replacing regression fixtures.
7. Revalidate regularization and collision behavior.
8. Update documentation and regenerate notebooks, figures, fixtures, the
   paper PDF, and Documenter output.

No partial patch that changes only the two signs in $\dot{\vec p}$ should be
merged. Those signs are consequences of the same canonical-momentum error and
cannot be corrected independently.

## Acceptance criteria

The finding is resolved only when:

1. theory, paper, verifier, and package all use
   $\vec p_i=\partial L/\partial\vec v_i$;
2. no Weber path obtains physical velocity by setting
   $\vec v_i=\vec p_i/m_i$;
3. the canonical Hamiltonian is the exact Legendre transform of the stated
   Weber Lagrangian;
4. the coordinate rate equals the recovered physical velocity;
5. the canonical momentum rate has the corrected signs;
6. the canonical equations reproduce the mechanical Weber force;
7. nonzero-radial initial conditions are converted to canonical momenta;
8. old trajectories cannot be mistaken for corrected Weber results;
9. independent two- and three-particle physics checks pass;
10. regularization, collision, regression, package, documentation, paper, and
    verifier checks all pass; and
11. affected notebooks, figures, fixtures, and the paper PDF have been
    regenerated from the corrected system.

Only after these criteria are satisfied should this planning note be removed.
