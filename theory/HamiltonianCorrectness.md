# Weber Hamiltonian Correctness Finding and Minimum Repair

> **Status (2026-08-20, branch `fix/canonical-weber-hamiltonian`)**
>
> | Area | State |
> | --- | --- |
> | Physics finding below | Confirmed, independently re-derived |
> | Companion paper | **Mostly corrected** — six sign fixes and the momentum definition landed; **Tier 0 below is still outstanding** |
> | Formula verifier | **Done** — 94 symbolic checks, all passing |
> | Theory notes (`theory/*.md`) | **Not started** — Tier 1 and Tier 2 below |
> | Julia package, tests, fixtures, notebooks | **Not started** — Phase 2 below, deliberately deferred |
>
> This note is the work order. Read Tier 0 first: everything else aligns *to*
> the paper, so the paper must be finished before the theory notes are
> rewritten to match it.
>
> The `_research/` directory is outside the scope of this repair.

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

## Tier 0 — remaining companion-paper items

Five commits on this branch (`5cc03c5`, `129c2da`, `bdd93f9`, `177ffc0`,
`f9b6eb7`) applied most of this note's original paper table; `129c2da`
narrowed `5cc03c5`, so the table below reflects the **verified current
state of the file**, not the sum of the commits. What landed, and what did not, is recorded here so the next agent does
not redo finished work or assume unfinished work is done.

All line numbers refer to
[`papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex`](../papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex)
at commit `87135d7`.

### Already applied — do not redo

| Item | Result |
| --- | --- |
| `xdot_two` | `ẋ₁ = (p_{x₁} + α_x)/m₁` ✔ |
| `xdot_two_p2` | `ẋ₂ = (p_{x₂} − α_x)/m₂` ✔ |
| `pdot_two`, `pdot_expanded` | both `1/c²` signs corrected ✔ |
| `xdot_expanded` | summed Weber correction now `+` ✔ |
| Canonical momentum `p_{x₁} = ∂L/∂ẋ₁` | inserted as `px1` ✔ |
| Appendix `p = m v` block | deleted ✔ |
| "replace each velocity by `p_{x_i}/m_i`" | deleted ✔ |
| Critical-radius wording | "At this critical radius" → "Below this critical radius" ✔ |

The critical-radius fix was not in this note's original table. It is recorded
because the reversal happens on *crossing* ρ: at `r = ρ` the effective inertia
`μ_eff = μ − q₁q₂/(c²r)` is exactly zero and `r̈` is singular, not merely
sign-reversed. Verifier checks `L.1`, `L.5`–`L.8` pin this down, including the
proof that the reversal holds at arbitrary angular momentum, not only for
radial motion.

### Required — the minimum that resolves the finding

Four items. Nothing here is optional; each one is either wrong or would lead a
reader to implement the wrong integrator.

| Line | Item | Required change |
| --- | --- | --- |
| 169 | "The *Legendre transform* then yields" | Followed only by an argument list; the transform is never performed. Insert the exact canonical Hamiltonian — for `n = 2` the closed form in "Explicit two-particle Hamiltonian" below is exact and also exhibits ρ as the pole. **Doing this one item largely resolves lines 211 and 258 as well**, which is why it is listed first. |
| 211 | `hamiltonian` | `H_n` is written in velocities (`v_i²`, `ṙ_ij`) but is immediately fed to `ẋ_i = ∂H_n/∂p_{x_i}` at line 220. Either relabel it as velocity-space energy `E_n(r,v)`, or state that the velocities are functions of `(q,p)`. Last surviving instance of the original error's presentation. |
| 258 | "are computed using `xdot_expanded` … at time `t`" | Circular as written: the right-hand side needs `ṙ_ij` at time `t`, built from the velocities at time `t` — the quantity the equation produces. State that velocities are recovered from the momenta each step. For two particles this is one scalar: `ṙ = p_r/(μ − q₁q₂/(rc²))`; no matrix notation is needed. |
| 306, 308 | complexity section | `n(n−1)/2` is right for assembling pair geometry, but "each evaluation … only needs to compute `r_ij` and `ṙ_ij`" hides the cost: obtaining `ṙ_ij` from `(q,p)` *is* the velocity solve. Name it as a separate per-step cost. |

### Optional — skip unless polishing

Not required to resolve the finding. Recorded so they are not re-discovered as
new defects; leave them alone if the goal is the minimum repair.

- Line 99: "It satisfies Newton's third law … **Therefore**, energy, linear
  momentum, and angular momentum are always conserved." The strong third law
  gives linear and angular momentum; energy conservation is separate, from `L`
  having no explicit time dependence. Over-claims via "therefore".
- Line 204: "`ṗ_{x₁} + ṗ_{x₂} = 0`, **but** `ẋ₁ − ẋ₂ ≠ 0` in general" is a
  non-sequitur — a sum contrasted with an unrelated difference. The real point
  is that `α_x` cancels in the sum but not in either momentum alone.
- Line 326: the non-separability display shows an isolated pair term; after the
  correct transform the non-separable object is the whole `½ Σ p_i · v_i(q,p)`.
- Line 355: "for `(q,p)` and `(x,y)` separately" contradicts the bullets below,
  which correctly use the mixed pairs `(x,p)` and `(q,y)`. The bullets are right.
- Lines 251–258: the Euler steps are legitimate — inside `Φ^A`/`Φ^B` the
  right-hand side is frozen, so Euler is the *exact* sub-flow. One sentence
  would say so.
- Line 165: `px1` uses bare `ṙ`, `r²` where neighbours use `ṙ₁₂`, `r₁₂²`;
  `r₁₂` is never defined, only `r` is.
- Line 376: `Z_n` uses `n` as a time-step index while `n` is the particle count.

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

### Companion paper (historical — now corrected)

v1.2 of
[`Computational-Weber-Electrodynamics.tex`](../papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex)
derived the velocity-space energy correctly and then said that every velocity
is replaced by momentum divided by mass. That sentence was the primary paper
error; all later expanded Weber Hamilton equations inherited it.

`git show main:` confirms v1.2 carried both halves of the defect together —
the `p/m` substitution **and** `ẋ₁ = (p_{x₁} − α_x)/m₁`. That is not a
coincidence: taking `∂H/∂p` of the naively-substituted Hamiltonian *produces*
the minus sign. The v1.2 paper and the package were therefore a consistent,
mutually-reinforcing pair — which is exactly why the error survived review.

Both halves are fixed as of commit `5cc03c5`. Tier 0 lists what remains.

### Julia package (still wrong)

[`src/hamiltonian/builders/weber.jl`](../src/hamiltonian/builders/weber.jl)
constructs the Hamiltonian as ordinary `Σ ‖p_i‖²/(2m_i)` plus pair terms whose
radial velocities are computed at line 60 from `p_i/m_i − p_j/m_j`.
Symbolics.jl differentiates that expression correctly, but the expression being
differentiated is not the Legendre transform of the Weber Lagrangian.

`git diff main..HEAD -- src/` is empty: the package still implements v1.2.

Two consequences are worth stating plainly, because they explain why the
existing test suite cannot catch this:

1. The package's Hamiltonian is a **perfectly good** Hamiltonian — autonomous,
   translation-invariant, rotation-invariant. A symplectic integrator conserves
   it, so **energy, linear momentum and angular momentum all come out flat**.
   The three headline conservation checks are structurally incapable of
   detecting the defect.
2. The discrepancy is `O(ρ/r)` in the momentum–velocity relation, so far from
   the critical radius it is numerically small. Near ρ it is not: the correct
   `H` has a pole at `μ_eff = 0`, while the package's `H` has `μ_eff` in the
   **numerator** and is finite there. The critical-radius physics — the
   paper's headline future direction — is simply absent from the equations the
   package integrates.

Energy-conservation tests therefore show that the integrator conserves the
Hamiltonian it was given. They do not establish that this Hamiltonian
represents the stated Weber Lagrangian or Weber force.

## Tier 1 — theory notes that carry wrong mathematics

These are the documents a reader would trust and copy from. Fix
`WeberElectrodynamics.md` first: it is the root note that both the paper and
the package descend from.

### `theory/WeberElectrodynamics.md`

| Location | Required correction |
| --- | --- |
| Lines 82–102, six blocks under "Momenta for Particle 1 and Particle 2" | Replace `p_{x₁} = m₁ẋ₁`, `p_{y₁} = m₁ẏ₁`, … with the canonical relation. The compact `p⃗₁ = m₁v⃗₁ − α⃗`, `p⃗₂ = m₂v⃗₂ + α⃗` pair from this note's "Two-particle form" replaces all six blocks. |
| Line 342, `ṗ_{x₁}` under "Equations of motion" | `1 − 3ṙ²/(2c²)` → `1 + 3ṙ²/(2c²)`, and the trailing `+ ṙ(ẋ₁−ẋ₂)/c²` → `−`. Both signs together; neither is independently correctable. |
| The `H = T + U` block under "Hamiltonian formulation" | Keep it, but label it explicitly as **velocity-space energy**. Add the exact canonical form `½ Σ p⃗_i · v⃗_i(r,p) + Σ_{i<j} q_i q_j / r_{ij}`, with velocities obtained from the canonical-momentum equations. |

Insert two derivational steps the note currently skips: `p⃗_i = ∂L/∂v⃗_i`, and
the two-particle scalar inverse `ṙ = p_r/(μ − q₁q₂/(rc²))`.

These formulas may remain unchanged: the position, velocity, acceleration and
radial-derivative definitions; the Weber pair energy and Weber force; the
Lagrangian; the Euler–Lagrange equation; the two-particle `ẋ₁` and `ẋ₂`
formulas, whose signs are already correct provided their `ṙ` is the physical
one; and `ṗ_{x₂} = −ṗ_{x₁}`.

### `theory/InitialConditions.md`

| Location | Required correction |
| --- | --- |
| Line 14 | `p⃗_i = m_i v⃗_i — canonical momentum` is the defect stated verbatim. Replace with the canonical definition. |
| Naive canonical Hamiltonian | Replace with the exact form from this note. |
| Every nonzero-radial construction assigning `p_r = μṙ` | Must become `p_r = (μ − q₁q₂/(rc²))ṙ`. |

State explicitly which constructions survive untouched: zero-radial and
rigid-rotation constructions impose `ṙ_ij(0) = 0`, so `α⃗ = 0` and canonical
equals kinetic momentum **at the initial instant**. That is the majority of
the file, and saying so prevents an over-broad rewrite.

## Tier 2 — documents whose status claims are stale

### `theory/NonZeroRadialVelocityBoundICs.md`

Line 21 is already correct — it distinguishes canonical from kinetic momentum.
Retain the forward canonical-momentum formula. Replace the approximate
canonical-Hamiltonian discussion with the exact implicit form used in this
note, and remove any claim that the current integrator performs the correct
inverse conversion. It does not.

### This note itself

After Tiers 0–1 land, update the status table at the top, and mark the Tier 0
"already applied" rows as they are completed. Do **not** delete this note when
the paper and theory notes are done: its acceptance criteria are the only
record that the package is still wrong. Retire it only when every criterion
passes.

## Tier 3 — housekeeping

- Bump [`VERSION`](../papers/Computational-Weber-Electrodynamics/VERSION) from `1.2` to `1.3`.
- Delete the review-markup block at `.tex` lines 54–65 once Tier 0 is accepted.
- Rebuild the paper PDF locally. It is **not** tracked: `.gitignore:97`
  matches `papers/**/*.pdf`, and commit `a111a74` removed it from the
  index. Any instruction to rebuild "the tracked PDF" is stale.
- Add a `CHANGELOG.md` entry describing the corrected Hamiltonian.

## Documents that need no change

- [`docs/src/theory.md`](../docs/src/theory.md) — twelve lines of links, no
  equations. It inherits correctness from Tier 1 automatically.
- [`docs/src/api/*.md`](../docs/src/api), [`docs/src/internals.md`](../docs/src/internals.md),
  [`docs/src/quickstart.md`](../docs/src/quickstart.md) — these describe the
  **API surface**, which cannot be written until the source settles. They
  belong to Phase 2, not to the documentation tiers.
- `docs/build/` — generated by Documenter; never hand-edited.
- `_research/**` — explicitly out of scope.

## Guard rails — code that looks wrong but is not

A `p_i/m_i` pattern is not by itself the defect. The next agent will grep for
it and must not "fix" these:

- [`src/regularization.jl`](../src/regularization.jl) lines 236–238 compute
  `rel_p[d] = μ (p_i/m_i − p_j/m_j)`. That is the **exact canonical relative
  momentum** `(m_j p⃗_i − m_i p⃗_j)/(m_i+m_j)`, a linear map on canonical
  momenta that is valid whether or not `p = mv`. Only the local variable
  names `vi_d`, `vj_d` are misleading; the mathematics is correct. Rename the
  variables, change nothing else.
- The same canonical centre-of-mass and relative transformations in
  [`src/solve.jl`](../src/solve.jl) are likewise correct.

The genuine defect is narrower: computing a **radial velocity** `ṙ_ij` from
`p_i/m_i − p_j/m_j`, or a **kinetic energy** from `Σ p_i²/(2m_i)`. Confirmed
occurrences are listed in the Phase 2 table below.

## Formula verifier — done

[`papers/Computational-Weber-Electrodynamics/verify_formulas.py`](../papers/Computational-Weber-Electrodynamics/verify_formulas.py)
has been rebuilt and is **complete**: 94 checks, all passing, re-based on the
single source of truth `L = T − S`. Nothing the paper asserts is assumed.

```
uv run --with sympy python3 verify_formulas.py     # 94/94, ~35 s
```

SymPy is not provisioned by the repository's Julia environment; `uv run --with
sympy` is the recorded invocation.

Groups, and what each settles:

| Group | Covers |
| --- | --- |
| A–C | potential identities; Euler–Lagrange → Weber force, in both relative and the paper's absolute coordinates; `H = T + U` |
| D–E | two-particle Hamilton equations; the `n`-particle equations at `n = 2` |
| F–G | Appendix A.1 radial-acceleration identities; conservation laws |
| H–I | `n = 3` equations of motion; rotational invariance |
| J | `dU/dt = −F ṙ`; Newton's third law in the strong form; centrality |
| K | the **closed-form** `H(q,p)` for `n = 2`, and the end-to-end chain Hamilton's equations → Weber's force law |
| L | critical radius: `μ_eff(ρ) = 0`, divergence at ρ, sign reversal below ρ **at arbitrary angular momentum** |
| M | extended phase space: `A Aᵀ = 2I`, the `1/4` Newton factor, the `Φ^A`/`Φ^B` flow structure, second order, symplecticity, constraint drift |
| N | the `n`-body momentum–velocity relation, and the `ρ/r` contraction factor of the naive fixed-point solve |
| O | units, dimensions, pair count |
| P | **negative controls** — each corrected sign is shown to be *necessary*, and the v1.2 formulation is shown not to reproduce Weber's force law |

Two results are load-bearing for the rest of this note:

- **`P.9` / `P.10`.** Feeding both Hamiltonians through Hamilton's equations
  and testing against Weber's force law gives residual **exactly 0** for the
  corrected `H`, and **non-zero** for the v1.2 / `weber.jl` form. This is the
  decisive test; it does not depend on any interpretation argument.
- **`N.7` / `N.8`.** Seeding the implicit velocity relation with `p/m` and
  iterating is a Jacobi solve whose contraction factor is exactly `ρ/r`. It
  converges only *above* the critical radius and diverges below it — the
  regime the paper calls its most promising direction. The velocity recovery
  must therefore be a real solve, not an iteration.

The suite is mutation-tested: flipping `L = T − S` to `L = T + S` fails 12
checks. Any Phase 2 work should extend this file rather than start a new one.

## Phase 2 — Julia source (deferred; not part of the documentation pass)

The directly incorrect production paths are:

| File | Minimum source change |
| --- | --- |
| [`src/hamiltonian/builders/weber.jl`](../src/hamiltonian/builders/weber.jl) (confirmed: line 60) | Replace the naive symbolic Hamiltonian and its $\vec p/m$ radial rates with the exact canonical Hamiltonian. Its pair decomposition must use physical velocities obtained from the canonical momenta. |
| [`src/hamiltonian_system.jl`](../src/hamiltonian_system.jl) | Route the default Weber system through corrected Hamiltonian and equation evaluators. Preserve the generic custom-Hamiltonian constructor. If exact symbolic elimination is not practical for general $n$, the default Weber constructor needs a specialized numerical path while retaining the public compiled-function signatures. |
| [`src/statistics/energy.jl`](../src/statistics/energy.jl) (confirmed: lines 128–129) | Compute physical kinetic energy from physical velocities, compute pair radial velocities from those velocities, and compare their decomposition with the corrected compiled Hamiltonian. |
| [`src/statistics/forces.jl`](../src/statistics/forces.jl) (confirmed: lines 206–207) | Stop constructing velocities and accelerations from $\vec p/m$; obtain physical velocities from the corrected coordinate-rate equation. |
| [`src/initial_conditions.jl`](../src/initial_conditions.jl) (confirmed: line 109, `rel_p[1] = reduced_mass * radial_velocity`) | Convert the `radial_velocity` keyword to canonical radial momentum using $p_r=(\mu-q_1q_2/(rc^2))\dot r$. This may require adding `c` to the helper input. |
| [`ext/WeberElectrodynamicsMakieExt.jl`](../ext/WeberElectrodynamicsMakieExt.jl) (confirmed: line 125) | Replace the live pair radial-rate calculation based on $\vec p/m$ with the corrected physical velocities. |
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

### Phase 2 — tests and generated artifacts

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

### Phase 2 — API documentation

Only after the source API settles:
[`docs/src/api/system.md`](../docs/src/api/system.md),
[`docs/src/api/statistics.md`](../docs/src/api/statistics.md),
[`docs/src/internals.md`](../docs/src/internals.md), and
[`docs/src/quickstart.md`](../docs/src/quickstart.md). Generated `docs/build/`
output is rebuilt through Documenter, never hand-edited.

## Work estimate

| Work area | Estimate | State |
| --- | --- | --- |
| Derivation, paper source, formula verifier | 1–2 days | **mostly done**; Tier 0 required items remain (~half a day) |
| Theory notes (Tiers 1–3) | ~1 day | not started |
| Correct default Hamiltonian and shared velocity solve | 3–5 days | not started |
| Statistics, initial conditions, Makie, archive handling | 2–3 days | not started |
| Independent tests, regularization/collision revalidation | 2–4 days | not started |
| Notebooks, fixtures, figures, Documenter, final validation | 1–3 days | not started |

Remaining total is roughly 9–16 focused days, almost all of it Phase 2. The
largest uncertainty is unchanged: whether the current symbolic representation
stays practical for general particle counts, or whether the default Weber
system needs specialised runtime equation evaluators.

## Implementation order

**Documentation pass (this branch, now):**

1. Tier 0 — finish the companion paper. Everything downstream aligns to it.
2. Tier 1 — `theory/WeberElectrodynamics.md`, then `theory/InitialConditions.md`.
3. Tier 2 — `theory/NonZeroRadialVelocityBoundICs.md`, then the status table
   at the top of this note.
4. Tier 3 — `VERSION`, review-markup removal, PDF rebuild, `CHANGELOG.md`.

Gate before moving on: every boxed equation in `theory/WeberElectrodynamics.md`
should name the `verify_formulas.py` check that proves it. That cross-reference
is the mechanism that stops the theory note and the verifier drifting apart
again — undetected drift between the paper and the package is precisely how
v1.2 and the current `weber.jl` became a mutually-consistent wrong pair.

**Package pass (Phase 2, deliberately deferred):**

5. Implement one tested canonical-momentum-to-velocity solve, shared by every
   consumer. Repeating it in the Hamiltonian, statistics and Makie paths would
   risk another divergence.
6. Implement the corrected canonical Hamiltonian and both canonical equations.
7. Update energy, force, initial-condition and Makie consumers.
8. Decide archive compatibility for trajectories generated by the old dynamics.
9. Add independent physics tests **before** replacing regression fixtures.
10. Revalidate regularization and collision behaviour.
11. Regenerate notebooks, figures, fixtures and Documenter output.

No partial patch that changes only the two signs in `ṗ⃗` should be merged.
Those signs are consequences of the same canonical-momentum error and cannot
be corrected independently.

## Acceptance criteria

**Documentation pass (this branch) is done when:**

1. the four Tier 0 required items are applied and the paper builds;
2. `theory/WeberElectrodynamics.md` and `theory/InitialConditions.md` use
   `p⃗_i = ∂L/∂v⃗_i` and contain no `p = m v` momentum definition;
3. no theory note obtains physical velocity by setting `v⃗_i = p⃗_i/m_i`;
4. `theory/NonZeroRadialVelocityBoundICs.md` no longer claims the current
   integrator performs the inverse conversion;
5. `VERSION` reads `1.3` and the review markup is gone; and
6. `verify_formulas.py` still passes 94/94.

**The finding as a whole is resolved only when, in addition:**

7. the canonical Hamiltonian in the package is the exact Legendre transform of
   the stated Weber Lagrangian;
8. the coordinate rate equals the recovered physical velocity, and the
   canonical momentum rate has the corrected signs;
9. the canonical equations reproduce the mechanical Weber force;
10. nonzero-radial initial conditions are converted to canonical momenta;
11. old trajectories cannot be mistaken for corrected Weber results;
12. independent two- and three-particle physics checks pass, including
    behaviour near `μ_eff = 0`; and
13. regularization, collision, regression, package, documentation and paper
    checks all pass, with notebooks, figures and fixtures regenerated.

Retire this note only when criteria 1–13 are satisfied. Until then it is the
only record that the package is still integrating the v1.2 equations.
