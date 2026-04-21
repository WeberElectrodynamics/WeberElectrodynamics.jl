# Agent 13 — Zöllner-enabled comparison

**Scope.** Re-run five ICs flagged by Agents 5, 6, and 8 with
`ZollnerOptions(enabled=true, a=a_val)` for `a ∈ {0.0, 0.1, 0.5, 1.0}`,
and quantify whether strengthening unlike-sign pair attraction
(κ_ij = 1+a for `q_i·q_j < 0`, κ_ij = 1 otherwise) stabilises the
4-body 2+/2− Weber problem or enlarges bound-state basins. Budget:
5 × 4 = 20 integrations, `tmax = 30`, `dt = 1e-3`, `c = 1`, no
regularization bounce. Script: `zollner_comparison.jl`, data:
`zollner_comparison.csv`.

## 1. IC inventory

| id | source | configuration |
|----|--------|---------------|
| `F1b_outbreath_s2.00_v0.50` | Agent 5, §2 (only closed orbit) | alternating square side 2, outward radial `v = 0.5` |
| `rhombus_a1.5_b1.15_eta0.75` | Agent 6, §3 (only Success) | rhombus `(a,b)=(1.5,1.15)`, η=0.75, rotating |
| `alt_square_s1.0_eta0.25` | Agents 6/7 chaotic baseline | alternating square side 1, η=0.25, rotating |
| `two_dimers_a0.2_R2.0_intra0.3` | Agent 5 F3 / Agent 6 dimer | two (+−) dyads, dyad 0.2, separation 2.0 |
| `nucleus_plus_orbiters_A` | Agent 8 Experiment A | (+,+) at ±0.25, p=±0.2; (−,−) at (±2.5, 5), p=±0.1 |

## 2. Results table

| IC | a | retcode | t_final | ΔE % | bound | r_max |
|----|--:|---------|--------:|-----:|:-----:|-----:|
| F1b_outbreath_s2.00_v0.50 | 0.0 | Failure | 12.91 | 0.181 | yes | 5.66 |
| F1b_outbreath_s2.00_v0.50 | 0.1 | Failure | 10.41 | 0.220 | yes | 5.30 |
| F1b_outbreath_s2.00_v0.50 | 0.5 | Failure |  5.91 | 0.062 | yes | 4.64 |
| F1b_outbreath_s2.00_v0.50 | 1.0 | Failure |  3.85 | 0.036 | yes | 4.34 |
| rhombus_a1.5_b1.15_eta0.75 | 0.0 | **Success** | 30.00 | 0.076 | no  | 29.6 |
| rhombus_a1.5_b1.15_eta0.75 | 0.1 | **Success** | 30.00 | 1.65  | no  | 28.6 |
| rhombus_a1.5_b1.15_eta0.75 | 0.5 | Failure | 1.21 | 9.85 | yes | 3.0 |
| rhombus_a1.5_b1.15_eta0.75 | 1.0 | Failure | 0.86 | 6.85 | yes | 3.0 |
| alt_square_s1.0_eta0.25 | 0.0 | Failure | 2.80 |  34.3 | yes | 1.42 |
| alt_square_s1.0_eta0.25 | 0.1 | Failure | 1.95 | 283   | yes | 1.44 |
| alt_square_s1.0_eta0.25 | 0.5 | Failure | 0.83 |  19.4 | yes | 1.41 |
| alt_square_s1.0_eta0.25 | 1.0 | Failure | 0.32 |   0.62 | yes | 1.41 |
| two_dimers_a0.2_R2.0_intra0.3 | 0.0 | Failure | 0.017 | 0.073 | yes | 2.01 |
| two_dimers_a0.2_R2.0_intra0.3 | 0.1 | Failure | 0.015 | 0.053 | yes | 2.01 |
| two_dimers_a0.2_R2.0_intra0.3 | 0.5 | Failure | 0.012 | 0.466 | yes | 2.01 |
| two_dimers_a0.2_R2.0_intra0.3 | 1.0 | Failure | 0.009 | 0.571 | yes | 2.01 |
| nucleus_plus_orbiters_A | 0.0 | Failure | 0.207 | 0.365 | yes | 5.71 |
| nucleus_plus_orbiters_A | 0.1 | Failure | 0.207 | 0.296 | yes | 5.71 |
| nucleus_plus_orbiters_A | 0.5 | Failure | 0.207 | 0.172 | yes | 5.71 |
| nucleus_plus_orbiters_A | 1.0 | Failure | 0.207 | 0.138 | yes | 5.71 |

"bound" = max pairwise distance at any snapshot stayed below escape
radius 20.0 up to whenever the run stopped. Every `:Failure` row is
still inside its initial envelope — boundedness at the failure time
says nothing about long-term confinement.

## 3. Answers to the four questions

**(a) Does `a` stabilise F1b, the breathing square periodic orbit?**
**No — it actively destabilises it.** At `a=0` the integrator survives
to `t = 12.9` (~1.1 periods of the Agent-5 orbit `T = 11.78`). At
`a = 0.1, 0.5, 1.0` the run fails earlier: `t = 10.4, 5.9, 3.9`. This
is consistent with the orbit being a `D₄ × T` brake-symmetric mode whose
shape is tuned to the bare Coulomb pair ratios; boosting κ on the four
unlike edges (but not the two like diagonals) breaks the exact force
balance the orbit relies on. The orbit's monodromy spectral radius
(228.6 per period) means even small mistuning shifts the trajectory off
the closed brake return and lets it drive into a near-collision. Local
energy drift is small (≤0.25 %) in every case — the failures are
geometric collisions, not energy blow-ups.

**(b) Does `a` enlarge the bound-state basin around the rhombus?**
**No — it catastrophically shrinks it.** `a = 0` and `a = 0.1` both
reach `tmax = 30` with `:Success`, but at `a = 0.5` the run already
fails at `t = 1.21` and at `a = 1.0` at `t = 0.86`, both with energy
drifts of 7–10 %. The rhombus Success at η=0.75 sits on an escape
trajectory (`r_max ≈ 29.6` in 30 t-units: the two pairs drift apart
quasi-Kepler). Turning Zöllner attraction up pulls the escaping
unlike pairs back together and straight into a collision before one
transit completes. So Zöllner converts an unbounded-but-clean run into
a bounded-but-singular run. This is a real effect of stronger binding,
not numerical — the energy-drift jump from 0.08 % (a=0) to 9.85 %
(a=0.5) reports the integrator hitting a near-singular event.

**(c) Does `a` help the (+,+) nucleus + (−,−) orbiters?**
**Not in the failure time, but in the local error.** All four `a`
values fail at exactly the same time (`t = 0.207`), which means the
projection step is hitting the same near-singular event driven by the
pre-existing (+,+) angular-momentum leak that Agent 8 diagnosed. Energy
drift at failure drops monotonically with `a` (0.365 → 0.296 → 0.172
→ 0.138 %), consistent with the theoretical note: Zöllner cannot
touch the (+,+) internal spiral (κ₁₂ = 1 for like pairs), so it cannot
regularise the topological obstruction. Its indirect effect is to
deepen the ambient well, which slows the leak of ℓ into pair (1,2) and
reduces the stress on the integrator, but this is quantitative
smoothing — the singular spiral still wins at the same clock time
because the dominant failure path is geometric (Frauenfelder–Weber
Thm 2.1), not energetic.

**(d) Does `a` change chaotic diffusion rates?**
Partially. Looking at the chaotic `alt_square_s1.0_eta0.25` baseline,
`t_final` falls monotonically with `a` (2.80 → 1.95 → 0.83 → 0.32).
The local energy drift signature is erratic (34 → 283 → 19 → 0.62 %),
which is the fingerprint of trajectories hitting different kinds of
close encounters. Stronger unlike attraction means *faster* close
encounters, i.e. an *accelerated* effective diffusion toward the
collision boundary rather than a calmer torus. The η=0.25 square is
far inside the chaotic zone, and Zöllner is pushing it deeper.

## 4. Theoretical reading

Zöllner's `a` modifies only the four unlike-sign κ_ij. For the 2+/2−
problem this means it multiplies the four edge attractions (1–3,
1–4, 2–3, 2–4) while leaving the two diagonal like-pair repulsions
(1–2, 3–4) untouched. Every tested IC has a closer-to-collision
unlike-pair approach than like-pair, so increasing `a` deepens exactly
the wells whose bottom is topologically non-regularisable (ℓ≠0 spirals
for unlike pairs are regularisable in the Coulomb limit, but the
Weber velocity term is unregularised under our two backends — see
CLAUDE.md "Regularization backends"). The net effect is that Zöllner
pushes 2+/2− ICs from marginally-stable (Success) into
numerically-singular (Failure) faster, because the failure mode is
not "too little binding" but "too-close unlike-pair encounter".
The like-pair (+,+) sub-ρ problem is, as the brief anticipated,
**untouched by Zöllner** — the (+,+) internal spiral remains the
terminal failure mechanism in Experiment A.

## 5. Bottom line

Across 20 integrations, no Zöllner setting produced a longer-lived or
more stable orbit than its `a = 0` baseline; in two cases
(F1b breathing square and the high-η rhombus) positive `a` demonstrably
*shortened* the Success window. The single interesting signal is the
monotone reduction in local energy drift for
`nucleus_plus_orbiters_A` with increasing `a`, suggesting that a
symmetric-spectator variant of Experiment A (Agent 8 §Open questions)
might benefit from a small `a ≈ 0.1` nudge — but we cannot test that
without redoing Agent 8's IC sweep, which is outside this agent's
budget.

**Conclusion.** For the 4-body 2+/2− Weber Hamiltonian as currently
implemented, Zöllner electrogravitational mismatch `a > 0` does **not**
stabilise periodic orbits, does **not** enlarge bound basins, and does
**not** regularise close encounters. Its effect is consistent with
"more attraction between unlike pairs" and nothing more — a useful
null result that confines the Agent-13 question.

## Files

- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/13_zollner_comparison/zollner_comparison.jl`
- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/13_zollner_comparison/zollner_comparison.csv`
