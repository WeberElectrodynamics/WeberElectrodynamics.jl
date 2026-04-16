# FW2024 Result-to-Codebase Mapping Table

Each row maps a specific result from Frauenfelder-Weber 2024 to its
manifestation in the WeberElectrodynamics codebase.

Status key:
- **Confirmed**: codebase reproduces or relies on the result
- **Extended**: codebase builds significantly beyond it
- **Implemented**: result is directly coded into the package
- **Open**: result identifies a gap not yet addressed
- **N/A**: result has no direct codebase counterpart

---

## Theorems and Propositions

| FW2024 Result | Statement (abbreviated) | Codebase File(s) | Status | Notes |
|---|---|---|---|---|
| **Theorem A** | Radial Schrodinger eq is limit circle at both r=0 and r=r_c | `CriticalRadiusAndLikeChargeAttraction.md` Section 9 | **Confirmed** | Codebase documents the classical-quantum parallel table. No quantum solver implemented. |
| **Theorem 2.1** | Complete classification of 2-body like-charge trajectories by (h, ell) | `CriticalRadiusAndLikeChargeAttraction.md` Sections 5-8; `TetheringImpossibility.md` eq. (FW) | **Confirmed + Extended** | Energy equation (2.3) is the workhorse for tethering impossibility proof. Classification table reproduced in CriticalRadius doc Section 8. |
| **Theorem 2.1, Case 1** | h <= 0: interior oscillation, no exterior solutions | `CriticalRadiusAndLikeChargeAttraction.md` Section 5 | **Confirmed** | Weber's "molecular state" -- documented as "Two States of Aggregation". |
| **Theorem 2.1, Case 2b** | 0 < h < h_c: interior oscillation in (0, r_+), exterior escape | `CriticalRadiusAndLikeChargeAttraction.md` Section 8 | **Confirmed** | The "impenetrability of the barrier" at r_c. |
| **Theorem 2.1, Case 2c** | h = h_c: unique energy allowing crossing of r_c | `CriticalRadiusAndLikeChargeAttraction.md` Section 8 | **Confirmed** | Documented as "the sole exception to Weber's dictum." |
| **Proposition 4.2** | ell = 0: limit circle at both endpoints | `AngularMomentumRegularization.md` Section "Connection to Quantum Theory" | **Confirmed** | Table mapping classical regularity to quantum natural BCs. |
| **Proposition 4.4** | ell != 0: limit circle at both endpoints | `AngularMomentumRegularization.md` Section "Connection to Quantum Theory" | **Confirmed** | Table mapping classical wild oscillation to quantum wild oscillation. |
| **Lemma 3.1** | Laplace-Beltrami operator formula on Weber plane | (none) | **N/A** | No quantum operator implemented in codebase. |

---

## Key Equations

| FW2024 Equation | Formula | Codebase Location | Status | Notes |
|---|---|---|---|---|
| **(1.1)** Weber Lagrangian | `L = (1/2)(g_rr v_r^2 + r^2 v_phi^2) - 1/r` | `CriticalRadiusAndLikeChargeAttraction.md` Section 7 | **Confirmed** | Metric form documented with identical g_rr = (r-rho)/r. |
| **(2.2)** Weber Hamiltonian | `H = (1/2)(r p_r^2/(r-r_c) + ell^2/r^2) + 1/r` | `src/weber_system.jl` line 125 | **Implemented** (velocity form) | Codebase uses `U = k/r (1 - rdot^2/(2c^2))` -- equivalent but not the metric form. |
| **(2.3)** Energy equation | `rdot^2 = (ell^2 + 2r - 2hr^2)/(r(r_c - r))` | `TetheringImpossibility.md` eq. (FW); `AngularMomentumRegularization.md` Section "The Problem" | **Confirmed** | Used as the primary analytical tool in both investigation documents. |
| **(2.5)** Turning points | `r_pm = (1 +/- sqrt(1 + 2h ell^2))/(2h)` | (not explicitly) | **Confirmed** (implicitly) | Turning point structure used in classification but not computed explicitly in code. |
| **(2.6)** Critical energy | `h_c = ell^2/(2 r_c^2) + 1/r_c` | `CriticalRadiusAndLikeChargeAttraction.md` Section 8 | **Confirmed** | Documented as `h_c = V_eff(rho)`. |
| **(3.7)** Weber metric/cometric | `g_ij = diag((r-r_c)/r, r^2)`, `g^ij = diag(r/(r-r_c), 1/r^2)` | `CriticalRadiusAndLikeChargeAttraction.md` Section 7 | **Confirmed** | Same diagonal metric, same notation. |
| **(3.10)** Weber-Schrodinger eq | `-(1/2) Delta_g psi + (1/r) psi = E psi` | `CriticalRadiusAndLikeChargeAttraction.md` Section 9 | **Confirmed** (described) | Equation stated; no solver implemented. |
| **(4.14)** Sturm-Liouville form (ell=0) | `(p Rdot)' + q R = w E R` | (none) | **N/A** | No ODE solver for S-L problem in codebase. |
| **(4.29)** Bessel parameters | `lambda^2 = (2k/c)^2`, `nu^2 = (1-32r_c)/4` | (none) | **N/A** | Bessel analysis is purely analytical, not implemented. |

---

## Qualitative Results and Observations

| FW2024 Statement | Location in FW2024 | Codebase File(s) | Status | Notes |
|---|---|---|---|---|
| Critical radius = metric signature flip | Section 1, eq. (1.1) | `CriticalRadiusAndLikeChargeAttraction.md` Sections 6-7 | **Confirmed** | mu_eff = mu(1 - rho/r) is exactly g_rr. |
| No periodic orbits inside Weber nucleus | Introduction, p.5 | `AngularMomentumRegularization.md`; Agent 5 numerical search | **Confirmed** | All ell != 0 trajectories terminate; ell = 0 oscillations are not C^1. |
| ell = 0 collision at speed sqrt(2) c | Section 2, p.7 | `CriticalRadiusAndLikeChargeAttraction.md` Section 5; `AngularMomentumRegularization.md` Experiment 0 | **Confirmed** | Numerical verification: final |rdot| = 1.414 ~ sqrt(2). |
| ell = 0 collision is C^0-continuable | Section 2, p.7 | `solve.jl` lines 1088-1141; `CollisionBounceRegularization.md` | **Implemented** | Collision bounce reflects q_rel -> -q_rel at bounce radius. |
| ell != 0 collision at infinite speed | Section 2, eq. (2.4) | `AngularMomentumRegularization.md` Experiment 0 | **Confirmed** | Numerical fit: rdot ~ r^{-0.48} approaching -0.5 theory. |
| ell != 0 collision has infinite winding | Section 2 (implicit via Fig.1-5) | `AngularMomentumRegularization.md` Section "The Topological Obstruction" | **Extended** | FW note it; codebase proves it rigorously (winding number integral). |
| Geometric regularization at r=0 is open | Section 2, p.7 | `AngularMomentumRegularization.md` (7 approaches); `CollisionBounceRegularization.md` | **Extended** | ell=0: Sundman/LC work. ell!=0: all 7 approaches fail. Partially resolves FW's open question. |
| Semiclassical interpretation is open | Introduction, p.5 | (none) | **Open** | Gutzwiller trace formula blocked by no periodic orbits. |
| Legendre transform degenerates at r_c | Implicit in eq. (2.2) | `src/weber_system.jl` | **Avoided** | Codebase uses velocity-dependent potential, bypassing metric formulation. |
| Classical-quantum parallel (Fig. 1) | Introduction, Fig. 1 | `AngularMomentumRegularization.md` Section "Connection to Quantum Theory"; `CriticalRadiusAndLikeChargeAttraction.md` Section 9 | **Confirmed** | Reproduced as a table in both documents. |

---

## Results NOT in FW2024 but developed in codebase

| Codebase Result | File(s) | Relationship to FW2024 |
|---|---|---|
| Non-regularizability theorem (topological: infinite winding number invariant) | `AngularMomentumRegularization.md` | Formalizes what FW state informally. Uses FW's Theorem 2.1 as foundation. |
| Tethering impossibility (external charges cannot stabilize ell!=0 sub-critical pairs) | `TetheringImpossibility.md` | 3 theorems extending FW's 2-body result to the perturbed case. Uses FW energy eq. (2.3). |
| Collision bounce implementation | `solve.jl`, `types.jl` | Practical answer to FW's open question on geometric regularization at r=0. |
| n-body Weber Hamiltonian (symbolic) | `src/weber_system.jl` | FW treats 2-body only; codebase generalizes to arbitrary n. |
| Zollner electrogravitational extension | `src/types.jl` (ZollnerOptions), `src/weber_system.jl` (kappa params) | Not mentioned in FW2024. Adds per-pair coupling constants. |
| 4-body 2+/2- periodic orbit search | `research/FourBodyTwoPlusTwoMinus/` | Far beyond FW2024 scope. Only 1 orbit found (breathing square, super-critical, unstable). |
| Floer homology framing for periodic orbit existence | `10_floer_symplectic/NOTES.md` | Uses FW's critical-radius obstruction as the key barrier in Floer theory. |
| Morse/Conley analysis of effective potential | `12_homology_morse/NOTES.md` | Alternating square is index-2 saddle. FW does not discuss effective potential topology. |
| McGehee blow-up analysis | `AngularMomentumRegularization.md` Experiment 4 | FW does not use blow-up coordinates. Codebase shows collision manifold is invariant only for ell=0. |

---

## Summary Statistics

- **FW2024 formal results**: 1 main theorem (A), 1 classification theorem (2.1), 2 propositions (4.2, 4.4), 1 lemma (3.1), 3 remarks (4.1, 4.3, 4.5)
- **Confirmed in codebase**: 6/6 theorems and propositions
- **Extended by codebase**: 3 results (non-regularizability, tethering, regularization)
- **Implemented in codebase**: 2 results (collision bounce, Hamiltonian)
- **Open/unaddressed**: 2 (semiclassical interpretation, explicit quantum spectrum)
- **Contradicted**: 0
