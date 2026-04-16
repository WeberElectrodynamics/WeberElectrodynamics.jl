# Technique Transfer Table: Frauenfelder's Works to Weber Hamiltonian

Agent 06 deliverable. See NOTES.md for full discussion of each entry.

## Legend

- **Applies?**: Yes / Partial / No / Conditional
- **Difficulty**: Easy / Moderate / Hard / Very Hard
- **Payoff**: Low / Moderate / High / Very High
- **Obstruction codes**:
  - PDM = Position-Dependent Metric (Weber's G(q) kinetic energy)
  - MD  = Metric Degeneracy at r = rho
  - NC  = Noncompactness of energy surfaces
  - NI  = Non-Integrability (n >= 3 body)
  - LOR = Lorentzian signature (r < rho)

---

## Table

| # | Paper | Year | Technique | Applies to Weber? | Obstruction | Difficulty | Payoff |
|---|-------|------|-----------|-------------------|-------------|------------|--------|
| 1 | Frauenfelder-J.Weber, "Fine structure of Weber's hydrogen atom" | 2019 | Bohr-Sommerfeld/EBK quantization of Weber H | **Yes** (unlike-charge) | Unlike-charge only; requires known orbits | Easy | High |
| 2 | Frauenfelder-J.Weber, "Weber nucleus" (FW2024) | 2024 | Lorentzian metric interpretation; Sturm-Liouville | **Yes** (negative result) | Proves no periodic orbits for r < rho, ell != 0 | N/A | Very High (constraint) |
| 3 | Frauenfelder, "Regularization of the Kepler problem" | 2012 | Moser regularization (stereographic -> geodesic on S^n) | **Partial** (far field only) | PDM near rho; LOR inside rho | Moderate | High |
| 4 | Frauenfelder-van Koert-Zhao, "Convex embedding for rotating Kepler" | 2016 | Ligon-Schaaf + Levi-Civita convex embedding | **Conditional** (large c) | PDM; MD at rho prevents global convexity | Hard | Very High |
| 5 | Albers-Fish-Frauenfelder-van Koert, "CZ indices of rotating Kepler" | 2013 | Explicit CZ index computation | **Yes** (for continued orbits) | Need orbits first | Moderate | High |
| 6 | Albers-Frauenfelder-van Koert-Paternain, "Contact geometry of R3BP" | 2012 | Contact type proof for energy surfaces | **Conditional** | PDM breaks standard Liouville form | Very Hard | Very High |
| 7 | Albers-Fish-Frauenfelder-Hofer-van Koert, "Global surfaces of section" | 2011 | Disk-like global surface of section for Reeb flow | **Conditional** (needs contact type) | Depends on #6 | Very Hard | Very High |
| 8 | Frauenfelder-van Koert, "R3BP and Holomorphic Curves" (book) | 2018 | Full holomorphic curve framework for celestial mechanics | **Partial** (framework) | PDM throughout | Hard | Very High |
| 9 | Frauenfelder-Kang, "Real holomorphic curves and invariant GSS" | 2016 | Symmetric periodic Reeb orbits via real J-curves | **Conditional** (needs dyn. convexity) | PDM; need convexity | Hard | High |
| 10 | Cieliebak-Frauenfelder-Paternain, "Mane critical values" | 2010 | RFH on twisted cotangent bundles; critical energy analysis | **Partial** (via magnetic analogy) | Weber is NOT a twisted cotangent bundle directly | Very Hard | Very High |
| 11 | Frauenfelder-Merry-Paternain, "Floer hom. for magnetic fields" | 2012 | Floer homology for magnetic cotangent bundles | **Conditional** (if magnetic reduction works) | Need to recast Weber as magnetic flow | Very Hard | High |
| 12 | Frauenfelder-Merry-Paternain, "Non-resonant magnetic fields on flat tori" | 2015 | Novikov Floer homology; multiplicity of periodic orbits | **Conditional** | Specific to flat tori; need magnetic reduction | Very Hard | Moderate |
| 13 | Frauenfelder-Schlenk, "Hamiltonian dynamics on convex symplectic mfds" | 2007 | Floer=Morse isomorphism; periodic points; Weinstein conj. | **Conditional** (if r > rho region is convex symplectic) | Need to establish convexity of Weber phase space | Hard | High |
| 14 | Cieliebak-Frauenfelder-Volkov, "Frozen planet orbits in helium" | 2022 | Nonlocal LC regularization; multi-time-scale; variational existence | **Partial** (structural analog) | Weber velocity-dependence not regularized by LC | Hard | Very High |
| 15 | Cieliebak-Frauenfelder-Volkov, "Nondegeneracy of frozen planets" | 2023 | Nondegeneracy proof for variational critical points | **Conditional** (needs #14 first) | Depends on #14 | Hard | High |
| 16 | Frauenfelder, "Compactness for frozen planets" | 2023 | Compactness of moduli space | **Conditional** (needs #14) | NC; depends on #14 | Hard | Moderate |
| 17 | Albers-Frauenfelder-Schlenk, "Hamiltonian delay equations" | 2020 | Graph trick + Lagrangian Floer hom.; Arnold conj. with delay | **No** (indirect inspiration) | Weber is not a delay equation | Very Hard | Low |
| 18 | Frauenfelder, "Helium and Hamiltonian delay equations" | 2021 | Multi-time-scale -> delay formulation of helium | **Partial** (multi-time idea) | Conceptual, not direct | Very Hard | Moderate |
| 19 | Frauenfelder-Moreno, "Doubly symmetric periodic orbits" | 2023 | No negative hyperbolicity; CZ parity determines stability | **Yes** | None (applies to any 4D Hamiltonian) | Easy | Moderate |
| 20 | Frauenfelder-Koh-Moreno, "Symplectic methods in numerical orbit search" | 2023 | Cell mapping; CZ index tracking; Krein-Moser signs | **Yes** | None (general numerical method) | Easy-Moderate | Very High |
| 21 | Moreno-Aydin-Frauenfelder-van Koert-Koh, "Bifurcation graphs for CR3BP" | 2024 | CZ index algorithm; bifurcation graph construction | **Yes** | None (general numerical method) | Moderate | Very High |
| 22 | Cieliebak-Frauenfelder-Oancea, "RFH and symplectic homology" | 2010 | Long exact sequence: SH -> SH -> RFH -> SH | **Conditional** | NC; PDM | Very Hard | High |
| 23 | Frauenfelder, "Stark problem as concave toric domain" | 2023 | Toric domain structure after LC regularization | **Conditional** (if Weber is integrable in some limit) | NI for n >= 3; PDM | Hard | Moderate |
| 24 | Cieliebak-Frauenfelder, "EBK quantization" | 2025 | Variational EBK from action spectrum on toric domains | **Conditional** (needs toric structure) | NI; PDM | Hard | Moderate |
| 25 | Cieliebak-Frauenfelder-van Koert, "Finsler geometry of rotating Kepler" | 2014 | Finsler metric interpretation of rotating Kepler | **Partial** | Weber's metric is more general than Finsler | Hard | Moderate |
| 26 | Frauenfelder-Zhao, "Periodic collisional orbit or infinitely many" | 2019 | Floer homology forced existence: collisional or infinitely many | **Conditional** | Need regularization framework | Hard | High |
| 27 | Cieliebak-Frauenfelder-Miranda-Wisniewska, "Two-boost problem" | 2025 | Lagrangian RFH for noncompact energy surfaces | **Partial** (noncompactness technique) | PDM; MD | Very Hard | High |
| 28 | Kang, "Generalized RFH and coisotropic intersections" | 2013 | RFH for coisotropic submanifolds | **Conditional** | Need Weber-adapted setup | Hard | Moderate |

---

## Summary Statistics

- **Directly applicable (Yes)**: Papers 1, 2, 5, 19, 20, 21 (6 papers)
- **Conditionally applicable**: Papers 4, 6, 7, 9, 11, 12, 13, 15, 16, 22, 23, 24, 26, 28 (14 papers)
- **Partially applicable**: Papers 3, 8, 10, 14, 18, 25, 27 (7 papers)
- **Not directly applicable**: Paper 17 (1 paper)

## Top 5 Most Actionable for Weber (ranked by payoff/difficulty ratio)

1. **Paper 20** (Frauenfelder-Koh-Moreno 2023): CZ-index numerical methods. Implement NOW.
2. **Paper 19** (Frauenfelder-Moreno 2023): Doubly symmetric orbit constraints. Apply to any found orbits.
3. **Paper 21** (Moreno et al. 2024): Bifurcation graph algorithm. Extend to Weber families.
4. **Paper 1** (FW2019): Bohr-Sommerfeld for Weber. Already done for unlike-charge; adapt to like-charge perturbatively.
5. **Paper 5** (CZ indices 2013): Track CZ indices of continued Coulomb orbits into Weber regime.

## Top 5 Highest Theoretical Payoff (regardless of difficulty)

1. **Paper 6** (Contact type for energy surfaces): Would unlock ALL holomorphic curve tools.
2. **Paper 10** (Mane critical values / magnetic analogy): Would connect Weber to rich magnetic flow theory.
3. **Paper 14** (Frozen planet variational method): Multi-body existence proofs via nonlocal regularization.
4. **Paper 4** (Convex embedding persistence): Dynamical convexity -> guaranteed periodic orbits.
5. **Paper 7** (Global surface of section): Reduces orbit existence to Brouwer fixed point theorem.
