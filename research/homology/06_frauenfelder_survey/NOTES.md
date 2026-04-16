# Frauenfelder Survey: Techniques Transferable to the Weber Hamiltonian

Agent 06 deliverable. Cross-refs: Agent 01 (FW2024 exegesis), Agent 08
(regularization/Floer), Agent 05 (two-body census).

## 0. Context: The Weber Problem

The Weber two-body Hamiltonian in polar relative coordinates:

```
H = (r / (2mu(r - rho))) p_r^2  +  ell^2 / (2mu r^2)  +  e^2/r
```

where rho = e^2/(mu c^2) is the critical radius. Key features:
- Kinetic metric G(q) is position-dependent: Riemannian for r > rho, Lorentzian for r < rho
- No periodic orbits inside the critical radius (spiraling to collision; FW2024)
- Only 1 known (unstable) periodic orbit in the 4-body system
- System reduces to Coulomb as c -> infinity

---

## 1. DIRECTLY WEBER-RELEVANT: Frauenfelder-Weber Papers

### 1a. "The fine structure of Weber's hydrogen atom" (FW2019)
- **Ref**: Frauenfelder & J. Weber, Z. Angew. Math. Phys. 70, 105 (2019). arXiv:1902.09612
- **Technique**: Bohr-Sommerfeld (EBK) quantization of Weber's Hamiltonian
- **Key result**: Energy levels to second order in fine structure constant match Wesley's
  Schrodinger calculation. Identifies quantized tori in phase space.
- **Relevance**: **Direct**. This paper works with exactly the Weber Hamiltonian we study.
  The quantized tori are the periodic/quasi-periodic orbits we seek. The perturbative
  approach (expanding in 1/c^2 around Coulomb) is a proven strategy for finding
  near-Coulombic bound orbits.
- **Transferable technique**: Action-angle variable analysis; perturbation from Coulomb limit.
  The paper confirms that for the hydrogen-like (unlike-charge, attractive) case, bound
  orbits exist and can be characterized semiclassically.
- **Limitation for our problem**: This treats UNLIKE charges (attractive Coulomb). Our
  hardest case is LIKE charges where the Coulomb repulsion must be overcome by Weber
  velocity-dependent attraction. The Bohr-Sommerfeld approach requires knowing the orbits
  first.

### 1b. "A mathematical description of the Weber nucleus" (FW2024)
- **Ref**: Frauenfelder & J. Weber, Anal. Math. Phys. 14, 31 (2024)
- **Technique**: Lorentzian metric interpretation of Weber's Lagrangian; Weyl singular
  Sturm-Liouville theory for quantum mechanics
- **Key result**: Inside the critical radius (r < rho), the kinetic energy has Lorentzian
  signature. There are NO periodic orbits for ell != 0; trajectories spiral into
  r = 0 (collision). The angular velocity diverges.
- **Relevance**: **Critical negative result**. This is the FW2024 paper that Agent 01
  analyzed. Establishes that sub-critical like-charge dynamics is fundamentally
  non-periodic. Any bound orbits must live in the super-critical region r > rho.
- **Transferable insight**: The Lorentzian-to-Riemannian transition at r = rho is a
  metric degeneracy that no standard regularization can remove. Agent 08's McGehee
  blow-up confirms the topological obstruction for ell != 0. Bound orbits, if they
  exist, must either (a) stay entirely in r > rho, or (b) have ell = 0 and undergo
  head-on collision bounces.

---

## 2. MOSER REGULARIZATION AND THE KEPLER PROBLEM

### 2a. "On the regularization of the Kepler problem" (FW, J. Sympl. Geom. 2012)
- **Ref**: Frauenfelder (solo, but building on Moser 1970). arXiv:1007.3695
- **Technique**: Moser regularization = stereographic projection mapping Kepler flow on
  fixed negative energy to geodesic flow on S^n. The Ligon-Schaaf map is shown to be
  an adaptation of Moser's construction, explaining the hidden SO(4) symmetry.
- **Applicability to Weber**: **Partially, with major caveats**. Moser regularization
  requires the Hamiltonian to be H = |p|^2/(2m) + V(q) with STANDARD kinetic energy.
  Weber's position-dependent kinetic metric G(q) breaks this. Specifically:
  - For r >> rho (far from critical radius): G(q) ~ I + O(rho/r), nearly flat.
    Moser regularization applies perturbatively. One could regularize the Coulomb
    part and treat the Weber correction as a perturbation of the geodesic flow on S^3.
  - Near r = rho: G(q) degenerates. Moser regularization fails completely.
  - Inside r < rho: Metric is Lorentzian. Moser's construction is fundamentally
    Riemannian and does not extend.
- **Difficulty**: Moderate (for perturbative far-field); impossible (near critical radius)
- **Payoff**: Would give rigorous existence of near-circular bound orbits in the
  large-r regime via continuation from Coulomb/Kepler.

### 2b. "A convex embedding for the rotating Kepler problem" (Frauenfelder-van Koert-Zhao, 2016)
- **Ref**: arXiv:1605.06981
- **Technique**: Combines Ligon-Schaaf and Levi-Civita regularizations to produce a
  convex symplectic embedding of rotating Kepler energy surfaces into R^4. Proves
  dynamical convexity.
- **Applicability to Weber**: **Moderate**. In the c -> infinity limit, Weber reduces
  to Coulomb/Kepler. For large c, the Weber energy surfaces are perturbations of
  Kepler energy surfaces. If one could show the convex embedding persists under the
  Weber perturbation, dynamical convexity would guarantee periodic orbits via
  holomorphic curve theory. The obstruction is that for finite c, the metric
  degeneracy at r = rho prevents global convexity.
- **Difficulty**: Hard
- **Key question**: For what range of c (or energy) does the convex embedding survive?

### 2c. "The Conley-Zehnder indices of the rotating Kepler problem" (Albers-Fish-Frauenfelder-van Koert, 2013)
- **Ref**: Math. Proc. Cambridge Phil. Soc. 154 (2013)
- **Technique**: Explicit CZ index computation for Kepler periodic orbits in rotating frame
- **Applicability**: If we find Weber periodic orbits by continuation from Coulomb, their
  CZ indices can be tracked. CZ index parity determines stability type (per
  Frauenfelder-Moreno 2023). This is directly useful for Agent 04's work.
- **Difficulty**: Moderate (once orbits are found)

---

## 3. RESTRICTED THREE-BODY PROBLEM AND CONTACT TOPOLOGY

### 3a. "Contact geometry of the restricted 3-body problem" (Albers-Frauenfelder-van Koert-Paternain, CPAM 2012)
- **Ref**: Comm. Pure Appl. Math. 65, 229-263 (2012). arXiv:1010.2140
- **Technique**: Proves the energy hypersurface of the planar CR3BP is of (restricted)
  contact type below the first critical energy value. This means the Reeb flow exists
  and standard contact-topological tools apply.
- **Applicability to Weber 2+1 problem**: **High potential, major obstruction**.
  The Weber "restricted" problem (test particle + 2 fixed/heavy charges) has energy
  surfaces that may be of contact type in the super-critical region. However:
  - The position-dependent kinetic metric means the standard Liouville form
    lambda = p dq does NOT give a contact form on energy surfaces.
  - One needs a MODIFIED Liouville form adapted to the metric G(q).
  - For the Jacobi metric formulation: H = ½ p^T G(q)^{-1} p + V(q) at energy h gives
    the Jacobi metric g_J = (h - V)G on configuration space. Contact type of the
    energy surface is related to convexity of g_J.
- **Difficulty**: Very hard
- **Payoff**: Would unlock the entire machinery of holomorphic curves and Floer homology

### 3b. "Global surfaces of section in the planar restricted 3-body problem" (Albers-Fish-Frauenfelder-Hofer-van Koert, ARMA 2011)
- **Ref**: Arch. Ration. Mech. Anal. 199, 411-436 (2011)
- **Technique**: Constructs a disk-like global surface of section for the Reeb flow on
  the regularized energy surface. Reduces 3D flow to 2D area-preserving map.
- **Applicability**: If contact type is established for Weber, a global surface of
  section would reduce the periodic orbit problem to finding fixed points of an
  area-preserving map (guaranteeing existence by Brouwer). This is the "royal road"
  to existence proofs.
- **Difficulty**: Very hard (depends on 3a)

### 3c. Frauenfelder-van Koert book: "The Restricted Three-Body Problem and Holomorphic Curves" (Springer, 2018)
- **Ref**: Pathways in Mathematics, Springer
- **Content**: Comprehensive treatment: Kepler problem, CR3BP, Moser regularization,
  holomorphic curves in symplectizations, numerical orbit search. Chapters 1-5 cover
  the classical mechanics; Chapters 6-10 develop holomorphic curve theory; Chapters
  11-14 apply it to CR3BP.
- **Applicability**: **Essential reference**. The entire framework developed in this book
  is what one would attempt to adapt to the Weber problem. The key question is always:
  does the position-dependent metric obstruct the construction?

### 3d. "Real holomorphic curves and invariant global surfaces of section" (Frauenfelder-Kang, Proc. LMS 2016)
- **Ref**: Proc. London Math. Soc. 112, 477-511 (2016)
- **Technique**: Uses real (anti-holomorphic-involution-invariant) pseudoholomorphic
  curves to find symmetric periodic orbits and invariant global surfaces of section.
  Key result: dynamically convex + symmetric => 2 or infinitely many symmetric
  periodic Reeb orbits.
- **Applicability to Weber**: **Moderate-High**. Weber's equations have time-reversal
  symmetry (t -> -t, p -> -p). This is exactly the kind of anti-symplectic involution
  that Frauenfelder-Kang exploit. If the energy surface is dynamically convex (big if),
  this guarantees multiple symmetric periodic orbits.
- **Difficulty**: Hard (need dynamical convexity first)
- **Payoff**: High. Symmetric orbits are the most natural candidates for bound states.

---

## 4. TWISTED COTANGENT BUNDLES AND MAGNETIC FLOWS

### 4a. "Symplectic topology of Mane's critical values" (Cieliebak-Frauenfelder-Paternain, Geom. & Topol. 2010)
- **Ref**: Geom. Topol. 14, 1765-1870 (2010). arXiv:0903.0700
- **Technique**: Studies energy hypersurfaces of mechanical Hamiltonians on TWISTED
  cotangent bundles: H = ½|p|^2 on (T*M, omega_can + pi*sigma) where sigma is a
  closed 2-form (magnetic field). Uses RFH to study periodic orbits above/below
  the Mane critical value.
- **KEY INSIGHT FOR WEBER**: The Weber Hamiltonian H = ½ p^T G(q)^{-1} p + V(q) is
  NOT on a twisted cotangent bundle (no magnetic term), but it IS a Hamiltonian with
  position-dependent kinetic energy. However, after Jacobi-Maupertuis reduction at
  energy h, the geodesic flow of the Jacobi metric g_J = (h-V)G on configuration
  space can sometimes be rewritten as a magnetic flow on a STANDARD metric with a
  twist term encoding the velocity dependence. This is the magnetic analogy!
- **How it could work**: Write Weber's kinetic metric as G(q) = I + delta G(q) where
  delta G encodes the Weber correction. In the Jacobi picture, this might decompose
  into a standard geodesic flow plus a "magnetic-like" correction.
- **Critical value analysis**: The Mane critical value separates contact-type energy
  surfaces (above) from non-contact-type (below). For Weber, the critical radius rho
  might play an analogous role to Mane's critical value.
- **Difficulty**: Very hard
- **Payoff**: Very high -- would connect Weber to the rich theory of magnetic flows

### 4b. "Floer homology for magnetic fields with at most linear growth on the universal cover" (Frauenfelder-Merry-Paternain, J. Funct. Anal. 2012)
- **Ref**: J. Funct. Anal. 262, 3062-3090 (2012). arXiv:1108.3044
- **Technique**: Extends Floer homology to magnetic cotangent bundles where the magnetic
  form has at most linear growth. Proves existence of periodic orbits on almost every
  energy level.
- **Applicability**: If Weber can be recast as a magnetic flow (see 4a), this gives
  periodic orbit existence on "almost every" energy level.

### 4c. "Floer homology for non-resonant magnetic fields on flat tori" (Frauenfelder-Merry-Paternain, Nonlinearity 2015)
- **Ref**: Nonlinearity 28, 1351-1370 (2015). arXiv:1305.3141
- **Technique**: Defines Novikov Floer homology for non-resonant magnetic fields.
  Proves >= 2N+1 contractible periodic solutions (generically 2^{2N}).
- **Applicability**: Specific to flat tori, but the Novikov Floer approach could extend
  to other configuration spaces if the magnetic analogy works.

### 4d. "Hamiltonian dynamics on convex symplectic manifolds" (Frauenfelder-Schlenk, Israel J. Math. 2007)
- **Ref**: Israel J. Math. 159, 1-56 (2007)
- **Technique**: Establishes isomorphism between Floer and Morse homology on convex
  symplectic manifolds. Applications: infinitely many periodic points, Weinstein
  conjecture, closed trajectories of charge in magnetic field on almost all small
  energy levels.
- **Applicability**: **Directly relevant IF** the Weber phase space (restricted to
  r > rho) can be given the structure of a convex symplectic manifold. The "charge in
  magnetic field" result is structurally analogous to Weber.
- **Difficulty**: Hard
- **Payoff**: High -- gives multiplicity results for periodic orbits

---

## 5. FROZEN PLANETS AND HELIUM: VARIATIONAL METHODS

### 5a. "A variational approach to frozen planet orbits in helium" (Cieliebak-Frauenfelder-Volkov, Ann. IHP 2022)
- **Ref**: Ann. Inst. H. Poincare Anal. Non Lineaire 40, 379-455 (2023). arXiv:2103.15485
- **Technique**: Nonlocal Levi-Civita regularization with DIFFERENT time reparametrizations
  for each electron. Deformation from full helium to mean-field interaction.
  Variational existence proof via Lusternik-Schnirelmann theory.
- **Applicability to Weber**: **Highly instructive**. The helium atom has two electrons
  interacting with a nucleus AND each other -- similar multi-body structure to our
  n-body Weber problem. The key innovation of using DIFFERENT time scales for
  different particles is directly relevant because:
  - In Weber's n-body problem, different pairs have different critical radii
  - Near-collision regularization requires pair-specific time reparametrization
  - The nonlocal functional approach bypasses the need for global regularization
- **Difficulty**: Hard
- **Payoff**: Very high for existence proofs in the 3+ body Weber problem

### 5b. "Nondegeneracy and integral count of frozen planet orbits in helium" (Cieliebak-Frauenfelder-Volkov, Tunisian J. Math 2023)
- **Ref**: Tunisian J. Math. 5, 195-242 (2023)
- **Technique**: Proves nondegeneracy of the critical points found in 5a, enabling an
  integral count of orbits.
- **Applicability**: If the variational approach of 5a is adapted to Weber, this paper
  provides the follow-up: proving the found orbits are nondegenerate.

### 5c. "A compactness theorem for frozen planets" (Frauenfelder, J. Topol. Anal. 2023)
- **Technique**: Compactness result for the moduli space of frozen planet orbits.
  Essential for the variational approach.
- **Applicability**: Technical prerequisite for extending the variational method to Weber.

---

## 6. HAMILTONIAN DELAY EQUATIONS

### 6a. "Hamiltonian delay equations -- examples and a lower bound for the number of periodic solutions" (Albers-Frauenfelder-Schlenk, Adv. Math. 2020)
- **Ref**: Adv. Math. 373, 107319 (2020). arXiv:1802.07453
- **Technique**: Extends Arnold conjecture to delay Hamiltonians using graph trick +
  Lagrangian Floer homology. Lower bound on periodic solutions = sum of Betti numbers.
- **Applicability to Weber**: **Indirect but suggestive**. Weber is NOT a delay equation
  (it uses instantaneous velocities, not retarded ones). However:
  - The original Weber electrodynamics of the 19th century was formulated as an
    action-at-a-distance theory. If one introduces finite propagation speed, one gets
    a retarded (delay) version.
  - The graph trick (embedding delay dynamics into a product symplectic manifold)
    could potentially be used to embed Weber's velocity-dependent dynamics into a
    larger phase space where standard tools apply.
  - Frauenfelder's "Helium and Hamiltonian delay equations" (Israel J. Math. 2021)
    explicitly connects helium (which has Coulombic interactions like Weber) to delay
    equations.
- **Difficulty**: Very hard (conceptual leap needed)
- **Payoff**: Moderate (indirect approach)

### 6b. "Helium and Hamiltonian delay equations" (Frauenfelder, Israel J. Math. 2021)
- **Ref**: Israel J. Math. 246, 323-347 (2021)
- **Technique**: Shows how the helium problem can be formulated as a Hamiltonian delay
  equation by using different time parametrizations for the two electrons.
- **Applicability**: The multi-time-scale approach is the same one used in the frozen
  planet work (5a) and is directly relevant to Weber's n-body problem.

---

## 7. DOUBLY SYMMETRIC ORBITS AND STABILITY

### 7a. "On doubly symmetric periodic orbits" (Frauenfelder-Moreno, Celest. Mech. Dyn. Astron. 2023)
- **Ref**: Celest. Mech. Dyn. Astron. 135, 20 (2023). arXiv:2301.01803
- **Technique**: In 4D Hamiltonian systems, doubly symmetric periodic orbits (invariant
  under two commuting anti-symplectic involutions) cannot be negative hyperbolic.
  Stable iff CZ-index is odd. No period-doubling bifurcation.
- **Applicability to Weber**: **Directly applicable**. The Weber 2-body problem in the
  plane has 4D phase space and two natural symmetries:
  (1) Time-reversal: (q, p, t) -> (q, -p, -t)
  (2) Reflection: (r, phi, p_r, p_phi) -> (r, -phi, -p_r, p_phi) (for zero total
      angular momentum configurations)
  Any periodic orbit symmetric under both involutions is guaranteed to be either
  elliptic or positive hyperbolic (never negative hyperbolic). This constrains the
  bifurcation behavior.
- **Difficulty**: Easy to apply once orbits are found
- **Payoff**: Moderate (constrains stability type, guides numerical search)

### 7b. "Symplectic methods in the numerical search of orbits" (Frauenfelder-Koh-Moreno, SIAM J. Appl. Dyn. Syst. 2023)
- **Ref**: SIAM J. Appl. Dyn. Syst. 22, 1533-1570 (2023). arXiv:2206.00627
- **Technique**: Practical numerical methods using CZ indices to: (1) predict orbit
  existence during bifurcations, (2) connect families of periodic orbits, (3) attach
  Krein-Moser signs to Floquet multipliers.
- **Applicability**: **Directly applicable to Weber numerics**. The cell-mapping method
  and CZ-index tracking can be applied to the Weber problem's Poincare map. This is
  probably the MOST IMMEDIATELY USEFUL paper for our computational work.
- **Difficulty**: Easy-Moderate (numerical implementation)
- **Payoff**: Very high (practical orbit search tool)

### 7c. "Bifurcation graphs for the CR3BP via symplectic methods" (Moreno-Aydin-Frauenfelder-van Koert-Koh, J. Astronaut. Sci. 2024)
- **Ref**: J. Astronaut. Sci. 71, 51 (2024)
- **Technique**: Algorithm for numerical CZ index computation; bifurcation graph analysis
  for families of periodic orbits.
- **Applicability**: Extension of 7b with more detailed numerical tools. The CZ-index
  algorithm is directly usable for Weber orbit families.

---

## 8. TORIC DOMAINS AND EBK QUANTIZATION

### 8a. "The Stark problem as a concave toric domain" (Frauenfelder, Geom. Dedicata 2023)
- **Technique**: After Levi-Civita regularization, the Stark problem's energy surfaces
  become boundaries of concave toric domains. This enables computation of symplectic
  capacities and EBK spectra.
- **Applicability**: If Weber's energy surfaces (after regularization in the super-critical
  region) can be expressed as toric domains, one gets quantitative capacity bounds.
  Obstruction: Weber's velocity-dependence may prevent the integrability needed for
  toric structure.

### 8b. "A variational characterization of EBK quantization" (Cieliebak-Frauenfelder, J. Sympl. Geom. 2025)
- **Technique**: Constructs EBK spectrum from the marked action spectrum. Minimax formula
  for concave toric domains.
- **Applicability**: If toric domain structure is established (big if), this gives a
  semiclassical quantization of Weber orbits -- connecting to FW2019's Bohr-Sommerfeld
  approach but more rigorously.

---

## 9. TWO-FIXED-CENTER PROBLEMS AND THE WEBER ANALOG

### 9a. J+ invariants for Stark-Zeeman systems (Cieliebak-Frauenfelder-Zhao, Ergod. Th. Dyn. Syst. 2022)
- **Technique**: Arnold's J+ knot invariant tracks how families of periodic orbits in
  two-center problems connect and bifurcate through homotopy.
- **Applicability**: The Weber 2+1 problem (test particle + 2 fixed charges) is a
  two-center problem. The J+ invariant formalism can classify orbit families.

### 9b. Euler problem convex embedding (Kim, building on Frauenfelder, 2018)
- **Technique**: Convex embedding for the Euler two-fixed-center problem.
- **Applicability**: The Weber two-center analog (with velocity-dependent potential)
  is a natural extension. If convexity survives the Weber perturbation, Floer theory
  applies.

---

## 10. THE TWO-BOOST PROBLEM AND LAGRANGIAN RFH

### 10a. "The two-boost problem and Lagrangian RFH" (Cieliebak-Frauenfelder-Miranda-Wisniewska, 2024/2025)
- **Ref**: arXiv:2412.08415
- **Technique**: Defines Lagrangian Rabinowitz Floer homology to solve the two-boost
  problem (connecting two phase-space points with bounded energy). Main technical
  challenge: noncompactness of energy hypersurfaces.
- **Applicability**: The Weber problem's noncompact energy surfaces (especially near
  the critical radius) are a major technical obstruction. This paper's methods for
  handling noncompactness are directly relevant.
- **Difficulty**: Very hard (cutting-edge technique)
- **Payoff**: High if adapted -- would solve existence questions for Weber connecting orbits

---

## 11. COLLABORATOR WORKS

### 11a. Jungsoo Kang: Generalized RFH and coisotropic intersections (IMRN 2013)
- Extends RFH to coisotropic submanifolds. The Weber energy surface intersected with
  angular momentum constraints is coisotropic. Could give leafwise intersection results
  (= periodic orbits on fixed angular momentum surfaces).
- **Difficulty**: Hard

### 11b. Jungsoo Kang: Magnetic geodesics on surfaces
- Studies closed magnetic geodesics via systolic inequalities. If Weber reduces to a
  magnetic-geodesic problem (see Section 4), Kang's bounds apply.
- **Difficulty**: Moderate (if the magnetic reduction works)

### 11c. Otto van Koert: Contact geometry of spatial CR3BP
- Extends the planar contact-type results to 3D. Relevant for 3D Weber simulations.
- **Difficulty**: Very hard

### 11d. Felix Schlenk: Symplectic capacities and embedding problems
- Capacities give quantitative bounds on orbit periods and actions. Schlenk's
  work on embedding problems could constrain which Weber energy surfaces admit
  periodic orbits.
- **Difficulty**: Hard

### 11e. Agustin Moreno: CR3BP modern symplectic viewpoint
- Recent survey (EMS Magazine) synthesizing the symplectic approach to CR3BP.
  Excellent road map for what to adapt to Weber.

---

## 12. SYNTHESIS: PRIORITY-ORDERED RESEARCH DIRECTIONS

### Tier 1: Immediately actionable (Easy-Moderate)
1. **CZ-index numerical methods** (7b, 7c): Apply Frauenfelder-Koh-Moreno's cell-mapping
   and CZ-tracking to Weber periodic orbit search. This can be implemented NOW.
2. **Doubly symmetric orbit constraints** (7a): Apply to any orbits found computationally.
3. **Perturbation from Coulomb** (2a, FW2019): Use Moser regularization in the large-r
   limit and continue orbits to finite c.

### Tier 2: Medium-term research (Hard)
4. **Convex embedding persistence** (2b, 3d): Can the rotating Kepler convex embedding
   survive the Weber perturbation for large c? This would give existence results.
5. **Magnetic flow reduction** (4a, 4b): Attempt to recast Weber as a magnetic geodesic
   problem. If successful, unlocks the entire Cieliebak-Frauenfelder-Paternain toolkit.
6. **Frozen planet variational approach** (5a): Adapt nonlocal Levi-Civita regularization
   with multi-time-scale to Weber 3-body problem.

### Tier 3: Long-term foundational (Very Hard)
7. **Contact type for Weber energy surfaces** (3a, 3b): The "holy grail" -- would
   guarantee periodic orbits and global surfaces of section.
8. **Lagrangian RFH for noncompact Weber surfaces** (10a): Handle the metric degeneracy
   at r = rho.
9. **Novikov Floer homology for Weber** (4c): Full Floer-theoretic treatment if the
   magnetic reduction is established.

### Key Obstructions (summarized)
- **Position-dependent kinetic metric**: The single biggest obstacle. Almost all of
  Frauenfelder's machinery assumes H = |p|^2/2 + V(q) or at most a twisted symplectic
  form. Weber's G(q)-dependent kinetic energy is more general.
- **Metric degeneracy at r = rho**: The Riemannian-to-Lorentzian transition is not
  handled by any existing technique in symplectic topology.
- **Noncompactness**: Weber's phase space is noncompact even after regularization.
- **No integrability**: The n >= 3 body Weber problem is not integrable, preventing
  toric domain methods.
