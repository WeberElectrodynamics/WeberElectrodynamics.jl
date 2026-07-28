### Breaking changes

- **Corrected the canonical Weber Hamiltonian.** The package previously treated
  `p_i = m_i v_i` as the canonical momentum of the Weber Lagrangian. It is not:
  the velocity-dependent pair term contributes a radial correction, so
  `p_i = ∂L/∂v_i = m_i v_i − Σ_{j≠i} (q_i q_j/c²) ṙ_ij (r_i − r_j)/r_ij²`.
  The integrated Hamiltonian was therefore not the Legendre transform of the
  stated Lagrangian, and both `1/c²` signs in the canonical momentum rate were
  reversed. **All Weber trajectories at finite `c` with nonzero pair radial
  velocity change.** Results agree with the previous version in the Coulomb
  limit `c → ∞` and at instants where every pair radial velocity vanishes.
  Upgrade: re-run any simulation whose conclusions depend on finite-`c`
  dynamics.
- `weber_term` has been removed. The default Weber system is now **analytic**:
  recovering physical velocities from canonical momenta requires a coupled
  linear solve over particle pairs at every evaluation, which has no practical
  closed symbolic form for general `n`. `HamiltonianSystem(n, dims)` therefore
  leaves `hamiltonian_symbolic`, `dq_dt_symbolic`, and `dp_dt_symbolic` as
  `nothing`; use `has_symbolic_hamiltonian(sys)` before reading them. Compiled
  function signatures are unchanged, so solving, statistics, plotting, and
  animation are unaffected. The generic symbolic constructor
  `HamiltonianSystem(H, q_vars, p_vars; …)` is unchanged and remains the path
  for custom Hamiltonians.
- `NamedTerm` pair-decomposition closures now take `(q, p, params)` and return
  per-pair **vectors** in `pair_indices` order, instead of `(i, j, q, p, params)`
  returning scalars. Taking the whole state at once matters for
  velocity-dependent Hamiltonians, where physical velocities come from one
  coupled solve over all pairs. `NamedTerm` also gained an optional
  `kinetic_energy(q, p, params)` closure returning physical kinetic energy.
- `EnergyData.kinetic_energy` is now the **physical** kinetic energy
  `Σ ½ mᵢ|vᵢ|²`, not `Σ |pᵢ|²/(2mᵢ)`. `PairEnergyData.radial_velocity` and the
  `PairForceData` velocity/acceleration decompositions are likewise physical.
- `compute_total_kinetic_energy(v, masses, dims)` now takes physical velocities
  instead of momenta; `compute_pair_weber_components(q, v, i, j, charges, c, dims)`
  now takes physical velocities and no longer takes `masses`.
- `two_body_initial_conditions` requires a `c` keyword when `radial_velocity`
  is nonzero, because the conjugate radial momentum is
  `p_r = (μ − q₁q₂/(r c²)) ṙ`, not `μ ṙ`. It throws if the pair sits exactly at
  Weber's critical radius. Zero-radial constructions are unchanged and still do
  not need `c`.
- JLD2 solution archives use format v3. Format-v2 archives were integrated with
  the pre-correction Hamiltonian and are rejected with an explanatory error
  rather than silently reinterpreted as corrected Weber dynamics. Re-run the
  simulation to regenerate them.

### Added

- `physical_velocities(prob, q, p)` and
  `physical_velocities(q, p, params; n_particles, dims)` — the single supported
  way to recover physical velocities from a canonical state. Every Weber path
  in the package now routes through the same internal solve, so no consumer can
  reintroduce `p/m` as a velocity.
- `has_symbolic_hamiltonian(sys)` to distinguish analytic from symbolic systems.
- `WeberCriticalRadiusError`, thrown when the canonical mass matrix is singular
  — for a like-charge pair, exactly at Weber's critical radius
  `ρ = q₁q₂/(μc²)`. Below `ρ` the effective radial inertia is negative but
  finite and integration proceeds normally; only the crossing is singular.
- New documentation pages: **The Weber Hamiltonian** (canonical formulation,
  velocity recovery, critical radius, and a correct/wrong quick-reference table)
  and **Custom Hamiltonians** (the analytic and symbolic construction paths,
  `NamedTerm` hooks, and how to supply your own analytic system).
- `test/test_canonical_weber.jl` — independent physics checks over 1D/2D/3D and
  n = 1…4: canonical momentum against `∂L/∂v`, the exact Legendre transform,
  both canonical equations against finite differences of `H`, recovery of the
  mechanical Weber force, Coulomb and zero-radial limits, translation and
  rotation invariance, and critical-radius behaviour. Every expectation is
  derived from the Lagrangian rather than from the routine under test.

### Changed

- The companion paper is bumped to v1.3 with the corrected canonical momentum,
  Hamiltonian, and both canonical equations, plus a corrected complexity
  analysis: pair geometry is still `O(n²)`, but obtaining physical velocities
  requires an `O(P³)` solve over the `P = n(n−1)/2` pairs.
- `papers/Computational-Weber-Electrodynamics/verify_formulas.py` rewritten
  against the corrected formulation — 40 checks, all passing, including an
  `n = 3` simultaneous velocity solve verified by high-precision finite
  differences. SymPy is now provisioned explicitly via `pyproject.toml` and
  `uv.lock` in the paper directory.
- Theory notes corrected: `WeberElectrodynamics.md` (canonical momenta, the
  exact canonical Hamiltonian, the general velocity solve, and both `ṗ` signs),
  `InitialConditions.md` (canonical vs kinetic momentum, and the corrected
  mid-oscillation radial momentum), and `NonZeroRadialVelocityBoundICs.md`
  (exact canonical Hamiltonian in place of the `O(c⁻⁴)` approximation).
  `theory/HamiltonianCorrectness.md`, which recorded the finding, is removed —
  its standing content now lives in the documentation.
- All seven regression fixtures regenerated from the corrected system, and the
  three example notebooks and three tracked figures re-executed.
- Under the corrected Hamiltonian the `:lifted_pair` backend no longer beats
  plain Cartesian on state error for mild encounters, because the Levi-Civita
  and KS charts regularize the Coulomb singularity but not the Weber
  modification to the pair's effective radial inertia. It does still improve
  energy conservation, and all backends retain second-order convergence. The
  regularization accuracy test now measures against an independent fine-step
  Cartesian reference rather than a lifted-pair one.
