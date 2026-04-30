### Added

- **3D KS Regularization.** The `:lifted_pair` backend now supports Kustaanheimo-Stiefel (KS) regularization in 3D, providing analytical resolution of close encounters for binary pairs in three dimensions.
- **1D Square-root Regularization.** Added a 1D lifted-pair backend for head-on collisions, completing the dimension-specific regularization suite (1D, 2D, 3D).
- **Parametric Type Refactor.** `HamiltonianProblem`, `HamiltonianSolution`, and `HamiltonianIntegrator` are now parametric on their system and problem types, enabling better compiler specialisation and significantly reducing dynamic dispatch in the integrator inner loops.
- **Sparse Saving and Streaming.** Added `save_stride` (or `save_every`) to `solve`/`init` to control the frequency of state snapshots, and a `stream_sink` callback to process states in real-time without storing the full trajectory.
- **Conservation Summary.** New `conservation_summary(sol)` utility providing a compact report of energy errors, momentum drift, and regularization diagnostics for quick verification.
- **Solution Archiving.** Integrated `save_solution` and `load_solution` helpers (requiring JLD2.jl) for easy persistence of simulation results including all metadata.

### Fixed

- **Regularization substepping logic.** Corrected incorrect time-offset accumulation in the `RegularizedIntegrator` substepping loop.
- **LC map singularity at origin.** The Levi-Civita projection `_lc_project!` now preserves momentum when `r ≈ 0`, preventing numerical NaN during exact head-on encounters.
- **Zöllner κ logic for neutral particles.** Charge pairs where at least one particle is neutral (`q = 0`) now correctly receive `κ = 1.0` instead of `1 + a`.
- **Animation dashboard stability.** Fixed several issues in `animate_weber` circular buffer handling and phase-space selection.
- **Kappa accessor robustness.** `kappa(prob, i, j)` now automatically handles index ordering (supports `j < i`).

### Changed

- **Regularization diagnostics validation.** Regression fixtures now include 16 fields of regularization diagnostics (substep counts, backend fallbacks, constraint violations), ensuring numerical equivalence of the regularization state.
- **Trail length slider.** The `animate_weber` trail-length slider now includes a log-linear range.
- **Release script safety.** `release.sh` now verifies the active `gh` account and uses a robust Python-based version bump.
