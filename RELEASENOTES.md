### Fixed
- Corrected `WeberProblem` docstring: `regularization_enabled` default was documented as `true` but has been `false` since v0.2.0.

### Changed
- Added inline code comments clarifying two known design choices in the regularization integrator: the frozen monitor (`r_eff`) in the Levi-Civita lifted-pair substep, and the one-pass (diagnostic-only) KS constraint projection in the 3D adaptive-Cartesian backend.
- Added note to `_compute_zollner_kappas` that neutral particles (`q = 0`) are treated as unlike any charged particle due to `sign(0.0)` semantics.
- Added a note to `docs/src/regularization.md` documenting that the 3D KS constraint is enforced by a one-pass projection (not iterative), and that `max_constraint_violation` tracks the residual.
