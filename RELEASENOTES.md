### Fixed

- `RegularizationDiagnostics.used_backend` now correctly reports `:disabled` when
  regularization is enabled but no encounter is detected during the run (previously it
  incorrectly reported the effective backend even when no regularized step was ever taken).
- Removed dead 2D chain-mode LC lifts in `_step_regularized_chain!`: per-substep
  Levi-Civita lifts for each chain edge in 2D were overwriting the same buffers on every
  loop iteration and never being read, wasting two compiled-RHS calls per substep per edge.
- Removed always-false `active_count < 2` guard in `_detect_regularization_component!`:
  both anchor particles are unconditionally seeded into the BFS, so the count is always ≥ 2.

### Changed

- `theory/RegularizedIntegrationDesign.md`: corrected and completed the spec with six
  additions — adaptive substep count formula, multi-substep ABA composition details,
  fixed-anchor hysteresis note, `active_steps` diagnostics field, `collision_bounce_radius`
  option section, and clarified `used_backend = :disabled` for zero-encounter enabled runs.
