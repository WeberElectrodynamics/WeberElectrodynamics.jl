### Breaking changes

- `WeberProblem` no longer accepts flat `regularization_*` / `zollner_*` kwargs.
  Pass `RegularizationOptions(...)` and `ZollnerOptions(...)` structs directly:

  ```julia
  # Before (0.2.x)
  WeberProblem(...;
      regularization_enabled = true,
      regularization_backend = :adaptive_cartesian,
      zollner_enabled = true,
      zollner_a = 0.05,
  )

  # After (0.3.0)
  WeberProblem(...;
      regularization = RegularizationOptions(enabled = true, backend = :adaptive_cartesian),
      zollner        = ZollnerOptions(enabled = true, a = 0.05),
  )
  ```

  Both structs default to their disabled state (`RegularizationOptions()`,
  `ZollnerOptions()`), so problems with no regularization or Zöllner options
  require no changes.

### Changed

- `WeberProblem` constructor signature reduced from 20+ kwargs to 6 core kwargs
  plus the two options structs.
- Makie animation extension no longer performs an internal reverse-mapping from
  struct fields back to flat kwargs — it passes the structs through directly.
