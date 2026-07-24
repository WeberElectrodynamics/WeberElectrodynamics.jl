### Breaking changes

- Removed the Zöllner electrogravitational extension outright, with no
  deprecations or compatibility aliases. `ZollnerOptions`, `zollner_term`,
  the `zollner=` problem keyword, all κ accessors and fields, κ-aware
  initial-condition keywords, Zöllner statistics fields, and the four
  Zöllner-specific plotting functions no longer exist.
- Compiled equations now use
  `dq_dt_compiled(out, q, p, t, params)` and
  `dp_dt_compiled(out, q, p, t, params)`; compiled Hamiltonians use
  `hamiltonian_compiled(q, p, t, params)`. Pair-decomposition closures now
  use `(i, j, q, p, params)`. Custom coefficients belong in
  `param_symbols`.
- Default systems now expose only the `:weber` named term. Standard pair
  coupling is directly `qᵢqⱼ`.
- JLD2 solution archives now use Weber-only format v2. Format-v1 archives
  from 0.5.x are rejected and must be exported or regenerated before
  upgrading.

### Changed

- Simplified system construction, problems, solvers, regularization buffers,
  statistics, animation recycling, and archives around the single packed
  parameter vector `[masses; charges; c]`.
- Moved the former Zöllner theory note into the non-authoritative `_research`
  sandbox and removed it from published documentation.
- Replaced the removed-model regression fixture with a pure-Weber 2D
  adaptive-Cartesian fixture. Existing pure-Weber trajectories remain
  unchanged to `1e-12`.

### Known issues

- The canonical-momentum/Hamiltonian correctness issue documented in
  `theory/HamiltonianCorrectness.md` remains confirmed and is **not fixed** by
  this release. Removing the secondary coupling model only narrows the future
  remediation to pure Weber dynamics.
