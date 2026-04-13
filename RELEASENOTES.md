<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Breaking changes
- `plot_momentum` removed. Replaced by `plot_momentum_errors(data::MomentumData)`,
  which plots conservation errors on a single panel: absolute drift
  `‖P(t) − P(0)‖` and `|L(t) − L(0)|` (log scale), with relative drift
  overlaid on a secondary axis where the initial magnitude is nonzero. The
  old raw-timeseries plot was visually flat for symplectic integrators in
  the COM frame and conveyed no integrator-quality signal. Migration:
  replace `plot_momentum(mom)` with `plot_momentum_errors(mom)`.

### Added
- Aqua.jl quality gate in the test suite — checks stale deps, unbound type
  parameters, undefined exports, compat bounds, and Project.toml formatting.
- CI matrix now includes `macos-latest` and `windows-latest` on Julia 1 in
  addition to the Ubuntu × {1.9, 1} coverage.
- `codecov.yml` with `project: auto` and `patch: 80%` targets, so unrelated
  PRs don't trip spurious coverage-drop failures.

### Changed
- `Project.toml` gains compat bounds for stdlibs (`LinearAlgebra`, `Printf`,
  `Random`, `Test`) and for the `Aqua` test extra, required by Aqua's
  `deps_compat` check.

### Chore
- Ignore `.claude/` (editor scratch directory) in `.gitignore`.
