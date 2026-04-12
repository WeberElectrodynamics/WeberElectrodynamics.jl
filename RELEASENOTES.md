<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

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
