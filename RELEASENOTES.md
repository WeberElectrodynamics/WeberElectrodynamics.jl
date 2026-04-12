<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Added
- New `WeberElectrodynamicsPythonPlotExt` weak-dependency extension — activate
  with `using PythonPlot` for publication-quality matplotlib figures as an
  alternative to the existing `using Plots` backend. Mirrors the full
  11-function plotting surface (`plot_trajectories`, `plot_energy`,
  `plot_pair_energy`, `plot_energy_errors`, `plot_pair_forces`,
  `plot_phase_space`, `plot_momentum`, and the four Zöllner variants), applies
  a global `rcParams` publication preset (serif font, Computer Modern mathtext,
  inward ticks, minor ticks, dashed grid, no top/right spines), and exposes
  `set_matplotlib_style!(; usetex=false)` for reapplying the preset or opting
  into LaTeX rendering. Load exactly one backend per session — the two
  extensions share the same stub functions and the last-loaded backend wins.
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
