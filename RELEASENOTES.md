<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Changed
- Plot labels in the `Plots` extension now use `LaTeXStrings` for axis labels
  and legend entries across all 11 plot functions (energy, pair energy,
  energy errors, pair forces, phase space, trajectories, momentum errors,
  and the four Zöllner plots). Titles remain prose. README figures
  regenerated from `examples/two_body_reference.ipynb`.
- `plot_energy` is now a single-panel plot of the energy components
  (`T`, `U`, `H`); the redundant relative-error panel has been removed.
- `plot_energy_errors` is now a two-panel plot (local error and relative
  energy error `|ΔE/E₀|`). The percentage scaling has been dropped in favour
  of the raw fraction, and the Hamiltonian-validation panel has been removed.
- `plot_momentum_errors` titles are now plain prose (`Linear Momentum Drift`
  / `Angular Momentum Drift`); the `‖ΔP‖` / `‖ΔL‖` symbols have moved from
  the panel titles to the y-axis labels where they belong.
- Absolute-value and norm labels across all plot functions
  (`plot_energy_errors`, `plot_momentum_errors`, `plot_pair_forces`,
  `plot_zollner_force_decomposition`) now use plain `|…|` / `||…||` bars
  instead of `\lvert…\rvert` / `\lVert…\rVert`, which the Plots.jl GR
  backend silently drops, leaving labels blank.
- README quickstart gained a third figure — `plot_energy_errors` — and the
  `examples/two_body_reference.ipynb` notebook now writes the three
  README figures to `examples/figures/` on every run.

### Added
- `LaTeXStrings` as a weak dependency of the `WeberElectrodynamicsPlotsExt`
  extension (compat `1`).
