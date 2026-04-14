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
  Absolute-value y-axis labels now use `\left\lvert…\right\rvert` so they
  render correctly in Jupyter/KaTeX.

### Added
- `LaTeXStrings` as a weak dependency of the `WeberElectrodynamicsPlotsExt`
  extension (compat `1`).
