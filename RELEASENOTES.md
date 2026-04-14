<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Changed
- Plot labels in the `Plots` extension now use `LaTeXStrings` for axis labels
  and legend entries across all 11 plot functions (energy, pair energy,
  energy errors, pair forces, phase space, trajectories, momentum errors,
  and the four Zöllner plots). Titles remain prose. README figures
  regenerated from `examples/two_body_reference.ipynb`.

### Added
- `LaTeXStrings` as a weak dependency of the `WeberElectrodynamicsPlotsExt`
  extension (compat `1`).
