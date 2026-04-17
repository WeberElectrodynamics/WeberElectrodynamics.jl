<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Added
- Research notebook `exhaust_nonzero_radial_velocity_ics.ipynb` exercising
  every IC recipe from `theory/NonZeroRadialVelocityBoundICs.md` across six
  scenarios: two-body mid-flight (2D outbound, 2D inbound, 3D tilted), hot
  binary with cold orbiter (3D), breathing square (N=4), and breathing
  hexagon (N=6). Each scenario builds physical ICs, applies the §6 forward
  map F to canonical momenta, runs the §7 verification checklist, integrates
  for ≥1 natural period, and produces trajectory + energy-error plots. The
  `verify_ic` helper amends the doc's §6.4 boxed identity with the exact
  analytic F-consistency check and also reports the Wesley-vs-Legendre H
  mismatch.
- Research notebook `scenario_3_makie_live.ipynb` demonstrating the
  redesigned `animate_weber` dashboard on the scenario-3 bound non-zero-ṙ
  IC (3D trajectory panel, phase-space sidebar, live energy-error readout).

### Changed
- `animate_weber` dashboard redesigned: the trajectory panel is now the
  dominant view, with the phase-space panel as a sidebar and a compact live
  energy-error readout below. The separate kinetic/potential,
  linear-momentum, and angular-momentum panels have been removed to reduce
  visual clutter and free up room for the trajectory view.
- The trajectory panel now uses `Axis3` when `prob.system.dims == 3`, giving
  a proper 3D view of the simulation alongside the existing 2D `Axis`
  rendering.
- Default `figure_size` reduced from `(1400, 900)` to `(1200, 800)` to
  better suit typical laptop displays.

The `animate_weber(prob)` / `animate_weber(sol)` public API (function names,
keyword arguments, return value) is unchanged.
