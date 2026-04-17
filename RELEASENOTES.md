<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Changed
- `animate_weber` dashboard redesigned: the trajectory panel is now the dominant
  view, with the phase-space panel as a sidebar and a compact live energy-error
  readout below. The separate kinetic/potential, linear-momentum, and
  angular-momentum panels have been removed to reduce visual clutter and free
  up room for the trajectory view.
- The trajectory panel now uses `Axis3` when `prob.system.dims == 3`, giving a
  proper 3D view of the simulation alongside the existing 2D `Axis` rendering.
- Default `figure_size` reduced from `(1400, 900)` to `(1200, 800)` to better
  suit typical laptop displays.

The `animate_weber(prob)` / `animate_weber(sol)` public API (function names,
keyword arguments, return value) is unchanged.
