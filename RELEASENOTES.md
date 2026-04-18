<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Changed

- `animate_weber` visual polish pass. Trajectories now use Wong-palette
  colors with an alpha-graded tail that fades older segments, rounded
  line joins, and stroked particle markers. 3D panels gain perspective
  (`perspectiveness = 0.3`) and a tight `viewmode = :fit`. Phase-space
  sidebar and info labels are restyled for readability. A shared
  `_weber_theme()` (LaTeX fonts, muted grid, clean spines) is applied
  via `with_theme` so all axes render consistently. No API change.
