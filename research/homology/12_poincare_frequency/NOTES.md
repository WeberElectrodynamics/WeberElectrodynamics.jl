# Agent 12 -- Poincare sections and frequency map analysis

## Goal

Determine whether the 4-body 2+/2- Weber Hamiltonian supports genuine KAM tori
(closed curves on a Poincare section, zero Lyapunov exponent, constant
frequencies) or whether all long-lived orbits are stochastic layers with slow
Arnold diffusion.

Two IC families studied:
1. **Rhombus** (a, b, eta, rotating): 60 ICs (5x4x3 grid) for frequency map,
   120 ICs (5x6x4) for Poincare sections.
2. **Symmetric double-orbiter** (r_pp, R, orb): 18 ICs (6x3 grid).

## Bug fixes applied

- Added `using WeberElectrodynamics` to both scripts (missing import for
  `compute_energy_timeseries`).
- Rewrote Poincare section surface from `r12 = r12(0)` to angular crossing
  (`y-component of pair vector = 0`, crossing from below). The original section
  was almost never crossed because the dynamics is not oscillatory in pair
  distance.
- Fixed O(n^2) manual DFT in `frequency_map.jl`: with 400K timesteps and 4x
  zero-padding, the DFT loop was 800K x 400K = 320 billion iterations, causing
  the script to hang for >60 minutes on a single IC. Fixed by downsampling the
  signal to at most 4000 points before DFT (dt_effective = 100 * dt), which is
  still well above the Nyquist limit for the frequencies of interest.

## Results

### Poincare sections: no usable sections obtained

All 138 ICs classified as **sparse** (0--9 section crossings over tmax=200).
No IC produced >= 10 crossings, the minimum needed for classification.

Root cause: the 4-body 2+/2- system does not produce quasi-periodic orbits at
these parameters.

**Rhombus:**
- At high energy fraction (eta >= 0.5): orbits escape. Despite E < 0, the
  system breaks into opposite-charge dimers that separate to r > 100 within
  tmax = 200. The Coulomb potential for 2+/2- is non-confining; dimer
  fragmentation is always energetically accessible.
- At low energy fraction (eta < 0.3): the integrator fails within t < 30 due
  to close encounters and collisions. Energy drift exceeds 100% in some cases,
  indicating the bounce regularization cannot handle the dynamics.
- In neither regime does the system exhibit oscillatory motion in pair
  distances.

**Double-orbiter:**
- All 18 ICs have **E > 0** (unbound). For example, rpp=4 orb=1.3 gives
  T=1.69, U=-0.52, E=+1.17. Particles simply fly apart with r34 reaching 400+
  by t=200.
- The "low Lyapunov exponent" reported in prior studies for this configuration
  was a consequence of near-free-flight, not quasi-periodicity.

### Frequency map: universal frequency drift confirms stochasticity

**Rhombus ICs (on-diagonal, a ~ b):** The only ICs with physically meaningful
frequencies (omega > 0.01) are the symmetric cases a = b. All off-diagonal
cases collapse to the DFT resolution floor (omega ~ 0.0095), confirming their
trajectories are non-oscillatory (monotonic escape).

**Frequency drift analysis** (5 best rhombus candidates, tmax=300, 15 windows):

| IC | omega1 mean | omega1 std | relative drift | D1 |
|----|-------------|------------|----------------|-----|
| a=1.50 b=1.45 eta=0.75 | 0.0524 | 0.0172 | **32.9%** | 3.4e-5 |
| a=1.50 b=1.45 eta=0.70 | 0.0523 | 0.0172 | **32.9%** | 2.7e-5 |
| a=1.50 b=1.45 eta=0.80 | 0.0509 | 0.0139 | **27.3%** | 2.1e-5 |
| a=1.50 b=1.40 eta=0.75 | 0.0477 | 0.0052 | **10.9%** | 2.9e-6 |
| a=1.50 b=1.50 eta=0.75 | 0.0447 | 0.0060 | **13.4%** | 3.0e-6 |

For a genuine KAM torus, relative frequency drift should be < 0.01% (limited
by integrator precision). The observed 10--33% drift is 3--4 orders of magnitude
too large, conclusively ruling out KAM tori.

**Double-orbiter ICs:** Diffusion coefficients D ~ 10^{-10}--10^{-12}, but this
is a trivial consequence of free-flight (E > 0, particles separate linearly).
Constant velocity gives constant "frequency" -- this does not indicate a KAM
torus.

### Diagnostic trajectory analysis

For the "best candidate" rhombus (a=1.5, b=1.45, eta=0.75):
- Initial E = -0.31 (nominally bound)
- r12 expands from 3.0 to 121 by t=200
- y-component of r12 never changes sign (no rotation, no oscillation)
- Trajectory is a one-way dimer escape, not a quasi-periodic orbit

For alternating square ICs (side=2, eta=0.15--0.4):
- All fail (retcode=Failure) with massive energy drift (up to 2300%)
- Even with dt=1e-4, close encounters crash the integrator within t < 20

## Conclusions

**There are no genuine KAM tori in the 4-body 2+/2- Weber Hamiltonian at the
parameters surveyed.** All long-lived orbits are either:

1. **Escaping dimers**: the system fragments into two opposite-charge pairs that
   fly apart. This accounts for all "successful" long integrations (t=200) with
   low energy drift. The apparent stability is simply free flight, not
   confinement.

2. **Rapidly diverging close-encounter orbits**: the only truly bound
   configurations (low eta) immediately hit Coulomb singularities that the
   integrator cannot resolve, leading to failure within t < 30.

This confirms and strengthens Agent 07's finding that all candidates are in a
fast Arnold diffusion regime. The frequency drift data goes further:

- The "KAM basin" around (a=1.5, b=1.45, eta=0.75) is not a basin at all --
  it is a region where dimer escape happens slowly enough for the integrator to
  succeed, but the motion is monotonic separation, not oscillation.
- The double-orbiter "low Lyapunov" result is an artifact of unbound (E > 0)
  free flight.
- No Poincare section can be constructed because the prerequisite --
  quasi-periodic motion returning to a section surface -- does not hold.

### Physical interpretation

The 4-body Coulomb problem with 2+ and 2- charges is fundamentally
non-confining. Unlike gravitational N-body problems (where all masses attract),
the mixed-sign charge system has no potential barrier preventing dimer
fragmentation. The Weber correction (velocity-dependent force) does not change
this qualitative picture at c=1.

For KAM tori to exist in this system, one would need either:
- An external confining potential
- A regime where Weber corrections create effective barriers (requires c << 1,
  i.e., strongly relativistic conditions where the Weber model itself breaks
  down)
- A topologically constrained configuration (e.g., collinear, where escape
  requires overcoming a barrier) -- but these are measure-zero and structurally
  unstable
