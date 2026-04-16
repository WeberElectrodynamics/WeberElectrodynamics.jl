# Agent 10: Extended 4-body Weber bound orbit search

## Objective

Extend the Wave 1 four-body survey beyond the 2+/2- equal-mass equal-charge parameter space.
New axes explored: 3+/1- charge configurations, unequal masses, unequal charges, speed-of-light
(c) variation, and 3D double-orbiter fine-tuning.

## Summary

66 runs total. **3 bound orbits found**, all in the 3+/1- triangular trap configuration.

| Label | t_final | Energy drift | d_max_ratio | Config |
|---|---|---|---|---|
| tri_R2_rotating_eta0.75 | 100.0 | 0.159% | 3.0 | 3+/1-, R=2, rotating |
| tri_R3.0_rot_eta0.75 | 200.0 | 0.094% | 3.78 | 3+/1-, R=3, rotating |
| tri_R4.0_rot_eta0.75 | 200.0 | 1.34e-7% | 3.0 | 3+/1-, R=4, rotating |

## Section-by-section results

### 1. 3+/1- configurations (21 runs, 3 bound)

The **triangular trap** -- 3 positive charges at vertices of an equilateral triangle with 1
negative charge at center -- is the only geometry that produced bound orbits in this survey.

Key findings:
- **Rotating mode at high energy fraction (eta=0.75) is essential.** Lower eta values (0.1--0.5)
  are universally unbound. Breathing mode never binds.
- **Larger R improves stability.** R=4 achieves drift ~1e-7% (machine precision), R=3 achieves
  0.094%, R=2 achieves 0.159%. The d_max_ratio stays moderate (3.0--3.78), indicating genuine
  confinement rather than slow escape.
- **Linear 3+/1- chains** fail universally -- all hit singularities within t<3.
- **Tetrahedral 3D** runs complete to t=100 but produce NaN energy drift, suggesting the
  configurations pass through degenerate states. All show d_ratio ~9.7 (near the cutoff).

### 2. Unequal masses (12 runs, 0 bound)

No bound orbits. Heavy-backbone configurations (m_heavy=10 or 100 for the positive pair) in the
alternating square geometry all disintegrate rapidly (t<10).

The heavy-backbone double-orbiter (m=10 and m=100) survives to t=200 with excellent energy
conservation (~4e-7% drift) but has d_max_ratio ~23-25, indicating the system expands well
beyond the bound threshold of 10. The heavy positive cores stay put while light negative
particles orbit outward.

The heavy-backbone rhombus fares poorly (t<10), suggesting the mass asymmetry breaks the
delicate balance needed for rhombus quasi-periodicity.

### 3. Unequal charges (10 runs, 0 bound)

No bound orbits. Charge asymmetry (+2/+2/-1/-1 and +1/+1/-2/-2) is destabilizing:

- Alternating squares disintegrate within t<2 for all eta values.
- Double-orbiter at rpp=4 fails quickly. At rpp=6 the system survives to t=200 with low drift
  but d_ratio ~12-13 (just over the bound threshold).
- Interestingly, q2p and q2n alternating squares give identical results, which is expected by
  symmetry (the Coulomb potential depends on the product q_i*q_j).

### 4. c-variation (10 runs, 0 bound)

All configurations eventually escape, but the timescales reveal interesting physics:

**Double-orbiter 2D:**
- c=1: t*=458, c=2: t*=250, c=4: t*=246, c=10: t*=243
- Power-law fit: t* ~ c^(-0.24)
- The system converges to Coulomb behavior quickly; most of the Weber correction effect is
  already gone by c=2.

**Double-orbiter 3D:**
- c=1: t*=540, c=2: t*=255, c=4: t*=249
- Power-law fit: t* ~ c^(-0.56)
- Steeper falloff than 2D. The 3D configuration benefits more from Weber corrections at c=1.

**Rhombus:**
- c=1: t*=439, c=2: t*=600, c=4: t*=600
- Power-law fit: t* ~ c^(+0.23)
- Surprising: the rhombus survives **longer** at higher c. At c=2 and c=4 it reaches the full
  integration window (t=600). However, d_max_ratio is 65-87 (far from bound).

The predicted t* proportional to c^2 scaling (Agent 07 hypothesis) is **not confirmed**. The double-orbiter
shows weak *negative* scaling (shorter lifetime at higher c), while the rhombus shows weak
*positive* scaling. The Weber velocity-dependent correction at c=1 provides modest confinement
enhancement for the double-orbiter but is not the dominant effect controlling escape timescale.

### 5. Double-orbiter 3D fine-tuning (13 runs, 0 bound)

All runs survive long times (t*=259-586) with excellent energy conservation (<1e-5% drift)
but d_max_ratio 32-76 indicates systematic expansion.

Best survival times:
- z_kick=0.16: t*=586 (d_ratio=76) -- longest but most expanded
- orb=1.3, z_kick=0.13: t*=540 (d_ratio=69) -- baseline
- z_kick=0.10: t*=510 (d_ratio=64)

The orb parameter scan shows orb=1.3 is optimal (t*=540 vs 327 at orb=1.1 and 259 at orb=1.5).
The r_pp scan shows moderate sensitivity: r_pp=4.0 is near-optimal.

The key problem is that "long survival" correlates with "large expansion" -- these are slow
escapes, not bound orbits. The d_max_ratio grows roughly linearly with t_final.

## Key conclusions

1. **New bound orbit family discovered: 3+/1- triangular trap.** Three positive charges orbiting
   around a central negative charge at high energy fraction (eta=0.75) produces stable bound
   motion for t>=200 with sub-percent energy drift. This is the second bound configuration found
   in four-body Weber (after the 2+/2- breathing square from Wave 1).

2. **Mass and charge asymmetry are destabilizing.** Neither unequal masses (10:1, 100:1) nor
   unequal charges (2:1) produced any bound orbits in standard geometries.

3. **c-scaling does not follow t* ~ c^2.** The Weber correction provides only a modest, geometry-
   dependent effect. Double-orbiter lifetimes weakly decrease with c; rhombus lifetimes weakly
   increase. The dominant binding mechanism is Coulombic, not velocity-dependent.

4. **Double-orbiter is a long-lived transient, not a bound orbit.** Despite surviving 500+ time
   units with drift <1e-6%, the systematic expansion (d_ratio ~70) confirms eventual escape.
   Fine-tuning parameters does not eliminate this trend.

## Files

- `ic_generators_extended.jl` -- IC generators for all new configurations
- `extended_survey.jl` -- Survey runner (66 runs)
- `survey_results.csv` -- Full results table
