# Conserved Quantities & Symmetries — 4-body 2+/2− Weber (2D)

All residuals computed numerically at 6 random phase-space points (seed 42).
A residual below 1e-10 certifies the identity to machine precision.

## Continuous generators — Poisson brackets with H

| F | definition | max\|{F,H}\| | conserved? |
|---|---|---|---|
| P_x | Σ p_{x,i} | ≈0 (4.4e-16) | **yes** |
| P_y | Σ p_{y,i} | ≈0 (1.9e-16) | **yes** |
| L | Σ (x_i p_{y,i} − y_i p_{x,i}) | ≈0 (1.6e-15) | **yes** |
| H | H itself (sanity) | ≈0 (0.0) | **yes** |
| D | Σ q_k p_k  (dilation) | 45.8 | no |
| d_x | Σ charge_i · x_i  (charge dipole, x) | 37.52 | no |
| d_y | Σ charge_i · y_i  (charge dipole, y) | 5.835 | no |

### Non-conserved residuals in closed form
- **{D,H} = 2T + U_c + 3 U_w** (derived from the q- and p-virial identities;
  verified numerically in virial.txt). D is conserved only on the zero-set
  of this expression, which is *not* a level set of H.
- **{d_x,H} = Σ charge_i · ∂H/∂p_{x,i} = Σ charge_i · (p_{x,i}/m_i + Weber-velocity terms)**
  — this is a charge-weighted current; non-zero in general.

## Discrete symmetries

| symmetry | action | residual max\|H∘g − H\| | invariant? |
|---|---|---|---|
| C (charge swap) | (1↔2),(3↔4) in q and p | ≈0 (1.1e-16) | **yes** |
| P (parity)      | q→−q, p→−p             | ≈0 (1.1e-16) | **yes** |
| T (time reverse)| p→−p, t→−t             | ≈0 (2.2e-16) | **yes** |

H depends on p only through p² (kinetic) and ṙ² (Weber), both even in p; so T
is manifest. P flips both r and v in each pair, leaving r² and ṙ² invariant.
C permutes identical-mass, identical-charge particles, leaving U pairwise invariant.

## Virial (Weber-modified)

    2⟨T⟩ + ⟨U_c⟩ + 3⟨U_w⟩ = 0   ⇔   2⟨T⟩ = −⟨U_c⟩ − 3⟨U_w⟩.
Reduces to 2⟨T⟩ = −⟨U⟩ when U_weber → 0.
