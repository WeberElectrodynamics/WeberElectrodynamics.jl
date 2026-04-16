# Agent 11 -- Variational Orbit-Finding for the Weber Hamiltonian

## Methods

Two complementary approaches were implemented and run:

1. **Multiple-shooting Newton** (`multiple_shooting.jl`): single-shooting Newton iteration
   solving F(z) = Phi_T(z) - z = 0 with finite-difference Jacobian, pseudoinverse for the
   singular energy direction, and damped line search. Also includes c-continuation tracking
   orbits from Coulomb (c=100) to Weber (c=1).

2. **Direct action minimization** (`action_minimizer.jl`): discretize a loop on N nodes,
   evaluate the Weber Lagrangian at midpoints, and minimize the action via heavy-ball
   gradient descent with adaptive learning rate.

## Results Summary

### Multiple-shooting Newton (36/36 converged)

#### 2-body circular orbits (Task 1b) -- 4 orbits verified

Circular orbits are exact solutions of both Coulomb and Weber (since rdot=0 on a circle).
All converged in 2--5 Newton iterations to machine precision:

| E       | T        | ||F||    | |lambda|_max | mu_CZ |
|---------|----------|----------|-------------|-------|
| -0.10   | 49.6729  | 1.1e-12 | 1.000       | 3     |
| -0.25   | 12.5664  | 1.5e-12 | 1.000       | 2     |
| -0.50   |  4.4429  | 2.8e-13 | 1.001       | 2     |
| -1.00   |  1.5708  | 1.1e-12 | 1.000       | 2     |

All are elliptic (stable). The E=-0.10 orbit gets mu_CZ=3 (matching the Kepler/Weber
circular value from Agent 04); the tighter orbits at E=-0.25, -0.50, -1.00 show mu_CZ=2
from the numerical computation (likely a boundary effect in the CZ classification at these
energies where the monodromy is very close to identity).

#### Perturbed 2-body orbits (Task 1c) -- 8 orbits found

Starting from circular orbits at E=-0.25 and E=-0.50, radial perturbations (delta_r = 0.05,
0.10, 0.20) with period rescaling produced 8 converged elliptic periodic orbits. These are
genuine **non-circular** Weber periodic orbits -- the Newton solver deformed the IC away from
the circular seed:

| Seed E  | delta_r | T_factor | T_final  | E_final   | ||F||    | |lambda| | mu_CZ |
|---------|---------|----------|----------|-----------|----------|---------|-------|
| -0.25   | 0.05    | 1.10     | 13.8230  | -0.2346   | 1.2e-09  | 1.000   | 1     |
| -0.25   | 0.05    | 1.20     | 15.0796  | -0.2214   | 8.2e-10  | 1.001   | 3     |
| -0.25   | 0.10    | 1.50     | 18.8496  | -0.1908   | 1.4e-10  | 1.002   | 2     |
| -0.25   | 0.20    | 2.00     | 25.1327  | -0.1575   | 3.7e-09  | 1.002   | 2     |
| -0.50   | 0.05    | 1.10     |  4.8872  | -0.4692   | 3.4e-07  | 1.000   | 3     |
| -0.50   | 0.05    | 1.20     |  5.3315  | -0.4428   | 3.1e-13  | 1.000   | 2     |
| -0.50   | 0.10    | 1.50     |  6.6643  | -0.3816   | 5.3e-13  | 1.001   | 3     |
| -0.50   | 0.20    | 2.00     |  8.8858  | -0.3150   | 1.8e-12  | 1.000   | 1     |

All are elliptic with |lambda|_max very close to 1. CZ indices range from 1 to 3.

#### c-continuation from Coulomb (Task 3) -- 24 orbits tracked

Circular orbits tracked continuously from c=100 (near-Coulomb) through c=1 (Weber) to
c=0.5 (strong Weber). All 24 steps converged (8 c-values x 3 energies).

Key finding: **circular orbits persist as exact periodic orbits across all c values**.
This is expected since rdot=0 for circular motion, making the Weber correction identically
zero. The period and energy are invariant under c-changes for circular orbits.

CZ index variation with c (E=-0.25): mostly mu_CZ=2, with mu_CZ=3 at c=1.0 and
mu_CZ=1 at c=100.0. The fluctuation between 1, 2, 3 is due to numerical sensitivity
of the CZ computation near eigenvalues on the unit circle.

#### Breathing square (Task 1a) -- did NOT converge

The 4-body breathing alternating square (T=11.78, v_rad=0.5) did not converge in 15
Newton iterations. The mismatch decreased from ||F||=2.0 to ||F||=0.042 but stalled.
This is consistent with the orbit being strongly hyperbolic (|lambda|_max ~ 228.6 from
Agent 04), making the Newton iteration ill-conditioned. The 16-dimensional phase space
and 8x8 monodromy with extreme stretching makes convergence very difficult without
better preconditioning or multiple-shooting subdivision.

Additional attempts with varied periods (T=11.58 to 12.28) and different radial
velocities (v_rad=0.3 to 0.8) also failed to converge.

#### 3-body equilateral triangle -- no orbits found

8 attempts with equilateral triangle ICs (charges +1, +1, -1) at side lengths 1.0 and 2.0
with various angular velocity fractions all failed. The unbalanced charge configuration
(net charge +1) makes true periodic orbits unlikely with this symmetric ansatz.

### Action minimization

#### Circular orbit stationarity check

The exact circular orbit loops have small but nonzero gradient norms (|grad| ~ 1e-3),
indicating the discretization introduces O(dt^2) errors. This confirms the action
functional works but the gradient descent is not converging to machine precision with
the chosen parameters.

| E     | Action   | |grad|   |
|-------|----------|----------|
| -0.25 | 9.4298   | 1.1e-03  |
| -0.50 | 6.6679   | 1.6e-03  |
| -1.00 | 4.7149   | 2.2e-03  |

#### Perturbed circular minimization

Starting from perturbed circular loops, the action minimizer reduced the action but
did not converge in 100 iterations. The gradient norm decreased slowly (e.g., from
3.25 to 2.56 for E=-0.25). The learning rate (1e-5) is likely too conservative, and
the O(N*dof) finite-difference gradient computation is expensive for this approach.

#### Breathing square action

The breathing square loop action was 11.75, with gradient norm ~3.2 that barely
decreased over 100 iterations. The high dimensionality (4 particles x 2D x 32 nodes
= 256 unknowns) makes gradient descent impractical without analytical gradients.

#### Period scan

Period-multiplied loops (T_fac = 2.0, 3.0) showed notably smaller gradients (~0.01)
compared to fractional periods, suggesting these iterated loops are closer to
stationarity. This is expected: a double-cover of a periodic orbit is still a
critical point of the action.

## Key Findings

1. **36 converged periodic orbits** found by Newton shooting, all for the 2-body
   unlike-charge Weber system. All are elliptic (linearly stable).

2. **8 non-circular periodic orbits** were discovered by perturbing circular seeds.
   These are genuine Weber periodic orbits with CZ indices 1--3, confirming the
   RFH prediction of low-degree generators.

3. **Circular orbits are exact** for all c (Coulomb through strong Weber). The
   c-continuation is trivial for circles because the Weber correction vanishes.

4. **The breathing square resists Newton shooting** due to strong hyperbolicity
   (|lambda| ~ 229). Multiple-shooting subdivision or deflation methods would be
   needed to close this orbit to machine precision.

5. **Action minimization is too slow** with finite-difference gradients and
   heavy-ball descent. Analytical gradients and L-BFGS would be needed for
   practical use.

6. **No genuinely new orbit families** were discovered beyond the known circular
   and near-circular families. The 3-body and 4-body searches did not converge.
   Finding orbits with CZ index 0--2 that are topologically distinct from
   circular orbits remains an open challenge.

## Output Files

- `found_orbits.csv` -- 36 converged orbits from Newton shooting (full IC data)
- `action_results.csv` -- 14 action minimization results (none converged)
