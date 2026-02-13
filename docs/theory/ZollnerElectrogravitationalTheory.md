# Zöllner Electrodynamic Theory of Matter and Gravitation (1872-1882)

## Scope

This document gives a self-contained theory formulation of Zöllner's electrogravitational program in a notation compatible with Weber electrodynamics. It is implementation-agnostic and focuses only on physical assumptions, mathematical structure, and resulting force laws.

## Core Postulates

Zöllner's program uses four linked assumptions:

1. Matter is composed of binary electrical constituents \((+e,-e)\) with inertial masses \((m_+,m_-)\).
2. The electrodynamic interaction of moving charges follows Weber's law.
3. For static interactions, there are three elementary potentials:
   - \(U_{++}\): repulsive (like-positive)
   - \(U_{--}\): repulsive (like-negative)
   - \(U_{+-}\): attractive (unlike)
4. The unlike attraction is slightly stronger than like repulsion by a small positive mismatch parameter \(a\).

The fourth postulate is the key move: gravitation appears as a residual after near-cancellation of large electrostatic terms in neutral matter.

## Static Interaction Model

Let the static pair potentials be

$$
U_{++}(r)=+\frac{k_{++}e^2}{r},\quad
U_{--}(r)=+\frac{k_{--}e^2}{r},\quad
U_{+-}(r)=-\frac{k_{+-}e^2}{r}.
$$

Zöllner's simplification (consistent with the derivation text) is

$$
k_{++}=k_{--}=k,\quad k_{+-}=k+a,\quad 0<a\ll k.
$$

So unlike attraction exceeds like repulsion by \(a\) at unit distance.

## Residual Attraction Between Two Neutral Dyads

Consider two neutral dyads \(A=(+_A,-_A)\) and \(B=(+_B,-_B)\), separated by distance \(R\), with internal dyad size negligible relative to \(R\).

Total static interaction:

$$
U_{AB}(R)=U(+_A,+_B)+U(-_A,-_B)+U(+_A,-_B)+U(-_A,+_B).
$$

Under the model above:

$$
U_{AB}(R)
=\frac{k e^2}{R}+\frac{k e^2}{R}-\frac{(k+a)e^2}{R}-\frac{(k+a)e^2}{R}
=-\frac{2ae^2}{R}.
$$

Hence

$$
F_{AB}(R)=-\frac{dU_{AB}}{dR}
=\frac{2ae^2}{R^2},
$$

an attractive inverse-square force. In units where \(e^2=1\), this is the historical \(2a\) at unit distance.

## Macroscopic Limit and Newton Form

Let body \(1\) contain \(N_1\) dyads and body \(2\) contain \(N_2\) dyads. With pairwise superposition and far-field approximation:

$$
U_{12}(R)=-\frac{2ae^2N_1N_2}{R},\quad
F_{12}(R)=\frac{2ae^2N_1N_2}{R^2}.
$$

If inertial mass scales with dyad count,

$$
M_i=\mu_i N_i,
$$

then

$$
F_{12}(R)=G_{12}\frac{M_1M_2}{R^2},\quad
G_{12}=\frac{2ae^2}{\mu_1\mu_2}.
$$

A composition-independent Newton constant requires an approximately universal \(\mu_i\approx\mu\), giving

$$
G=\frac{2ae^2}{\mu^2}.
$$

This is Zöllner's intended bridge from electrical microstructure to universal gravitation and equal free-fall acceleration.

## Relation to Weber Electrodynamics

For moving charges, Weber's pair force is

$$
\vec{F}^{(W)}_{ij}
=\frac{q_iq_j}{r_{ij}^2}\hat{r}_{ij}
\left(1-\frac{\dot r_{ij}^2}{2c^2}+\frac{r_{ij}\ddot r_{ij}}{c^2}\right).
$$

Using the same notation as `docs/theory/WeberElectrodynamics.md`, this is equivalently

$$
\vec{F}^{(W)}_{ij}
=\frac{q_iq_j}{r_{ij}^2}\hat{r}_{ij}
\left(1+\frac{1}{c^2}\left(\vec v\cdot\vec v+\vec r\cdot\vec a-\frac{3}{2}(\hat r\cdot\vec v)^2\right)\right).
$$

In Zöllner's gravitation argument, the Newton-like term is static and comes from the small mismatch \(a\), not from Weber's velocity/acceleration corrections. In this sense:

- electrostatic near-cancellation in neutral matter sets the gravitational residual;
- Weber dynamics governs electromagnetic and magnetic effects of moving charges.

## Effective n-Body Theory Form

At body scale, the induced gravitational potential has Newton form:

$$
U_g=-\sum_{i<j}\frac{G m_i m_j}{r_{ij}},
$$

with \(G\) interpreted as an emergent constant from \((a,e,\mu)\).

The combined theoretical Hamiltonian (theory-level decomposition) is

$$
H=\sum_i\frac{\|p_i\|^2}{2m_i}+U_W+U_g,
$$

where \(U_W\) is the Weber electrodynamic interaction and \(U_g\) is the Zöllner residual static attraction.

## Assumptions and Validity Domain

This formulation depends on:

1. Binary electrical constitution of matter.
2. Small static asymmetry \(a\) with \(a/k\ll1\).
3. Pairwise superposition of effective dyad interactions.
4. Far-field neutral-body approximation (internal structure scale \(\ell\ll R\)).
5. Approximate universality of mass-per-dyad ratio for composition-independent \(G\).

Within these assumptions, the model reproduces inverse-square attraction for neutral matter and provides a direct microphysical interpretation of Newton's constant.

## References

- Zöllner, J. K. F. *Über die Natur der Cometen: Beiträge zur Geschichte und Theorie der Erkenntnis*. Engelmann, Leipzig (1872; 3rd ed. 1883). Digital facsimile: <https://archive.org/details/berdienaturdeco00zolgoog>
- Zöllner, J. K. F. "Über die physikalischen Beziehungen zwischen hydrodynamischen und elektrodynamischen Erscheinungen." *Annalen der Physik und Chemie* 158 (1876), 497-539. Digital volume facsimile: <https://archive.org/details/annalenderphysi97unkngoog>
- Zöllner, J. K. F. *Principien einer elektrodynamischen Theorie der Materie*. Wilhelm Engelmann, Leipzig (1876). Digital facsimile: <https://archive.org/details/principieneinere00zgoog>
- Zöllner, J. K. F. *Erklärung der universellen Gravitation aus den statischen Wirkungen der Elektrizität und die allgemeine Bedeutung des Weberschen Gesetzes*. L. Staackmann, Leipzig (1882). Digital facsimile: <https://archive.org/details/erklrungderuniv00zgoog>
- Zöllner, J. C. F. *Wissenschaftliche Abhandlungen*, no. 2 (contains section "Ueber die Ableitung der Newton'schen Gravitation aus den statischen Wirkungen der Elektricität", pp. 417-449 in volume pagination). Facsimile used for derivation check: <https://archive.org/details/wissenschaftlic01unkngoog>
