# Weber Electrodynamics

## Constants and units

We use absolute system of units throughout this document, unless otherwise specified.

|              | Description                        | SI Unit                             | Gauss-Weber absolute units         |
| ------------ | ---------------------------------- | ----------------------------------- | ---------------------------------- |
| $q^2$        | electric charge squared            | $C^2$                               | $\frac{mg \cdot mm^3}{s^2}$        |
| $\epsilon_0$ | permittivity of free space         | $\frac{C^2\cdot s^2}{kg \cdot m^3}$ | dimensionless                      |
| $c$          | speed of light                     | $\frac{\text{m}}{\text{s}}$         | $\frac{mm}{s}$                     |
| $q_i$        | electric charge of particle i      | $C$                                 | $\sqrt{\frac{mg \cdot mm^3}{s^2}}$ |
| $m_i$        | mass of particle i                 | $kg$                                | $mg$                               |

## 2-Particle system notation

### Mass and Charge of Particle 1 and Particle 2

$$
m_1, q_1
$$

$$
m_2, q_2
$$

### Positions for Particle 1 and Particle 2

#### Position

$$
x_1(t), y_1(t), z_1(t)
$$

$$
x_2(t), y_2(t), z_2(t)
$$

#### Velocity

$$
\dot{x}_1(t), \dot{y}_1(t), \dot{z}_1(t)
$$

$$
\dot{x}_2(t), \dot{y}_2(t), \dot{z}_2(t)
$$

#### Acceleration

$$
\ddot{x}_1(t), \ddot{y}_1(t), \ddot{z}_1(t)
$$

$$
\ddot{x}_2(t), \ddot{y}_2(t), \ddot{z}_2(t)
$$

### Vector notation for Particle 1 and Particle 2

#### Position

$$
\vec{r}_1 = (x_1(t), y_1(t), z_1(t)), \quad \vec{r}_2 = (x_2(t), y_2(t), z_2(t))
$$

#### Velocity

$$
\vec{v}_1 = (\dot{x}_1(t), \dot{y}_1(t), \dot{z}_1(t)), \quad \vec{v}_2 = (\dot{x}_2(t), \dot{y}_2(t), \dot{z}_2(t))
$$

#### Acceleration

$$
\vec{a}_1 = (\ddot{x}_1(t), \ddot{y}_1(t), \ddot{z}_1(t)), \quad \vec{a}_2 = (\ddot{x}_2(t), \ddot{y}_2(t), \ddot{z}_2(t))
$$

### Momenta for Particle 1 and Particle 2

The Weber Lagrangian defined below depends on velocities through both the
kinetic energy and the relative radial velocity $\dot{r}$, so the canonical
momentum $\vec{p}_i = \partial L/\partial \vec{v}_i$ is **not** the kinetic
momentum $m_i\vec{v}_i$. Writing the Weber momentum correction as

$$
\vec{\alpha} = \frac{q_1 q_2}{c^2}\frac{\dot{r}\,(\vec{r}_1 - \vec{r}_2)}{r^2}
= \frac{q_1 q_2}{r c^2}\,\dot{r}\,\hat{r},
$$

where $q_1$ and $q_2$ are the charges, $c$ is the speed of light, $r$ is the
separation, $\hat{r}$ is the unit separation vector, and $\dot{r}$ is the
physical relative radial velocity, the canonical momenta are

$$
\vec{p}_1 = m_1\vec{v}_1 - \vec{\alpha},
\qquad
\vec{p}_2 = m_2\vec{v}_2 + \vec{\alpha}.
$$

Here $\vec{p}_1$ and $\vec{p}_2$ are canonical momenta, $m_1\vec{v}_1$ and
$m_2\vec{v}_2$ are kinetic momenta. The opposite signs of $\vec{\alpha}$
preserve total momentum, $\vec{p}_1 + \vec{p}_2 = m_1\vec{v}_1 + m_2\vec{v}_2$.

The correction is purely radial, so only the radial component is affected.
Solving for the physical velocities gives

$$
\vec{v}_1 = \frac{\vec{p}_1 + \vec{\alpha}}{m_1},
\qquad
\vec{v}_2 = \frac{\vec{p}_2 - \vec{\alpha}}{m_2}.
$$

These are implicit because $\vec{\alpha}$ itself contains $\dot{r}$. For two
particles the implicit relation collapses to a single scalar equation. With
$\mu = m_1m_2/(m_1+m_2)$ the reduced mass and

$$
p_r = \mu\,\hat{r}\mathbin{\cdot}
\left(\frac{\vec{p}_1}{m_1} - \frac{\vec{p}_2}{m_2}\right)
$$

the canonical momentum conjugate to the relative radial coordinate, we have

$$
p_r = \left(\mu - \frac{q_1 q_2}{r c^2}\right)\dot{r},
\qquad
\dot{r} = \frac{p_r}{\mu - q_1 q_2/(r c^2)}.
$$

Here the second form is valid whenever its denominator is non-zero. For like
charges the denominator vanishes at Weber's critical radius
$\rho = q_1 q_2/(\mu c^2)$, where a canonical momentum no longer determines a
finite physical velocity. Below $\rho$ the effective radial inertia is negative
but finite.

Canonical and kinetic momenta coincide exactly when $\dot{r} = 0$, and
asymptotically as $c \to \infty$.

### Relative Coordinates and Vectors for Particle 1 and Particle 2

#### Relative position

$$
x = x_1 - x_2
$$

$$
y = y_1 - y_2
$$

$$
z = z_1 - z_2
$$

#### Relative velocity

$$
\dot{x} = \dot{x}_1 - \dot{x}_2
$$

$$
\dot{y} = \dot{y}_1 - \dot{y}_2
$$

$$
\dot{z} = \dot{z}_1 - \dot{z}_2
$$

#### Relative acceleration

$$
\ddot{x} = \ddot{x}_1 - \ddot{x}_2
$$

$$
\ddot{y} = \ddot{y}_1 - \ddot{y}_2
$$

$$
\ddot{z} = \ddot{z}_1 - \ddot{z}_2
$$

#### Relative vectors

$$
\vec{r} = (x, y, z)
$$

$$
\vec{v} = (\dot{x}, \dot{y}, \dot{z})
$$

$$
\vec{a} = (\ddot{x}, \ddot{y}, \ddot{z})
$$

### Relative distance

$$
r = |\vec{r}| = \sqrt{x^2 + y^2 + z^2}
$$

### Relative radial velocity

$$
\dot{r} = \frac{x\dot{x} + y\dot{y} + z\dot{z}}{r} = \frac{\vec{r} \cdot \vec{v}}{r}
$$

### Relative radial acceleration

$$
\begin{aligned}
\ddot{r} &= \frac{x\ddot{x} + y\ddot{y} + z\ddot{z}}{r} + \frac{\dot{x}^2 + \dot{y}^2 + \dot{z}^2}{r} - \frac{(x\dot{x} + y\dot{y} + z\dot{z})^2}{r^3} \\
&= \frac{\vec{r} \cdot \vec{a}}{r} + \frac{\vec{v} \cdot \vec{v}}{r} - \frac{(\vec{r} \cdot \vec{v})^2}{r^3} \\
&= \frac{1}{r}\left(\vec{r} \cdot \vec{a} + \vec{v} \cdot \vec{v} - (\hat{r} \cdot \vec{v})^2\right)
\end{aligned}
$$

### Unit relative position vector

$$
\hat{r} = \frac{\vec{r}}{r}
$$

## Weber Potential Energy

For a two-particle system, we denote the distance between particles 1 and 2 as $r_{12}$. To keep expressions clean, we write $r = r_{12}$, $\dot{r} = \dot{r}_{12}$, $\ddot{r} = \ddot{r}_{12}$ when the context is clear.

In SI units:

$$
U_{SI} = \frac{1}{4\pi\epsilon_0} \frac{q_1 q_2}{r}\left(1 - \frac{\dot{r}^2}{2c^2}\right)
$$

In absolute units:

$$
U = \frac{q_1 q_2}{r}\left(1 - \frac{\dot{r}^2}{2c^2}\right)
$$

## Weber Force

The force can be derived from the potential in two ways. The first assumes $r$ is an explicit function of time. The second assumes $r$ is only a function of the coordinates.

### Explicit Time Parameterization of $r$ and $\dot{r}$

If $r$ and $\dot{r}$ are both functions of time, that is, $r = r(t)$ and $\dot{r} = \dot{r}(t)$, then the force on particle 1 due to particle 2 is:

$$
\vec{F}_{12} = -\frac{dU}{dr}\hat{r}
$$

Using the chain rule:

$$
\frac{d\dot{r}^2}{dr} = 2\dot{r} \frac{d\dot{r}}{dr} = 2\dot{r} \frac{d\dot{r}/dt}{dr/dt} = 2\dot{r} \frac{\ddot{r}}{\dot{r}} = 2\ddot{r}
$$

Computing the derivative:

$$
\frac{dU}{dr} = -\frac{q_1 q_2}{r^2}\left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right)
$$

This yields:

$$
\vec{F}_{12} = \frac{q_1 q_2}{r^2}\hat{r}\left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right)
$$

### Implicit coordinate-based formulation of $r$ and $\dot{r}$

If we take $\nabla_1$ as the gradient operator with respect to the coordinates of particle 1, that is $\nabla_1 = \left(\frac{\partial}{\partial x_1}, \frac{\partial}{\partial y_1}, \frac{\partial}{\partial z_1}\right)$, then we can derive Weber's force as follows:

$$
\vec{F}_{12} = -\nabla_1 U
$$

where the first component of the force is given by:

$$
-\frac{\partial U}{\partial x_1} = q_1 q_2 \frac{(x_1 - x_2)}{r^3} \left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right)
$$

Weber's force can be rewritten in vector form as:

$$
\vec{F}_{12} = \frac{q_1 q_2}{r^2}\hat{r}\left(1 + \frac{1}{c^2} \left(\vec{v} \cdot \vec{v} + \vec{r} \cdot \vec{a} - \frac{3}{2} (\hat{r} \cdot \vec{v})^2 \right)\right)
$$

## Lagrangian formulation

$$
T = \frac{1}{2} m_1 (\vec{v}_1 \cdot \vec{v}_1) + \frac{1}{2} m_2 (\vec{v}_2 \cdot \vec{v}_2)
$$

We define $S$ as the velocity-dependent interaction term with a positive sign in the velocity-dependent term (note this differs from the potential energy $U$ which has a negative sign). Again, we write $r = r_{12}$ and $\dot{r} = \dot{r}_{12}$ for cleaner notation:

$$
S = \frac{q_1 q_2}{r}\left(1 + \frac{\dot{r}^2}{2c^2}\right)
$$

$$
L = T - S
$$

In the following, we use $q_i$ to denote generalized coordinates. For a 2-particle system in 3D:

$$
\begin{pmatrix} \dot{q}_1 \\ \dot{q}_2 \\ \dot{q}_3 \\ \dot{q}_4 \\ \dot{q}_5 \\ \dot{q}_6 \end{pmatrix} = \begin{pmatrix} \dot{x}_1 \\ \dot{y}_1 \\ \dot{z}_1 \\ \dot{x}_2 \\ \dot{y}_2 \\ \dot{z}_2 \end{pmatrix}
$$

The generalized coordinates should not be confused with the electric charges. The Hamiltonian is defined as:

$$
H = \left(\sum_{i=1}^{6} \dot{q}_{i}\frac{\partial L}{\partial \dot{q}_{i}}\right) - L
$$

Performing the algebraic expansion yields:

$$
H = T + U
$$

This is the conserved energy expressed in **positions and physical
velocities**. It becomes a canonical Hamiltonian only after the physical
velocities have been eliminated in favour of the canonical momenta using the
relations above. In particular, it must not be turned into a function of
$(\vec{r}, \vec{p})$ by substituting $\vec{v}_i = \vec{p}_i/m_i$.

## Hamiltonian formulation

For an n-body system of particles interacting via Weber electrodynamics, the kinetic energy for each particle is

$$
T_i = \frac{1}{2} m_i (\vec{v}_i \cdot \vec{v}_i)
$$

The potential energy for a pair of particles is given by (in SI units):

$$
U_{ij} = \frac{1}{4\pi\epsilon_0} \frac{q_i q_j}{r_{ij}}\left(1 - \frac{\dot{r}_{ij}^2}{2c^2}\right)
$$

where $r_{ij}$ is the distance between particles $i$ and $j$, and $\dot{r}_{ij}$ is its time derivative. Here $q_i$ and $q_j$ denote the electric charges of particles $i$ and $j$.

The conserved energy for n charges, in positions and physical velocities, is

$$
E = \sum_{i=1}^{n} T_i + \sum_{i\lt j}^{n} U_{ij}
$$

where $E$ is the total energy, $T_i$ is the kinetic energy of particle $i$, and
$U_{ij}$ is the Weber pair potential.

### Canonical momenta and the canonical Hamiltonian

Differentiating the n-particle Weber Lagrangian with respect to a particle
velocity gives the canonical momentum

$$
\vec{p}_i = m_i\vec{v}_i
- \sum_{j\ne i}\frac{q_i q_j}{c^2}\frac{\dot{r}_{ij}}{r_{ij}^2}
(\vec{r}_i - \vec{r}_j)
$$

where $\vec{p}_i$ is the canonical momentum of particle $i$, $m_i$ its mass,
$\vec{v}_i$ its physical velocity, and each term of the sum is the Weber
correction contributed by a pair containing particle $i$.

These equations are linear in all particle velocities, but different pairs
share particle velocities, so the velocities must be recovered together. It is
convenient to write $k_{ij} = q_i q_j/(c^2 r_{ij})$ and introduce, for each
pair, the *naive* radial rate

$$
s_{ij} = \hat{r}_{ij}\mathbin{\cdot}
\left(\frac{\vec{p}_i}{m_i} - \frac{\vec{p}_j}{m_j}\right)
$$

which is what one would wrongly call $\dot{r}_{ij}$ if canonical momentum were
kinetic momentum. The physical radial velocities $\dot{r}_{ij}$ then solve the
$n(n-1)/2$ dimensional linear system

$$
\dot{r}_{ab} = s_{ab}
+ \sum_{\text{pairs } cd} G_{ab,cd}\,k_{cd}\,\dot{r}_{cd},
\qquad
G_{ab,cd} = (\hat{r}_{ab}\mathbin{\cdot}\hat{r}_{cd})
\left(
\frac{\delta_{ac} - \delta_{ad}}{m_a}
- \frac{\delta_{bc} - \delta_{bd}}{m_b}
\right)
$$

where $G$ couples pairs through the particles they share and $\delta$ is the
Kronecker delta. Once $\dot{r}_{ij}$ is known the physical velocities follow
explicitly,

$$
\vec{v}_i = \frac{\vec{p}_i}{m_i}
+ \frac{1}{m_i}\sum_{j\ne i} k_{ij}\,\dot{r}_{ij}\,\hat{r}_{ij}.
$$

The exact canonical Hamiltonian, obtained by the Legendre transform, is

$$
H(\vec{r},\vec{p}) = \frac{1}{2}\sum_i \vec{p}_i\mathbin{\cdot}\vec{v}_i
+ \sum_{i\lt j}\frac{q_i q_j}{r_{ij}}
= \sum_i \frac{\lVert\vec{p}_i\rVert^2}{2m_i}
+ \frac{1}{2}\sum_{i\lt j} k_{ij}\,\dot{r}_{ij}\,s_{ij}
+ \sum_{i\lt j}\frac{q_i q_j}{r_{ij}}
$$

where the second form is the one evaluated by the package. Both equal the
velocity-space energy $E$ above once $\vec{v}(\vec{r},\vec{p})$ is substituted.
In the Coulomb limit $c\to\infty$, and at any instant where every pair radial
velocity vanishes, $k_{ij}\to 0$ and $H$ reduces to
$\sum_i\lVert\vec{p}_i\rVert^2/(2m_i) + \sum_{i<j} q_iq_j/r_{ij}$.

## Euler-Lagrange equations

The Euler-Lagrange equation for Weber electrodynamics yields the force on particle 1 in the $x$-direction:

$$
F_{12}^x = \frac{d}{dt}\frac{\partial S}{\partial \dot{x}_1} - \frac{\partial S}{\partial x_1} = q_1 q_2 \frac{(x_1 - x_2)}{r^3} \left(1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right)
$$

with analogous equations for the $y$ and $z$ directions and for particle 2.

## Equations of motion

Hamilton's equations give the time evolution of positions and momenta:

$$
\dot{x}_i = \frac{\partial H}{\partial p_{x_i}}, \quad \dot{p}_{x_i} = -\frac{\partial H}{\partial x_i}
$$

For the two-particle Weber system, these expand to:

$$
\dot{x}_1 = \frac{1}{m_1} \left(p_{x_1} + \frac{q_1 q_2}{c^2} \frac{\dot{r} (x_1 - x_2)}{r^2}\right)
$$

$$
\dot{x}_2 = \frac{1}{m_2} \left(p_{x_2} - \frac{q_1 q_2}{c^2} \frac{\dot{r} (x_1 - x_2)}{r^2}\right)
$$

$$
\dot{p}_{x_1} = \frac{q_1 q_2}{r^2}\left[\frac{(x_1 - x_2)}{r}\left(1 + \frac{3\dot{r}^2}{2 c^2}\right) - \frac{\dot{r}(\dot{x}_1 - \dot{x}_2)}{c^2}\right]
$$

$$
\dot{p}_{x_2} = -\dot{p}_{x_1}
$$

with analogous equations for the $y$ and $z$ directions. In all four equations
$\dot{r}$ and $\dot{x}_1 - \dot{x}_2$ are the **physical** radial and relative
velocities obtained from the canonical momenta by the scalar inverse given
earlier, never $p/m$ surrogates.

For $n$ particles the momentum rate is

$$
\dot{\vec{p}}_i = \sum_{j\ne i}\frac{q_i q_j}{r_{ij}^2}
\left[
\hat{r}_{ij}\left(1 + \frac{3\dot{r}_{ij}^2}{2c^2}\right)
- \frac{\dot{r}_{ij}}{c^2}(\vec{v}_i - \vec{v}_j)
\right]
$$

where every $\dot{r}_{ij}$ and $\vec{v}_i - \vec{v}_j$ comes from the
simultaneous velocity solve.

Canonical momentum rate is not itself the mechanical Weber force. Because the
canonical momentum contains the velocity-dependent pair correction
$\vec{\alpha}$, the mechanical acceleration satisfies
$m_i\dot{\vec{v}}_i = \dot{\vec{p}}_i + \mathrm{d}\vec{\alpha}_i/\mathrm{d}t$;
combining that identity with the equation above reproduces the Weber force law
stated earlier in this document.
