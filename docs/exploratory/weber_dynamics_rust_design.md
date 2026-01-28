# Rust Crate Design: `weber-dynamics`

A Rust crate implementing symbolic Weber electrodynamics with semi-explicit symplectic integration, powered by Symbolica.

---

## 1. Architecture Overview

```
weber-dynamics/
├── Cargo.toml
├── src/
│   ├── lib.rs                    # Public API re-exports
│   ├── symbolic/
│   │   ├── mod.rs
│   │   ├── hamiltonian.rs        # Weber Hamiltonian construction
│   │   ├── equations.rs          # Equations of motion via differentiation
│   │   └── latex.rs              # LaTeX formatting
│   ├── numeric/
│   │   ├── mod.rs
│   │   ├── state.rs              # Phase space state representation
│   │   ├── evaluator.rs          # Compiled expression evaluators
│   │   └── integrator.rs         # Semi-explicit symplectic integrator
│   └── particles/
│       ├── mod.rs
│       └── system.rs             # N-body particle system definition
└── examples/
    ├── two_body.rs
    ├── latex_output.rs
    └── jupyter_integration.rs
```

---

## 2. Core Types

### 2.1 Particle System Definition

```rust
/// Physical parameters for a charged particle
#[derive(Clone, Debug)]
pub struct Particle {
    pub mass: f64,      // m_i in absolute units (mg)
    pub charge: f64,    // q_i in absolute units
    pub label: String,  // e.g., "1", "2", "electron"
}

/// N-body system configuration
pub struct ParticleSystem {
    particles: Vec<Particle>,
    speed_of_light: f64,  // c in absolute units (mm/s)
}

impl ParticleSystem {
    pub fn new(particles: Vec<Particle>, c: f64) -> Self;
    pub fn n_particles(&self) -> usize;
    pub fn n_dof(&self) -> usize;  // 3 * n_particles (3D)
}
```

### 2.2 Phase Space State

```rust
/// State vector in phase space: (q, p) ∈ ℝ^{2d} where d = 3n
pub struct PhaseState {
    pub positions: Vec<[f64; 3]>,   // [r_1, r_2, ..., r_n]
    pub momenta: Vec<[f64; 3]>,     // [p_1, p_2, ..., p_n]
}

/// Extended phase space state: (q, x, p, y) ∈ ℝ^{4d}
pub struct ExtendedState {
    pub q: Vec<[f64; 3]>,
    pub x: Vec<[f64; 3]>,
    pub p: Vec<[f64; 3]>,
    pub y: Vec<[f64; 3]>,
}
```

---

## 3. Symbolic Module

### 3.1 Weber Hamiltonian Construction

The Weber Hamiltonian for $n$ particles:

$$H = \sum_{i=1}^{n} T_i + \sum_{i<j}^{n} U_{ij}$$

where kinetic energy:

$$T_i = \frac{|\vec{p}_i|^2}{2m_i}$$

and Weber potential energy:

$$U_{ij} = \frac{q_i q_j}{r_{ij}} \left(1 - \frac{\dot{r}_{ij}^2}{2c^2}\right)$$

**Key insight**: The radial velocity $\dot{r}_{ij}$ depends on velocities, making this a **non-separable Hamiltonian** where $H = H(q, p)$ cannot be split into $T(p) + V(q)$.

```rust
use symbolica::{atom::Atom, symbol, parse};

pub struct WeberHamiltonian {
    /// The full symbolic Hamiltonian H(q, p)
    pub expr: Atom,

    /// Position symbols: x_i, y_i, z_i for each particle
    pub position_symbols: Vec<[Atom; 3]>,

    /// Momentum symbols: px_i, py_i, pz_i for each particle
    pub momentum_symbols: Vec<[Atom; 3]>,

    /// System parameters
    pub system: ParticleSystem,
}

impl WeberHamiltonian {
    /// Build the symbolic Hamiltonian for n particles
    pub fn build(system: &ParticleSystem) -> Self {
        let n = system.n_particles();

        // Create position symbols: x_1, y_1, z_1, x_2, ...
        let positions = (0..n).map(|i| {
            let idx = i + 1;
            [
                symbol!(format!("x_{}", idx)),
                symbol!(format!("y_{}", idx)),
                symbol!(format!("z_{}", idx)),
            ]
        }).collect::<Vec<_>>();

        // Create momentum symbols: px_1, py_1, pz_1, ...
        let momenta = (0..n).map(|i| {
            let idx = i + 1;
            [
                symbol!(format!("px_{}", idx)),
                symbol!(format!("py_{}", idx)),
                symbol!(format!("pz_{}", idx)),
            ]
        }).collect::<Vec<_>>();

        // Build kinetic terms: T_i = (px_i² + py_i² + pz_i²) / (2 * m_i)
        let kinetic = build_kinetic_terms(&momenta, &system.particles);

        // Build potential terms: U_ij for all pairs
        let potential = build_weber_potential(&positions, &momenta, system);

        WeberHamiltonian {
            expr: kinetic + potential,
            position_symbols: positions,
            momentum_symbols: momenta,
            system: system.clone(),
        }
    }
}
```

### 3.2 Helper: Building Weber Potential Symbolically

For each pair $(i, j)$, we need:

$$r_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$$

$$\dot{r}_{ij} = \frac{\vec{r}_{ij} \cdot \vec{v}_{ij}}{r_{ij}}$$

where velocities come from momenta: $\vec{v}_i = \vec{p}_i / m_i$.

```rust
fn build_weber_potential(
    positions: &Vec<[Atom; 3]>,
    momenta: &Vec<[Atom; 3]>,
    system: &ParticleSystem,
) -> Atom {
    let mut potential = parse!("0");
    let c = system.speed_of_light;
    let c2 = c * c;

    for i in 0..system.n_particles() {
        for j in (i+1)..system.n_particles() {
            let q_i = system.particles[i].charge;
            let q_j = system.particles[j].charge;
            let m_i = system.particles[i].mass;
            let m_j = system.particles[j].mass;

            // Relative position: r_ij = r_i - r_j
            let dx = &positions[i][0] - &positions[j][0];
            let dy = &positions[i][1] - &positions[j][1];
            let dz = &positions[i][2] - &positions[j][2];

            // r_ij = sqrt(dx² + dy² + dz²)
            let r_ij_sq = &dx * &dx + &dy * &dy + &dz * &dz;
            let r_ij = r_ij_sq.clone().sqrt();

            // Relative velocity from momenta
            let dvx = &momenta[i][0] / m_i - &momenta[j][0] / m_j;
            let dvy = &momenta[i][1] / m_i - &momenta[j][1] / m_j;
            let dvz = &momenta[i][2] / m_i - &momenta[j][2] / m_j;

            // rdot_ij = (r · v) / r
            let r_dot_v = &dx * &dvx + &dy * &dvy + &dz * &dvz;
            let rdot_ij = &r_dot_v / &r_ij;
            let rdot_sq = &rdot_ij * &rdot_ij;

            // U_ij = (q_i * q_j / r_ij) * (1 - rdot²/(2c²))
            let coulomb = q_i * q_j / &r_ij;
            let weber_factor = parse!("1") - &rdot_sq / (2.0 * c2);

            potential = potential + coulomb * weber_factor;
        }
    }

    potential
}
```

### 3.3 Equations of Motion via Differentiation

Hamilton's equations:

$$\dot{q}_i = \frac{\partial H}{\partial p_i}, \quad \dot{p}_i = -\frac{\partial H}{\partial q_i}$$

```rust
pub struct EquationsOfMotion {
    /// ∂H/∂p for each momentum component
    pub q_dot: Vec<Atom>,  // length = 3n

    /// -∂H/∂q for each position component
    pub p_dot: Vec<Atom>,  // length = 3n
}

impl WeberHamiltonian {
    /// Derive equations of motion by symbolic differentiation
    pub fn equations_of_motion(&self) -> EquationsOfMotion {
        let n = self.system.n_particles();

        // ∂H/∂p_i for each momentum component
        let q_dot: Vec<Atom> = self.momentum_symbols
            .iter()
            .flat_map(|[px, py, pz]| {
                [
                    self.expr.derivative(px),
                    self.expr.derivative(py),
                    self.expr.derivative(pz),
                ]
            })
            .collect();

        // -∂H/∂q_i for each position component
        let p_dot: Vec<Atom> = self.position_symbols
            .iter()
            .flat_map(|[x, y, z]| {
                [
                    -self.expr.derivative(x),
                    -self.expr.derivative(y),
                    -self.expr.derivative(z),
                ]
            })
            .collect();

        EquationsOfMotion { q_dot, p_dot }
    }
}
```

### 3.4 LaTeX Output

```rust
use symbolica::printer::{PrintOptions, PrintMode};

pub trait ToLatex {
    fn to_latex(&self) -> String;
    fn to_latex_with_options(&self, opts: &LatexOptions) -> String;
}

pub struct LatexOptions {
    /// Use \frac{}{} instead of /
    pub use_frac: bool,
    /// Simplify before printing
    pub simplify: bool,
    /// Maximum terms before collapsing with "..."
    pub max_terms: Option<usize>,
}

impl ToLatex for WeberHamiltonian {
    fn to_latex(&self) -> String {
        // Use Symbolica's printer with LaTeX mode
        let opts = PrintOptions {
            mode: PrintMode::Latex,
            ..Default::default()
        };
        self.expr.to_string_with_options(&opts)
    }
}

impl ToLatex for EquationsOfMotion {
    fn to_latex(&self) -> String {
        let mut output = String::new();

        // Format each equation: \dot{x}_i = ...
        for (i, expr) in self.q_dot.iter().enumerate() {
            let coord = match i % 3 {
                0 => "x",
                1 => "y",
                _ => "z",
            };
            let particle = i / 3 + 1;
            output.push_str(&format!(
                "\\dot{{{}}}_{} &= {} \\\\\n",
                coord, particle, expr.to_latex()
            ));
        }

        output
    }
}
```

---

## 4. Numeric Module

### 4.1 Compiled Evaluators

Symbolica can compile expressions to efficient native code:

```rust
use symbolica::evaluate::{ExpressionEvaluator, CompiledRealEvaluator};

/// Compiled evaluator for the Hamiltonian and its derivatives
pub struct CompiledWeberSystem {
    /// H(q, p) evaluator
    hamiltonian: CompiledRealEvaluator,

    /// ∂H/∂p evaluators (for q_dot)
    dh_dp: Vec<CompiledRealEvaluator>,

    /// ∂H/∂q evaluators (for p_dot, note: returns -∂H/∂q)
    neg_dh_dq: Vec<CompiledRealEvaluator>,

    /// Parameter mapping: symbol name -> index in input array
    var_indices: HashMap<String, usize>,

    /// Dimension d = 3n
    dim: usize,
}

impl CompiledWeberSystem {
    pub fn from_hamiltonian(h: &WeberHamiltonian) -> Self {
        let eom = h.equations_of_motion();

        // Build variable map: [x_1, y_1, z_1, ..., px_1, py_1, pz_1, ...]
        let vars: Vec<Atom> = h.position_symbols
            .iter()
            .flat_map(|s| s.iter().cloned())
            .chain(h.momentum_symbols.iter().flat_map(|s| s.iter().cloned()))
            .collect();

        // Compile H
        let hamiltonian = h.expr
            .evaluator(&vars)
            .compile();

        // Compile ∂H/∂p (q_dot)
        let dh_dp = eom.q_dot
            .iter()
            .map(|expr| expr.evaluator(&vars).compile())
            .collect();

        // Compile -∂H/∂q (p_dot)
        let neg_dh_dq = eom.p_dot
            .iter()
            .map(|expr| expr.evaluator(&vars).compile())
            .collect();

        CompiledWeberSystem {
            hamiltonian,
            dh_dp,
            neg_dh_dq,
            var_indices: build_var_indices(&vars),
            dim: vars.len() / 2,
        }
    }

    /// Evaluate H at given state
    pub fn energy(&self, state: &PhaseState) -> f64 {
        let input = self.state_to_input(state);
        self.hamiltonian.evaluate(&input)
    }

    /// Evaluate q_dot = ∂H/∂p at given state
    pub fn q_dot(&self, state: &PhaseState, out: &mut [f64]) {
        let input = self.state_to_input(state);
        for (i, eval) in self.dh_dp.iter().enumerate() {
            out[i] = eval.evaluate(&input);
        }
    }

    /// Evaluate p_dot = -∂H/∂q at given state
    pub fn p_dot(&self, state: &PhaseState, out: &mut [f64]) {
        let input = self.state_to_input(state);
        for (i, eval) in self.neg_dh_dq.iter().enumerate() {
            out[i] = eval.evaluate(&input);
        }
    }
}
```

### 4.2 Semi-Explicit Symplectic Integrator

Implementation following the algorithm from `SemiExplicitIntegrator.md`:

```rust
/// Configuration for the integrator
pub struct IntegratorConfig {
    /// Integration order (2, 4, 6)
    pub order: usize,
    /// Composition method for higher orders
    pub composition: Composition,
    /// Projection tolerance ε
    pub projection_tol: f64,
    /// Maximum projection iterations
    pub max_projection_iter: usize,
    /// Use warm start for μ from previous step
    pub warm_start: bool,
}

pub enum Composition {
    /// 3-stage triple jump: 3^(n/2) stages
    TripleJump,
    /// 5-stage Suzuki composition
    Suzuki,
    /// 7-stage Yoshida (6th order only)
    Yoshida6,
}

/// Semi-explicit symplectic integrator for non-separable Hamiltonians
pub struct SemiExplicitIntegrator {
    /// Compiled system for fast evaluation
    system: CompiledWeberSystem,
    /// Configuration
    config: IntegratorConfig,
    /// Cached μ for warm start
    cached_mu: Vec<f64>,
    /// Dimension d = 3n
    dim: usize,
}
```

#### 4.2.1 Extended Phase Space Operations

```rust
impl SemiExplicitIntegrator {
    /// Embed z = (q, p) into extended space Z = (q, q, p, p)
    fn embed(&self, state: &PhaseState) -> ExtendedState {
        ExtendedState {
            q: state.positions.clone(),
            x: state.positions.clone(),
            p: state.momenta.clone(),
            y: state.momenta.clone(),
        }
    }

    /// Project: apply A·Z to get constraint violation
    /// A = [I -I 0 0; 0 0 I -I]
    /// Returns (q - x, p - y) ∈ ℝ^{2d}
    fn constraint_violation(&self, z: &ExtendedState) -> Vec<f64> {
        let d = self.dim;
        let mut result = vec![0.0; 2 * d];

        // First d components: q - x
        for i in 0..z.q.len() {
            for j in 0..3 {
                result[3*i + j] = z.q[i][j] - z.x[i][j];
            }
        }

        // Second d components: p - y
        for i in 0..z.p.len() {
            for j in 0..3 {
                result[d + 3*i + j] = z.p[i][j] - z.y[i][j];
            }
        }

        result
    }

    /// Apply A^T·μ shift to extended state
    /// A^T = [I 0; -I 0; 0 I; 0 -I]
    fn apply_shift(&self, z: &mut ExtendedState, mu: &[f64]) {
        let d = self.dim;
        let mu_q = &mu[0..d];     // Shift for position components
        let mu_p = &mu[d..2*d];   // Shift for momentum components

        for i in 0..z.q.len() {
            for j in 0..3 {
                let idx = 3*i + j;
                z.q[i][j] += mu_q[idx];
                z.x[i][j] -= mu_q[idx];
                z.p[i][j] += mu_p[idx];
                z.y[i][j] -= mu_p[idx];
            }
        }
    }
}
```

#### 4.2.2 Flow Maps

```rust
impl SemiExplicitIntegrator {
    /// Flow A: evolves (x, p) using H(q, y), freezes (q, y)
    /// x += dt * ∂H/∂p(q, y)
    /// p += dt * (-∂H/∂q)(q, y)  [note: p_dot already includes minus sign]
    fn flow_a(&self, z: &mut ExtendedState, dt: f64) {
        // Create state from (q, y) for evaluation
        let eval_state = PhaseState {
            positions: z.q.clone(),
            momenta: z.y.clone(),
        };

        // Get derivatives at (q, y)
        let mut q_dot = vec![0.0; self.dim];
        let mut p_dot = vec![0.0; self.dim];
        self.system.q_dot(&eval_state, &mut q_dot);
        self.system.p_dot(&eval_state, &mut p_dot);

        // Update x and p
        for i in 0..z.x.len() {
            for j in 0..3 {
                let idx = 3*i + j;
                z.x[i][j] += dt * q_dot[idx];
                z.p[i][j] += dt * p_dot[idx];
            }
        }
    }

    /// Flow B: evolves (q, y) using H(x, p), freezes (x, p)
    fn flow_b(&self, z: &mut ExtendedState, dt: f64) {
        // Create state from (x, p) for evaluation
        let eval_state = PhaseState {
            positions: z.x.clone(),
            momenta: z.p.clone(),
        };

        // Get derivatives at (x, p)
        let mut q_dot = vec![0.0; self.dim];
        let mut p_dot = vec![0.0; self.dim];
        self.system.q_dot(&eval_state, &mut q_dot);
        self.system.p_dot(&eval_state, &mut p_dot);

        // Update q and y
        for i in 0..z.q.len() {
            for j in 0..3 {
                let idx = 3*i + j;
                z.q[i][j] += dt * q_dot[idx];
                z.y[i][j] += dt * p_dot[idx];
            }
        }
    }

    /// 2nd-order Strang splitting: Φ = A_{dt/2} ∘ B_{dt} ∘ A_{dt/2}
    fn strang_step(&self, z: &mut ExtendedState, dt: f64) {
        self.flow_a(z, dt / 2.0);
        self.flow_b(z, dt);
        self.flow_a(z, dt / 2.0);
    }
}
```

#### 4.2.3 Higher-Order Composition

```rust
impl SemiExplicitIntegrator {
    /// Triple jump coefficients for order n
    fn triple_jump_coeffs(n: usize) -> [f64; 3] {
        let w = 2.0_f64.powf(1.0 / (n as f64 - 1.0));
        let gamma1 = 1.0 / (2.0 - w);
        let gamma2 = -w / (2.0 - w);
        [gamma1, gamma2, gamma1]  // γ₁, γ₂, γ₃
    }

    /// Recursive composition: builds order n from order n-2
    fn composed_step(&self, z: &mut ExtendedState, dt: f64, order: usize) {
        if order == 2 {
            self.strang_step(z, dt);
            return;
        }

        let coeffs = Self::triple_jump_coeffs(order);
        self.composed_step(z, coeffs[0] * dt, order - 2);
        self.composed_step(z, coeffs[1] * dt, order - 2);
        self.composed_step(z, coeffs[2] * dt, order - 2);
    }
}
```

#### 4.2.4 Symmetric Projection

The key algorithm step: find $\mu$ such that the evolved state lands back on the invariant subspace $\mathcal{N}$ where $q = x$ and $p = y$.

```rust
impl SemiExplicitIntegrator {
    /// Solve for μ using simplified Newton iteration
    /// f(μ) = A·(Φ(Z + A^T·μ) + A^T·μ) = 0
    fn solve_projection(&self, z0: &ExtendedState, dt: f64) -> Vec<f64> {
        let d2 = 2 * self.dim;
        let mut mu = if self.config.warm_start {
            self.cached_mu.clone()
        } else {
            vec![0.0; d2]
        };

        for iter in 0..self.config.max_projection_iter {
            // Compute f(μ) = A·(Φ(Z + A^T·μ) + A^T·μ)
            let mut z = z0.clone();
            self.apply_shift(&mut z, &mu);
            self.composed_step(&mut z, dt, self.config.order);
            self.apply_shift(&mut z, &mu);

            let f = self.constraint_violation(&z);

            // Check convergence
            let f_norm: f64 = f.iter().map(|x| x * x).sum::<f64>().sqrt();
            if f_norm < self.config.projection_tol {
                break;
            }

            // Simplified Newton: μ_{k+1} = μ_k - (1/4)·f(μ_k)
            // Since Df₀ = 4·I_{2d}
            for i in 0..d2 {
                mu[i] -= 0.25 * f[i];
            }
        }

        mu
    }

    /// Full integration step: z_n → z_{n+1}
    pub fn step(&mut self, state: &PhaseState, dt: f64) -> PhaseState {
        // 1. Embed to extended space
        let z0 = self.embed(state);

        // 2. Solve for projection shift μ
        let mu = self.solve_projection(&z0, dt);

        // 3. Shift starting point
        let mut z = z0;
        self.apply_shift(&mut z, &mu);

        // 4. Evolve with composed integrator
        self.composed_step(&mut z, dt, self.config.order);

        // 5. Shift back to subspace
        self.apply_shift(&mut z, &mu);

        // 6. Cache μ for warm start
        if self.config.warm_start {
            self.cached_mu = mu;
        }

        // 7. Extract result (q, p) from (q, q, p, p)
        PhaseState {
            positions: z.q,
            momenta: z.p,
        }
    }
}
```

---

## 5. High-Level API

### 5.1 Builder Pattern

```rust
pub struct WeberSimulation {
    system: ParticleSystem,
    hamiltonian: WeberHamiltonian,
    compiled: CompiledWeberSystem,
    integrator: SemiExplicitIntegrator,
}

impl WeberSimulation {
    pub fn builder() -> WeberSimulationBuilder {
        WeberSimulationBuilder::default()
    }

    /// Get symbolic Hamiltonian for inspection
    pub fn hamiltonian(&self) -> &WeberHamiltonian {
        &self.hamiltonian
    }

    /// Compute energy at current state
    pub fn energy(&self, state: &PhaseState) -> f64 {
        self.compiled.energy(state)
    }

    /// Integrate for given time with given step size
    pub fn integrate(
        &mut self,
        initial: PhaseState,
        total_time: f64,
        dt: f64,
    ) -> Vec<(f64, PhaseState)> {
        let mut trajectory = Vec::new();
        let mut state = initial;
        let mut t = 0.0;

        trajectory.push((t, state.clone()));

        while t < total_time {
            state = self.integrator.step(&state, dt);
            t += dt;
            trajectory.push((t, state.clone()));
        }

        trajectory
    }
}

pub struct WeberSimulationBuilder {
    particles: Vec<Particle>,
    speed_of_light: f64,
    integrator_order: usize,
    composition: Composition,
    projection_tol: f64,
}

impl WeberSimulationBuilder {
    pub fn add_particle(mut self, mass: f64, charge: f64) -> Self {
        self.particles.push(Particle {
            mass,
            charge,
            label: format!("{}", self.particles.len() + 1),
        });
        self
    }

    pub fn speed_of_light(mut self, c: f64) -> Self {
        self.speed_of_light = c;
        self
    }

    pub fn order(mut self, order: usize) -> Self {
        self.integrator_order = order;
        self
    }

    pub fn build(self) -> WeberSimulation {
        let system = ParticleSystem::new(self.particles, self.speed_of_light);
        let hamiltonian = WeberHamiltonian::build(&system);
        let compiled = CompiledWeberSystem::from_hamiltonian(&hamiltonian);
        let integrator = SemiExplicitIntegrator::new(
            compiled.clone(),
            IntegratorConfig {
                order: self.integrator_order,
                composition: self.composition,
                projection_tol: self.projection_tol,
                max_projection_iter: 10,
                warm_start: true,
            },
        );

        WeberSimulation {
            system,
            hamiltonian,
            compiled,
            integrator,
        }
    }
}
```

### 5.2 Example Usage

```rust
use weber_dynamics::prelude::*;

fn main() {
    // Build a two-body Weber system
    let mut sim = WeberSimulation::builder()
        .add_particle(1.0, 1.0)   // m=1, q=1 (proton-like)
        .add_particle(1.0, -1.0)  // m=1, q=-1 (electron-like)
        .speed_of_light(299792.458)  // mm/s
        .order(4)
        .build();

    // Print Hamiltonian as LaTeX
    println!("Hamiltonian:");
    println!("{}", sim.hamiltonian().to_latex());

    // Print equations of motion
    let eom = sim.hamiltonian().equations_of_motion();
    println!("\nEquations of motion:");
    println!("{}", eom.to_latex());

    // Set initial conditions
    let initial = PhaseState {
        positions: vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        momenta: vec![[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]],
    };

    // Integrate
    let trajectory = sim.integrate(initial, 100.0, 0.01);

    // Check energy conservation
    let e0 = sim.energy(&trajectory[0].1);
    let ef = sim.energy(&trajectory.last().unwrap().1);
    println!("\nEnergy conservation: ΔE/E = {:.2e}", (ef - e0).abs() / e0.abs());
}
```

---

## 6. Jupyter Integration

For interactive exploration in Jupyter notebooks via evcxr:

```rust
// In evcxr Jupyter kernel
:dep weber-dynamics = "0.1"
:dep evcxr_jupyter

use weber_dynamics::prelude::*;
use evcxr_jupyter::display;

let sim = WeberSimulation::builder()
    .add_particle(1.0, 1.0)
    .add_particle(1.0, -1.0)
    .speed_of_light(299792.458)
    .build();

// Display Hamiltonian as rendered LaTeX
display::latex(&sim.hamiltonian().to_latex());

// Display equations of motion
display::latex(&format!(
    r"\begin{{align}} {} \end{{align}}",
    sim.hamiltonian().equations_of_motion().to_latex()
));
```

---

## 7. Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| **Symbolica for CAS** | Mature Rust-native CAS with differentiation, pattern matching, and compiled evaluators |
| **Compiled evaluators** | Avoid symbolic overhead during integration; Symbolica compiles to native code |
| **Extended phase space** | Enables explicit flow maps for non-separable Hamiltonians |
| **Symmetric projection** | Ensures symplecticity on original phase space; fast convergence |
| **Warm start for μ** | Reduces projection iterations from 2-3 to 1-2 per step |
| **Builder pattern** | Ergonomic API for configuring particle systems and integrator |
| **Separate symbolic/numeric modules** | Clear separation: build once symbolically, run fast numerically |

---

## 8. Verification Plan

1. **Two-body Coulomb limit**: Set $c \to \infty$, verify Kepler ellipses
2. **Energy conservation**: Check $|H(t) - H(0)| / |H(0)| < \varepsilon$ over long integrations
3. **Symplecticity**: Verify phase space volume preservation
4. **Order verification**: Confirm $O(\Delta t^n)$ global error scaling
5. **Compare with analytical**: Test radial motion with known solutions

---

## 9. Dependencies

```toml
[dependencies]
symbolica = "1.2"
nalgebra = "0.33"        # Linear algebra for projections
num-traits = "0.2"
thiserror = "1.0"

[dev-dependencies]
criterion = "0.5"        # Benchmarking
approx = "0.5"           # Float comparisons in tests
```
