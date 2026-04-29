#!/usr/bin/env julia
# api_showcase.jl
#
# Plain-script mirror of `examples/api_showcase.ipynb` — kept in sync so the
# notebook can be regenerated and the same code can be run headlessly. The
# notebook lives alongside `examples/two_body_reference.ipynb` and exercises
# the v0.5.0 API surfaces that the reference notebook does not touch:
#
#   1. Custom Hamiltonian via `kinetic_term + coulomb_term`
#   2. Term introspection (`term_names`, `has_term`, `pair_decomposition`)
#   3. Accessor tour
#   4. N=3 mixed signs with a neutral particle (Zöllner)
#   5. 3D regularization with `:adaptive_cartesian`
#   6. CollisionBounce callback
#   7. Combined RegularizedIntegrator + CollisionBounce
#
# Run from the project root:
#   julia --project=. examples/api_showcase.jl

using WeberElectrodynamics
using LinearAlgebra: norm
using Symbolics: @variables, Num
using Printf: @printf

println("=" ^ 72)
println("WeberElectrodynamics.jl v0.5.0 API showcase")
println("=" ^ 72)

# ---------------------------------------------------------------------------
# Section 1 — Custom Hamiltonian via term builders
# ---------------------------------------------------------------------------
println("\n[1] Custom Hamiltonian — kinetic_term + coulomb_term (pure Coulomb)")

@variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 qq1 qq2 cc tt
q_syms = [x1, y1, x2, y2]
p_syms = [px1, py1, px2, py2]

H_pure_coulomb =
    kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2) +
    coulomb_term(q_syms; charges = [qq1, qq2], n_particles = 2, dims = 2)

sys_coulomb = HamiltonianSystem(
    H_pure_coulomb,
    q_syms,
    p_syms;
    param_symbols = [m1, m2, qq1, qq2, cc],
    kappa_symbols = Num[],
    t = tt,
    n_particles = 2,
    dims = 2,
)

@printf "  H = %s\n" string(sys_coulomb.hamiltonian_symbolic)

# Bound elliptic orbit with Kepler ICs.
m1_v, m2_v = 1.0, 1.0
q1_v, q2_v = 1.0, -1.0
r0 = 2.0
M = m1_v + m2_v
μ_red = m1_v * m2_v / M
v_circ = sqrt(abs(q1_v * q2_v) / (μ_red * r0))
v_init = 0.95 * v_circ
q0 = [-m2_v / M * r0, 0.0, m1_v / M * r0, 0.0]
p0 = [0.0, m1_v * (-m2_v / M * v_init), 0.0, m2_v * (m1_v / M * v_init)]
E_keplerian = 0.5 * (sum(p0[1:2] .^ 2) / m1_v + sum(p0[3:4] .^ 2) / m2_v) + q1_v * q2_v / r0
T_kepler = 2π * sqrt(μ_red * (-q1_v * q2_v / (-2 * E_keplerian))^3 / (-q1_v * q2_v))

prob_coulomb = HamiltonianProblem(
    sys_coulomb,
    (0.0, T_kepler),
    q0,
    p0;
    masses = [m1_v, m2_v],
    charges = [q1_v, q2_v],
    c = 1.0,                     # unused by pure Coulomb but required by the params layout
    dt = T_kepler / 10_000,
)
sol_coulomb = solve(prob_coulomb)
@printf "  retcode=%s  steps=%d  Δq(end-start)=%.2e\n" sol_coulomb.retcode length(
    sol_coulomb.t,
) maximum(abs, sol_coulomb.q[end] .- q0)

# ---------------------------------------------------------------------------
# Section 2 — Term introspection
# ---------------------------------------------------------------------------
println("\n[2] Term introspection — default Weber+Zöllner system")

sys_default = HamiltonianSystem(2, 2)
@printf "  term_names(sys) = %s\n" term_names(sys_default)
@printf "  has_term(sys, :weber) = %s\n" has_term(sys_default, :weber)
@printf "  has_term(sys, :zollner) = %s\n" has_term(sys_default, :zollner)

prob_default = HamiltonianProblem(
    sys_default,
    (0.0, 0.5),
    [1.0, 0.0, -1.0, 0.0],
    [0.0, 0.2, 0.0, -0.2];
    masses = [1.0, 1.0],
    charges = [1.0, -1.0],
    c = 5.0,
    dt = 0.005,
    zollner = ZollnerOptions(enabled = true, a = 0.05),
)
weber = get_term(sys_default, :weber)
decomp = weber.pair_decomposition(
    1,
    2,
    prob_default.q_initial,
    prob_default.p_initial,
    params(prob_default),
    kappas(prob_default),
)
@printf "  Weber pair decomposition (1,2): coulomb=%.4e velocity=%.4e r=%.3f rdot=%.3e\n" decomp.coulomb decomp.velocity decomp.r decomp.rdot

# ---------------------------------------------------------------------------
# Section 3 — Accessor tour
# ---------------------------------------------------------------------------
println("\n[3] Accessor tour")
@printf "  n_particles(prob) = %d\n" n_particles(prob_default)
@printf "  dims(prob) = %d\n" dims(prob_default)
@printf "  masses(prob) = %s\n" string(collect(masses(prob_default)))
@printf "  charges(prob) = %s\n" string(collect(charges(prob_default)))
@printf "  speed_of_light(prob) = %g\n" speed_of_light(prob_default)
@printf "  kappas(prob) = %s\n" string(kappas(prob_default))
@printf "  kappa(prob, 1, 2) = %g\n" kappa(prob_default, 1, 2)
@printf "  parent(masses(prob)) === params(prob): %s\n" (
    parent(masses(prob_default)) === params(prob_default)
)

# ---------------------------------------------------------------------------
# Section 4 — N=3 mixed signs with a neutral particle
# ---------------------------------------------------------------------------
println("\n[4] N=3 with charges [+1, 0, -1] and Zöllner enabled")
sys_n3 = HamiltonianSystem(3, 2)
prob_n3 = HamiltonianProblem(
    sys_n3,
    (0.0, 0.4),
    [1.0, 0.0, 0.0, 0.0, -1.0, 0.0],
    zeros(6);
    masses = [1.0, 1.0, 1.0],
    charges = [1.0, 0.0, -1.0],
    c = 5.0,
    dt = 0.002,
    zollner = ZollnerOptions(enabled = true, a = 0.07),
)
@printf "  kappas(prob) = %s   (neutral pairs keep κ=1; only opposite nonzero signs get κ=1+a)\n" string(
    kappas(prob_n3),
)
sol_n3 = solve(prob_n3)
@printf "  retcode=%s  steps=%d\n" sol_n3.retcode length(sol_n3.t)

# ---------------------------------------------------------------------------
# Section 5 — 3D regularization with :adaptive_cartesian
# ---------------------------------------------------------------------------
println("\n[5] 3D close approach with :adaptive_cartesian regularization")
sys_3d = HamiltonianSystem(2, 3)
m1_3d, m2_3d = 1.0, 0.1
q1_3d, q2_3d = sqrt(0.1), -sqrt(0.1)
c_3d = 4.0
r0_3d = 2.0
M_3d = m1_3d + m2_3d
v_circ_3d = sqrt(abs(q1_3d * q2_3d) * M_3d / (m1_3d * m2_3d * r0_3d))
v_3d = 0.5 * v_circ_3d
q0_3d = [-m2_3d / M_3d * r0_3d, 0.0, 0.0, m1_3d / M_3d * r0_3d, 0.0, 0.0]
p0_3d = [0.0, m1_3d * (-m2_3d / M_3d * v_3d), 0.0, 0.0, m2_3d * (m1_3d / M_3d * v_3d), 0.0]
prob_3d = HamiltonianProblem(
    sys_3d,
    (0.0, 8.0),
    q0_3d,
    p0_3d;
    masses = [m1_3d, m2_3d],
    charges = [q1_3d, q2_3d],
    c = c_3d,
    dt = 2e-3,
)
alg_3d = RegularizedIntegrator(
    SymmetricProjectionIntegrator();
    backend = :adaptive_cartesian,
    r_on_factor = 0.3,
    r_off_factor = 0.45,
    warn_on_fallback = false,
)
sol_3d = solve(prob_3d, alg_3d)
diag = sol_3d.regularization
@printf "  retcode=%s  steps=%d\n" sol_3d.retcode length(sol_3d.t)
@printf "  regularization: backend=%s activations=%d active_steps=%d total_substeps=%d min_r=%.3e\n" diag.used_backend diag.activation_count diag.active_steps diag.total_substeps diag.min_encounter_distance

# ---------------------------------------------------------------------------
# Section 6 — CollisionBounce callback (sub-critical like-charge oscillation)
# ---------------------------------------------------------------------------
println("\n[6] Sub-critical like-charge bounce (unregularized + CollisionBounce)")
# Frauenfelder–Weber sub-critical regime: like-sign charges fall in through
# near-zero separation (Weber's velocity term overcomes Coulomb repulsion).
# CollisionBounce reflects the relative coordinate when r < bounce_r,
# C⁰-continuing the head-on (ℓ=0) collision. Bounded energy error requires
# a dt small enough that the bounce fires before the projection diverges
# (validated: dt=1e-4 with bounce_r=0.02 → ~0.01% drift over 100 periods).
m_b = 1.0
q_b = 1.0
c_b = 4.0
r0_b = 0.05
mu_b = m_b * m_b / (m_b + m_b)
rho_b = q_b * q_b / (mu_b * c_b^2)
T_b = (2π / c_b) * rho_b
sys_b = HamiltonianSystem(2, 2)
prob_b = HamiltonianProblem(
    sys_b,
    (0.0, 3 * T_b),
    [r0_b / 2, 0.0, -r0_b / 2, 0.0],
    [0.0, 0.0, 0.0, 0.0];
    masses = [m_b, m_b],
    charges = [q_b, q_b],
    c = c_b,
    dt = 1e-4,
)
sol_b = solve(prob_b, SymmetricProjectionIntegrator(); callbacks = CollisionBounce(0.02))
en_b = compute_energy_timeseries(sol_b)
ΔE_b = abs(en_b.statistics.global_error_percent_max)
@printf "  retcode=%s  steps=%d  max global energy drift = %.3e %%\n" sol_b.retcode length(
    sol_b.t,
) ΔE_b

# ---------------------------------------------------------------------------
# Section 7 — RegularizedIntegrator synthesises a CollisionBounce automatically
# ---------------------------------------------------------------------------
println("\n[7] collision_bounce_radius kwarg synthesises an equivalent callback")
alg_synth = RegularizedIntegrator(
    SymmetricProjectionIntegrator();
    backend = :adaptive_cartesian,
    r_on = 0.001,
    r_off = 0.002,
    warn_on_fallback = false,
    collision_bounce_radius = 0.02,
)
integrator = init(prob_b, alg_synth)
@printf "  integrator.callbacks length = %d (auto-synthesised)\n" length(
    integrator.callbacks,
)
@printf "  integrator.callbacks[1] isa CollisionBounce = %s\n" (
    integrator.callbacks[1] isa CollisionBounce
)
@printf "  radius = %g\n" integrator.callbacks[1].radius

# Synthesised vs explicit CollisionBounce on a regularized integrator should
# match bit-for-bit (proven by `equivalence: collision_bounce_radius
# synthesises CollisionBounce` testset in test/test_callbacks.jl).
alg_explicit = RegularizedIntegrator(
    SymmetricProjectionIntegrator();
    backend = :adaptive_cartesian,
    r_on = 0.001,
    r_off = 0.002,
    warn_on_fallback = false,
)
sol_synth = solve(prob_b, alg_synth)
sol_explicit = solve(prob_b, alg_explicit; callbacks = CollisionBounce(0.02))
max_dq = maximum(
    maximum(abs, sol_synth.q[i] .- sol_explicit.q[i]) for i in eachindex(sol_synth.q)
)
max_dp = maximum(
    maximum(abs, sol_synth.p[i] .- sol_explicit.p[i]) for i in eachindex(sol_synth.p)
)
@printf "  bit-exact match against explicit CollisionBounce: max|Δq|=%.2e max|Δp|=%.2e\n" max_dq max_dp

println("\n" * "=" ^ 72)
println("All sections completed cleanly.")
