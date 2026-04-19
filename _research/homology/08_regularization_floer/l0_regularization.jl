# l0_regularization.jl
#
# Agent 08: Test whether a radial (ell=0) like-charge trajectory can be
# integrated through the critical radius rho.
#
# Physics: Two like charges, ell=0, reduced mass mu, charges e1=e2=e.
# Critical radius rho = e^2/(mu c^2). The radial Hamiltonian:
#
#   H = (1/(2mu)) * (r/(r-rho)) * p_r^2 + e^2/r
#
# The effective mass m_eff(r) = mu * (1 - rho/r)^{-1} diverges at r=rho
# and is negative for r<rho.
#
# The energy equation (Frauenfelder-Weber 2024):
#   r_dot^2 = (2(e^2/mu) r - 2 h r^2) / (r (rho - r))   [for ell=0, r<rho]
#           = (2(e^2/mu) - 2 h r) / (rho - r)
#
# We test three things:
# 1. Super-critical approach: start above rho, fall inward. Does r reach rho?
# 2. Sub-critical dynamics: start below rho, integrate the ODE.
# 3. Through-critical matching: can we join the two at r=rho?

using Printf

# ============================================================================
# Physical parameters (dimensionless: e^2/mu = 1, c = 1, rho = 1)
# ============================================================================
const RHO = 1.0      # critical radius (= e^2/(mu c^2) in natural units)
const E2_MU = 1.0     # e^2/mu
const C = 1.0

# ============================================================================
# Part 1: Energy equation analysis
# ============================================================================
function radial_velocity_squared(r, h)
    # r_dot^2 = (ell^2 + 2(e^2/mu) r - 2 h r^2) / (r (rho - r))
    # For ell=0:
    # r_dot^2 = (2 E2_MU r - 2 h r^2) / (r (rho - r))
    #         = (2 E2_MU - 2 h r) / (rho - r)
    # Valid for r != rho
    num = 2 * E2_MU - 2 * h * r
    den = RHO - r
    return num / den
end

println("=" ^ 72)
println("Part 1: Radial velocity analysis at the critical sphere (ell = 0)")
println("=" ^ 72)

# At r = rho, the denominator vanishes. We need the numerator to also
# vanish for r_dot^2 to have a finite limit.
# num(rho) = 2 E2_MU - 2 h rho = 2(1 - h)  [with our units]
# For this to vanish: h = E2_MU / rho = 1 (the "critical energy")
h_critical = E2_MU / RHO
@printf("Critical energy h* = e^2/(mu rho) = %g\n", h_critical)
@printf("At r = rho with h = h*: numerator = %g (both vanish => 0/0)\n",
    2 * E2_MU - 2 * h_critical * RHO)

# L'Hopital: lim_{r->rho} (2 E2_MU - 2h r)/(rho - r) = 2h/1 = 2h
# So r_dot^2 -> 2h at r = rho (when h = h*)
println()
@printf("L'Hopital limit: r_dot^2 -> 2 h* = %g at r = rho\n", 2 * h_critical)
@printf("So |r_dot| -> sqrt(2 h*) = %g (finite!)\n", sqrt(2 * h_critical))

println()
println("Key finding: For the SPECIFIC energy h = h* = e^2/(mu rho) = c^2,")
println("the radial velocity is finite at r = rho. For all other energies,")
println("r_dot^2 diverges (h < h*) or the trajectory does not reach rho (h > h*).")

# ============================================================================
# Part 2: Numerical integration of the radial ODE
# ============================================================================
println()
println("=" ^ 72)
println("Part 2: Numerical integration through the critical radius")
println("=" ^ 72)

# Use the energy equation directly: dr/dt = -sqrt(r_dot^2) for infall
# We need to handle the 0/0 at r = rho carefully.

# Strategy: integrate from r0 < rho (sub-critical) toward r = 0,
# and from r0 > rho (super-critical) toward r = rho.

function integrate_radial_supercritical(r0, h, dt, n_steps)
    # Start at r0 > rho, integrate inward
    r = r0
    t = 0.0
    trajectory = [(t, r)]

    for _ in 1:n_steps
        if r <= RHO * (1 + 1e-12)
            break
        end
        v2 = radial_velocity_squared(r, h)
        if v2 < 0
            println("  Turning point reached at r = $r")
            break
        end
        rdot = -sqrt(v2)  # infall
        r_new = r + rdot * dt
        if r_new < RHO
            # Would cross rho; find the crossing time
            # Linear interpolation
            dt_cross = (r - RHO) / (-rdot)
            t += dt_cross
            r = RHO
            push!(trajectory, (t, r))
            break
        end
        r = r_new
        t += dt
        push!(trajectory, (t, r))
    end
    return trajectory
end

function integrate_radial_subcritical(r0, h, dt, n_steps)
    # Start at r0 < rho, integrate inward toward r = 0
    # In the sub-critical region, r_dot^2 = (2 E2_MU - 2h r)/(rho - r)
    # For r < rho, (rho - r) > 0, and if 2 E2_MU - 2h r > 0, then r_dot^2 > 0.
    r = r0
    t = 0.0
    trajectory = [(t, r)]

    for _ in 1:n_steps
        if r < 1e-10
            println("  Reached r ~ 0 at t = $t")
            break
        end
        v2 = radial_velocity_squared(r, h)
        if v2 < 0
            @printf("  v^2 = %.6e < 0 at r = %.6e -- trajectory forbidden\n", v2, r)
            break
        end
        rdot = -sqrt(v2)  # infall
        r_new = r + rdot * dt
        if r_new < 0
            r_new = 0.0
        end
        r = r_new
        t += dt
        push!(trajectory, (t, r))
    end
    return trajectory
end

# --- Experiment 2a: Super-critical infall at h = h* ---
println()
println("Experiment 2a: Super-critical infall at critical energy h = h*")
r0_super = 2.0
dt_super = 1e-5
traj_super = integrate_radial_supercritical(r0_super, h_critical, dt_super, 10_000_000)
t_end, r_end = traj_super[end]
@printf("  r0 = %.3f, h = %.6f\n", r0_super, h_critical)
@printf("  Final: t = %.6f, r = %.10f\n", t_end, r_end)
@printf("  Reached rho? %s\n", r_end <= RHO * (1 + 1e-6) ? "YES" : "NO")

# --- Experiment 2b: Sub-critical infall at h = h* ---
println()
println("Experiment 2b: Sub-critical infall at critical energy h = h*")
r0_sub = RHO * 0.99  # start just inside rho
dt_sub = 1e-6
traj_sub = integrate_radial_subcritical(r0_sub, h_critical, dt_sub, 10_000_000)
t_end_sub, r_end_sub = traj_sub[end]
@printf("  r0 = %.6f (= 0.99 * rho), h = %.6f\n", r0_sub, h_critical)
@printf("  Final: t = %.6f, r = %.10f\n", t_end_sub, r_end_sub)

# Check velocity near r = 0
if length(traj_sub) > 2
    t2, r2 = traj_sub[end]
    t1, r1 = traj_sub[end-1]
    if t2 > t1
        v_final = abs((r2 - r1) / (t2 - t1))
        @printf("  Final |r_dot| ~ %.6f (theory: sqrt(2) = %.6f)\n", v_final, sqrt(2.0))
    end
end

# --- Experiment 2c: Through-critical trajectory at h = h* ---
println()
println("Experiment 2c: Through-critical matching at h = h*")
println()
println("At r = rho with h = h*:")
println("  Super-critical side: r_dot -> -sqrt(2h*) = -sqrt(2) (infall)")
println("  Sub-critical side:   r_dot -> -sqrt(2h*) = -sqrt(2) (infall continues)")
println()
println("The velocity is CONTINUOUS through rho! The trajectory passes")
println("through the critical sphere with finite speed sqrt(2) * c.")
println()

# Verify by integrating from well above rho, through rho, to r = 0
function integrate_through_critical(r0, h, dt, n_steps)
    # Adaptive integration that handles the 0/0 at r = rho
    r = r0
    t = 0.0
    trajectory = [(t, r)]

    for step in 1:n_steps
        if r < 1e-10
            @printf("  Reached r ~ 0 at t = %.6f, step %d\n", t, step)
            break
        end

        # Near r = rho, use L'Hopital limit
        if abs(r - RHO) < 1e-8
            v2 = 2 * h  # L'Hopital limit
        else
            v2 = radial_velocity_squared(r, h)
        end

        if v2 < 0
            @printf("  v^2 < 0 at r = %.10f, step %d\n", r, step)
            break
        end

        rdot = -sqrt(v2)
        r_new = r + rdot * dt
        if r_new < 0
            r_new = 0.0
        end
        r = r_new
        t += dt
        push!(trajectory, (t, r))
    end
    return trajectory
end

println("Experiment 2c: Full through-critical integration")
r0_through = 1.5
dt_through = 1e-6
traj_through = integrate_through_critical(r0_through, h_critical, dt_through, 10_000_000)
t_end_through, r_end_through = traj_through[end]
@printf("  r0 = %.3f (super-critical), h = h* = %.6f\n", r0_through, h_critical)
@printf("  Final: t = %.6f, r = %.10f\n", t_end_through, r_end_through)

# Track the moment of crossing rho
local crossed = false
for i in 2:length(traj_through)
    t_prev, r_prev = traj_through[i-1]
    t_curr, r_curr = traj_through[i]
    if r_prev > RHO && r_curr <= RHO
        @printf("  Crossed rho at t ~ %.6f\n", t_curr)
        # Check velocity at crossing
        if t_curr > t_prev
            v_cross = abs((r_curr - r_prev) / (t_curr - t_prev))
            @printf("  |r_dot| at crossing ~ %.6f (theory: sqrt(2) = %.6f)\n",
                v_cross, sqrt(2.0))
        end
        crossed = true
        break
    end
end
if !crossed
    println("  Did NOT cross rho")
end

# ============================================================================
# Part 3: Physical interpretation
# ============================================================================
println()
println("=" ^ 72)
println("Part 3: Physical interpretation and Floer-theoretic consequences")
println("=" ^ 72)
println()
println("RESULTS SUMMARY:")
println()
println("1. For ell = 0, the critical sphere r = rho IS traversable, but ONLY")
println("   at the critical energy h* = e^2/(mu rho) = c^2.")
println()
println("2. The effective mass m_eff = mu/(1 - rho/r) passes through infinity")
println("   at r = rho and becomes negative for r < rho. At h = h*, the")
println("   trajectory has r_dot = sqrt(2) c at the crossing -- the particle")
println("   passes through the infinite-mass barrier with finite velocity.")
println()
println("3. This is possible because the MOMENTUM p_r = m_eff * r_dot also")
println("   passes through infinity and changes sign:")
println("   - Super-critical: p_r = mu(1-rho/r)^{-1} * r_dot > 0 (large)")
println("   - At rho: p_r -> +infinity")
println("   - Sub-critical: p_r = mu(1-rho/r)^{-1} * r_dot < 0 (negative mass!)")
println()
println("4. The energy h = h* is special because the numerator and denominator")
println("   of r_dot^2 vanish simultaneously at r = rho, giving a 0/0 that")
println("   resolves to a finite value by L'Hopital.")
println()
println("5. For h != h*, the trajectory either:")
println("   - h < h*: r_dot^2 -> +infinity at rho (momentum diverges)")
println("   - h > h*: r_dot^2 -> -infinity (turning point before rho)")

# ============================================================================
# Part 4: Energy scan -- which energies allow crossing?
# ============================================================================
println()
println("=" ^ 72)
println("Part 4: Energy scan near h*")
println("=" ^ 72)

for delta_h in [-0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1]
    h = h_critical + delta_h
    num_at_rho = 2 * E2_MU - 2 * h * RHO
    @printf("  h = h* + %.3f = %.3f: numerator at rho = %.4f, ",
        delta_h, h, num_at_rho)
    if abs(num_at_rho) < 1e-10
        @printf("0/0 -> r_dot^2 = %.4f (L'Hopital)\n", 2*h)
    elseif num_at_rho > 0
        @printf("r_dot^2 -> +inf (momentum blowup)\n")
    else
        @printf("r_dot^2 -> -inf (turning point)\n")
    end
end

# ============================================================================
# Part 5: Sub-critical ell=0 oscillation (collision bounce verification)
# ============================================================================
println()
println("=" ^ 72)
println("Part 5: Sub-critical ell=0 collision (verification of sqrt(2) speed)")
println("=" ^ 72)

# Start deep in the sub-critical region with a specific energy
r0_deep = 0.5 * RHO
h_deep = E2_MU / r0_deep  # energy such that r_dot = 0 at r0 (turning point)
# Actually: r_dot^2 = (2 - 2h r)/(1 - r) at r0 = 0.5:
#   r_dot^2 = (2 - h) / (0.5) => for r_dot = 0, need 2 - h = 0 => h = 2
h_turn = E2_MU / r0_deep  # = 2.0
@printf("Starting at r0 = %.3f with h = %.3f (turning point energy)\n", r0_deep, h_turn)

v2_check = radial_velocity_squared(r0_deep, h_turn)
@printf("Verification: r_dot^2 at r0 = %.6e (should be ~0)\n", v2_check)

# Give a small inward kick
h_kick = h_turn - 0.1
@printf("Using h = %.3f (slightly below turning point, particle falls inward)\n", h_kick)

traj_deep = integrate_radial_subcritical(r0_deep, h_kick, 1e-6, 5_000_000)
if length(traj_deep) > 10
    t_f, r_f = traj_deep[end]
    t_p, r_p = traj_deep[end-1]
    if t_f > t_p
        v_f = abs((r_f - r_p) / (t_f - t_p))
        @printf("Final r = %.10f, |r_dot| ~ %.6f\n", r_f, v_f)
        @printf("Theory predicts |r_dot| -> sqrt(2) c = %.6f as r -> 0\n", sqrt(2.0) * C)
    end
end

println()
println("=" ^ 72)
println("CONCLUSION")
println("=" ^ 72)
println()
println("For ell = 0 (head-on) like-charge Weber dynamics:")
println()
println("(A) The sub-critical region r < rho has well-defined dynamics with")
println("    negative effective mass. Collisions at r = 0 arrive at finite")
println("    speed sqrt(2) c and are C^0-continuable (collision bounce).")
println()
println("(B) The critical sphere r = rho can be traversed radially ONLY at the")
println("    specific energy h* = c^2. At this energy, the 0/0 indeterminacy")
println("    resolves to r_dot^2 = 2c^2, giving continuous passage.")
println()
println("(C) For Floer theory: the critical sphere is NOT a 'stop' in the")
println("    Ganatra-Pardon-Shende sense (it is codimension 1, not 2). It is")
println("    a degenerate boundary where the kinetic metric changes signature.")
println("    The ell = 0 sector provides a 1-parameter family of orbits that")
println("    can cross this boundary (at h = h*), suggesting the critical sphere")
println("    is 'permeable' in a measure-zero sense, not impermeable.")
println()
println("(D) Practical regularization: for ell = 0, the existing collision")
println("    bounce at r ~ 0 handles the sub-critical dynamics. Adding a")
println("    'critical sphere bounce' (reflecting r -> 2*rho - r at r = rho)")
println("    would handle super-critical orbits that reach rho, but this is")
println("    only needed for the measure-zero set at h = h*.")
