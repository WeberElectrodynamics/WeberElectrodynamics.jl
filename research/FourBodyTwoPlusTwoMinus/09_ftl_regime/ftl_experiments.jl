"""
Agent 9 — Faster-Than-Light Relative Velocity Regime (4-body 2+/2−).

Builds initial conditions where one or more pair radial velocities |ṙ_ij|
exceed c (and in some cases √2·c), runs a short integration, and plots
|ṙ_ij|(t) for all 6 pairs together with the c and √2·c thresholds.

Run from project root:
    julia --project=. research/FourBodyTwoPlusTwoMinus/09_ftl_regime/ftl_experiments.jl
"""

using LinearAlgebra
using Printf
using WeberElectrodynamics
using Plots

const HERE = @__DIR__
const SHARED = joinpath(HERE, "..", "shared")
include(joinpath(SHARED, "ic_generators.jl"))
include(joinpath(SHARED, "metrics.jl"))
include(joinpath(SHARED, "run_survey.jl"))
using .SharedSurvey

const FIGDIR = joinpath(HERE, "figures")
isdir(FIGDIR) || mkpath(FIGDIR)

const C_SPEED = 1.0      # natural units
const SQRT2C  = sqrt(2) * C_SPEED

# --------------------------------------------------------------------------
# IC builders. Charges layout: [+q, +q, -q, -q]; particles 1,2 positive.
# We construct positions then assign momenta so that the relevant pair
# radial velocity ṙ_ij = (x_i - x_j)·(v_i - v_j)/r equals a target value.
# --------------------------------------------------------------------------

"""
IC-F1: alternating square (+,+,-,-) with breathing radial velocity targeted
to ṙ ≈ v_target on all four edge pairs simultaneously. We achieve this
with a uniform inward radial velocity per particle of magnitude v_target/2,
since for any opposing pair the relative radial velocity is 2·(v_radial).
"""
function ic_f1_breathing_square(; side = 1.5, q = 1.0, m = 1.0, v_target = 1.5)
    s = side / 2
    # vertices: 1(+,+y), 2(+,-y), 3(-,-x at +x), 4(-,+x at -x)? Use an
    # alternating layout — opposite corners share charge sign so all 4 edges
    # are unlike pairs.
    X = [
         s   s;   # +
        -s  -s;   # +
         s  -s;   # -
        -s   s;   # -
    ]
    masses  = fill(m, 4)
    charges = [q, q, -q, -q]
    # inward radial velocity v_in per particle => for any pair with opposite
    # position vectors, ṙ_rel = -2*v_in. For adjacent pairs, ṙ_rel = -v_in*√2.
    v_in = v_target / 2
    P = zeros(4, 2)
    for i = 1:4
        r = X[i, :]
        P[i, :] = -m * v_in * r / norm(r)
    end
    # Already symmetric, COM=0, total p=0
    return (_flatten(X), _flatten(P), masses, charges)
end

"""
IC-F2: two (+,−) dimers with a large intra-dyad radial velocity > √2·c.
Particles 1,3 form right dyad (+ at top, − at bottom of right column);
particles 2,4 form left dyad. We give each + an outward y push and each −
a matching inward push so ṙ for both dimers is +v_target.
"""
function ic_f2_fast_dimers(; dyad_len = 0.4, sep = 2.5, q = 1.0, m = 1.0,
                           v_target = 1.8)
    a = dyad_len
    R = sep
    X = [
         R/2   a/2;   # +1
        -R/2   a/2;   # +2
         R/2  -a/2;   # -1
        -R/2  -a/2;   # -2
    ]
    masses  = fill(m, 4)
    charges = [q, q, -q, -q]
    # For pair (1,3): r = (0, a), so ṙ = (v1y - v3y). Set v1y = +v_target/2,
    # v3y = -v_target/2. Same for (2,4). Total p_y = 0 by symmetry.
    P = zeros(4, 2)
    P[1, 2] =  m * v_target / 2
    P[3, 2] = -m * v_target / 2
    P[2, 2] =  m * v_target / 2
    P[4, 2] = -m * v_target / 2
    return (_flatten(X), _flatten(P), masses, charges)
end

"""
IC-F3: one (+,+) like-pair receding with ṙ ≈ 2c, the two negatives stationary
far away. Pair (1,2) is the like pair; 3,4 are placed off-axis at large
separation with zero momentum.
"""
function ic_f3_fast_like_pair(; pair_sep = 1.0, far = 5.0, q = 1.0, m = 1.0,
                              v_target = 2.0)
    X = [
        -pair_sep/2  0.0;   # +
         pair_sep/2  0.0;   # +
         0.0         far;   # -
         0.0        -far;   # -
    ]
    masses  = fill(m, 4)
    charges = [q, q, -q, -q]
    # Pair (1,2): r = (pair_sep, 0). Set v1x = -v_target/2, v2x = +v_target/2
    # so ṙ = +v_target. Recoil on negatives to zero total p (none needed,
    # already symmetric).
    P = zeros(4, 2)
    P[1, 1] = -m * v_target / 2
    P[2, 1] =  m * v_target / 2
    return (_flatten(X), _flatten(P), masses, charges)
end

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
const PAIR_LABELS = ["(1,2) ++", "(1,3) +-", "(1,4) +-",
                     "(2,3) +-", "(2,4) +-", "(3,4) --"]

function compute_pair_rdot_series(sol, masses)
    n, d = 4, 2
    nt = length(sol.t)
    out = zeros(nt, 6)
    for k = 1:nt
        rdots = pair_rdots(sol.q[k], sol.p[k], masses, n, d)
        out[k, :] = rdots
    end
    return out
end

function plot_rdot(sol, masses, title_str, fname)
    rdots = compute_pair_rdot_series(sol, masses)
    p = plot(size=(900, 520), xlabel="t", ylabel="|ṙ_ij|",
             title=title_str, legend=:outerright)
    for k = 1:6
        plot!(p, sol.t, abs.(rdots[:, k]); label=PAIR_LABELS[k], lw=1.4)
    end
    hline!(p, [C_SPEED];   ls=:dash, color=:black, label="c")
    hline!(p, [SQRT2C];    ls=:dash, color=:red,   label="√2·c")
    savefig(p, joinpath(FIGDIR, fname))
    return p, rdots
end

function report_initial(label, q0, p0, masses, charges, c)
    println("\n=== $label ===")
    rdots0 = pair_rdots(q0, p0, masses, 4, 2)
    rs0    = pair_distances(q0, 4, 2)
    rho    = critical_radii(charges, masses, c)
    rho_pairs = [rho[1,2], rho[1,3], rho[1,4], rho[2,3], rho[2,4], rho[3,4]]
    for k = 1:6
        weber_factor = 1 - rdots0[k]^2 / (2 * c^2)
        @printf("  pair %s  r=%.3f  ṙ=%+.3f  (1-ṙ²/2c²)=%+.3f  ρ=%+.3f\n",
                PAIR_LABELS[k], rs0[k], rdots0[k], weber_factor, rho_pairs[k])
    end
end

# --------------------------------------------------------------------------
# Run the three experiments
# --------------------------------------------------------------------------
function run_all()
    tmax = 5.0
    dt   = 1e-4
    c    = C_SPEED
    results = Dict{String,Any}()

    cases = [
        ("IC-F1 breathing square (v≈1.5c)",
         "ftl_f1_breathing_square.png",
         ic_f1_breathing_square(v_target = 1.5)),
        ("IC-F2 fast dimers (v≈1.8c)",
         "ftl_f2_fast_dimers.png",
         ic_f2_fast_dimers(v_target = 1.8)),
        ("IC-F3 fast (++) pair (v≈2.0c)",
         "ftl_f3_fast_like_pair.png",
         ic_f3_fast_like_pair(v_target = 2.0)),
    ]

    for (label, fname, ic) in cases
        q0, p0, masses, charges = ic
        report_initial(label, q0, p0, masses, charges, c)
        sol = SharedSurvey.run(q0, p0, masses, charges;
                               tmax = tmax, dt = dt, c = c, dims = 2)
        s = SharedSurvey.summarize(sol)
        @printf("  retcode=%s  steps=%d  t_final=%.3f  E_drift_max=%.3e %%\n",
                s.retcode, s.n_steps, s.t_final, s.E_drift_pct)
        _, rdots = plot_rdot(sol, masses, label, fname)
        # Threshold-crossing diagnostics
        absmax = maximum(abs, rdots)
        cross_c     = any(any(abs.(rdots[:, k]) .> c)         for k = 1:6)
        cross_s2c   = any(any(abs.(rdots[:, k]) .> SQRT2C)    for k = 1:6)
        passes_s2c  = false
        for k = 1:6
            v = abs.(rdots[:, k])
            for i = 2:length(v)
                if (v[i-1] - SQRT2C) * (v[i] - SQRT2C) < 0
                    passes_s2c = true
                    break
                end
            end
            passes_s2c && break
        end
        @printf("  |ṙ|_max=%.3f   any>c:%s   any>√2c:%s   crosses √2c smoothly:%s\n",
                absmax, cross_c, cross_s2c, passes_s2c)
        results[label] = (sol = sol, summary = s, rdots = rdots,
                          absmax = absmax, crosses_sqrt2c = passes_s2c)
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all()
end
