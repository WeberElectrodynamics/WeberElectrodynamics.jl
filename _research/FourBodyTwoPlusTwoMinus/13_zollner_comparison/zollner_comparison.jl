"""
Agent 13 — Zöllner-enabled comparison.

Re-runs a small set of ICs flagged by Agents 5, 6, 8 with
`zollner_enabled=true` and `a ∈ {0.0, 0.1, 0.5, 1.0}`, and writes
`zollner_comparison.csv` summarizing retcode, final time, energy drift,
and a bound/unbound flag.

Budget: 5 ICs × 4 a-values = 20 integrations, tmax=30, dt=1e-3.
"""

using LinearAlgebra
using Printf

const SHARED = joinpath(@__DIR__, "..", "shared")
include(joinpath(SHARED, "ic_generators.jl"))
include(joinpath(SHARED, "run_survey.jl"))
using .SharedSurvey

using WeberElectrodynamics

# ---------------------------------------------------------------------------
# IC builders (return NamedTuple with q, p, m, q_, label)
# ---------------------------------------------------------------------------

function ic_F1b_breathing_square()
    # Agent 5's only closed periodic orbit: alternating square side=2.0,
    # outward radial velocity 0.5. Not the shared generator breathing mode
    # (which uses energy_fraction), so we build it directly.
    s = 2.0
    vrad = 0.5
    X = [
         s/2  s/2;
        -s/2 -s/2;
        -s/2  s/2;
         s/2 -s/2;
    ]
    masses = fill(1.0, 4)
    charges = [1.0, 1.0, -1.0, -1.0]
    # Outward radial momentum per particle.
    P = zeros(4, 2)
    for i in 1:4
        r = X[i, :]
        P[i, :] = vrad * r / norm(r)
    end
    # COM / total-p zero by symmetry; flatten.
    q0 = vec(permutedims(X))
    p0 = vec(permutedims(P))
    return (q=q0, p=p0, m=masses, q_=charges, label="F1b_outbreath_s2.00_v0.50")
end

function ic_rhombus_075()
    q, p, m, c = rhombus(a=1.5, b=1.15, energy_fraction=0.75, velocity_mode=:rotating)
    return (q=q, p=p, m=m, q_=c, label="rhombus_a1.5_b1.15_eta0.75")
end

function ic_square_025()
    q, p, m, c = alternating_square(side=1.0, energy_fraction=0.25, velocity_mode=:rotating)
    return (q=q, p=p, m=m, q_=c, label="alt_square_s1.0_eta0.25")
end

function ic_two_dimers()
    q, p, m, c = two_dimers(dyad_length=0.2, separation=2.0, intra_fraction=0.3, inter_velocity=0.05)
    return (q=q, p=p, m=m, q_=c, label="two_dimers_a0.2_R2.0_intra0.3")
end

function ic_nucleus_orbiters()
    # Agent 8 Experiment A: (+,+) tight nucleus at ±0.25 with head-on momenta
    # ±0.2; (−,−) at (±2.5, +5) with small tangential momenta ±0.1.
    X = [
        -0.25  0.0;
         0.25  0.0;
        -2.5   5.0;
         2.5   5.0;
    ]
    P = [
         0.2   0.0;
        -0.2   0.0;
         0.1   0.0;
        -0.1   0.0;
    ]
    masses = fill(1.0, 4)
    charges = [1.0, 1.0, -1.0, -1.0]
    # COM / total-p not exactly zero but fine; HamiltonianProblem does not require it.
    q0 = vec(permutedims(X))
    p0 = vec(permutedims(P))
    return (q=q0, p=p0, m=masses, q_=charges, label="nucleus_plus_orbiters_A")
end

const ICS = [
    ic_F1b_breathing_square(),
    ic_rhombus_075(),
    ic_square_025(),
    ic_two_dimers(),
    ic_nucleus_orbiters(),
]

const A_VALUES = [0.0, 0.1, 0.5, 1.0]
const TMAX = 30.0
const DT = 1e-3

# ---------------------------------------------------------------------------
# Boundedness check: max pairwise distance over trajectory ≤ escape_radius.
# ---------------------------------------------------------------------------
function is_bound(sol, escape_radius::Float64 = 20.0)
    n = 4
    d = 2
    maxr = 0.0
    for k in 1:length(sol.t)
        q = sol.q[k]
        for i in 1:n, j in (i+1):n
            dx = q[(i-1)*d+1] - q[(j-1)*d+1]
            dy = q[(i-1)*d+2] - q[(j-1)*d+2]
            r = sqrt(dx*dx + dy*dy)
            maxr = max(maxr, r)
        end
    end
    return maxr < escape_radius, maxr
end

# ---------------------------------------------------------------------------
# Sweep.
# ---------------------------------------------------------------------------
rows = Vector{NamedTuple}()

for ic in ICS
    for a_val in A_VALUES
        enabled = a_val > 0.0
        @printf("▶ %-36s  a=%.2f (enabled=%s) ... ", ic.label, a_val, enabled)
        local sol
        local retcode, t_final, drift, bound, rmax
        try
            sol = SharedSurvey.run(
                ic.q, ic.p, ic.m, ic.q_;
                tmax=TMAX, dt=DT, c=1.0, dims=2,
                zollner_enabled=enabled, zollner_a=a_val,
            )
            summ = SharedSurvey.summarize(sol)
            retcode = string(summ.retcode)
            t_final = summ.t_final
            drift = summ.E_drift_pct
            bound, rmax = is_bound(sol)
        catch err
            retcode = "Error:" * sprint(showerror, err)[1:min(60, end)]
            t_final = NaN
            drift = NaN
            bound = false
            rmax = NaN
        end
        @printf("rc=%s t=%.3f drift=%.3g bound=%s rmax=%.3g\n",
                retcode, t_final, drift, bound, rmax)
        push!(rows, (
            ic_name = ic.label,
            a = a_val,
            retcode = retcode,
            t_final = t_final,
            E_drift_pct = drift,
            bound = bound,
            r_max = rmax,
        ))
    end
end

# ---------------------------------------------------------------------------
# CSV output.
# ---------------------------------------------------------------------------
const OUT_CSV = joinpath(@__DIR__, "zollner_comparison.csv")
open(OUT_CSV, "w") do io
    println(io, "ic_name,a,retcode,t_final,E_drift_pct,bound,r_max")
    for r in rows
        @printf(io, "%s,%.2f,%s,%.6f,%.6g,%s,%.6g\n",
                r.ic_name, r.a, r.retcode, r.t_final, r.E_drift_pct,
                r.bound, r.r_max)
    end
end

println("\nWrote ", OUT_CSV, " (", length(rows), " rows)")
