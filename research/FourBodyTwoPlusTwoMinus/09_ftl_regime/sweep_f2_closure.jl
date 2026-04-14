# Autonomous follow-up to REPORT §16 item 3: shoot for a closed Weber loop
# starting from IC-F2 (fast dimers). Scan (v_target, dyad_len) and measure
# phase-space return distance to detect near-periodic closure.
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "ftl_experiments.jl"))
using .SharedSurvey

function return_metric(sol; t_skip = 1.0)
    ts = sol.t
    n = length(ts)
    n < 20 && return (Inf, 0.0)
    z0 = vcat(sol.q[1], sol.p[1])
    best = Inf; tbest = 0.0
    for k in 2:n
        ts[k] < t_skip && continue
        zk = vcat(sol.q[k], sol.p[k])
        d = norm(zk - z0)
        if d < best
            best = d; tbest = ts[k]
        end
    end
    return (best, tbest)
end

function rdot_13(sol)
    rdots = Float64[]
    for k in 1:length(sol.t)
        q = sol.q[k]; p = sol.p[k]
        r13 = [q[5]-q[1], q[6]-q[2]]
        v1 = [p[1], p[2]]; v3 = [p[5], p[6]]
        push!(rdots, dot(r13, v3 - v1) / norm(r13))
    end
    return rdots
end

results = []
for v in (-1.90, -1.75, -1.55, 1.55, 1.75, 1.90)
    for a in (0.30, 0.40, 0.50)
        q0, p0, m, ch = ic_f2_fast_dimers(dyad_len=a, sep=2.5, v_target=v)
        sol = try
            SharedSurvey.run(q0, p0, m, ch; tmax=20.0, dt=5e-4, c=1.0)
        catch e
            nothing
        end
        if sol === nothing
            push!(results, (v, a, :Err, 0.0, Inf, 0.0, 0.0))
            continue
        end
        s = SharedSurvey.summarize(sol)
        d, tb = return_metric(sol; t_skip=1.5)
        rd = rdot_13(sol)
        min_rd = minimum(rd); max_rd = maximum(rd)
        push!(results, (v, a, s.retcode, s.t_final, d, tb, s.E_drift_pct))
        @printf("v=%.2f a=%.2f  rc=%s  tf=%.3f  ret=%.4f@t=%.3f  drift=%.2e  rd∈[%.2f,%.2f]\n",
                v, a, s.retcode, s.t_final, d, tb, s.E_drift_pct, min_rd, max_rd)
    end
end

# Save
open(joinpath(@__DIR__, "sweep_f2_closure.log"), "w") do io
    println(io, "v,a,retcode,t_final,ret_dist,t_ret,drift_pct")
    for r in results
        println(io, join(r, ","))
    end
end
println("done")
