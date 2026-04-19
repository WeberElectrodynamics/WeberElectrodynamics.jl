# Scan R at the 3D resonance-tongue optimum to locate the full peak.
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

results = []
for R in (3.0, 3.5, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 5.0, 5.5)
    q0, p0, m, ch = ic_sym_orbiter_3d(; r_pp=4.0, R=R, orb=1.30, z_kick=0.13)
    sol = SharedSurvey.run(q0, p0, m, ch; tmax=2000.0, dt=1e-3, c=1.0, dims=3)
    s = SharedSurvey.summarize(sol)
    @printf("R=%.2f  rc=%s  tf=%.2f  drift=%.2e %%\n", R, s.retcode, s.t_final, s.E_drift_pct)
    push!(results, (R, string(s.retcode), s.t_final, s.E_drift_pct))
end

open(joinpath(@__DIR__, "sym_orbiter_3d_R_scan.log"), "w") do io
    println(io, "R,retcode,t_final,drift_pct")
    for r in results; println(io, join(r, ",")); end
end
