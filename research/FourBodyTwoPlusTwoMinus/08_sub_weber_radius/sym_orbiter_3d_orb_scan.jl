# Push the 3D record further: vary orb near 1.3 at (rpp=R=4, z_kick=0.13).
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

results = []
for orb in (1.15, 1.20, 1.25, 1.28, 1.30, 1.32, 1.35, 1.40, 1.45, 1.50)
    q0, p0, m, ch = ic_sym_orbiter_3d(; r_pp=4.0, R=4.0, orb=orb, z_kick=0.13)
    sol = SharedSurvey.run(q0, p0, m, ch; tmax=2000.0, dt=1e-3, c=1.0, dims=3)
    s = SharedSurvey.summarize(sol)
    @printf("orb=%.2f  rc=%s  tf=%.2f  drift=%.2e %%\n", orb, s.retcode, s.t_final, s.E_drift_pct)
    push!(results, (orb, string(s.retcode), s.t_final, s.E_drift_pct))
end

open(joinpath(@__DIR__, "sym_orbiter_3d_orb_scan.log"), "w") do io
    println(io, "orb,retcode,t_final,drift_pct")
    for r in results; println(io, join(r, ",")); end
end
