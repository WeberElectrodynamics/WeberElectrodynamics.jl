# Locate the 3D survival optimum by fine-scanning z_kick near 0.10 at (4,4,1.3).
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

results = []
for z in (0.06, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.17)
    q0, p0, m, ch = ic_sym_orbiter_3d(; r_pp=4.0, R=4.0, orb=1.3, z_kick=z)
    sol = SharedSurvey.run(q0, p0, m, ch; tmax=1500.0, dt=1e-3, c=1.0, dims=3)
    s = SharedSurvey.summarize(sol)
    @printf("z=%.3f  rc=%s  tf=%.2f  drift=%.2e %%\n", z, s.retcode, s.t_final, s.E_drift_pct)
    push!(results, (z, string(s.retcode), s.t_final, s.E_drift_pct))
end

open(joinpath(@__DIR__, "sym_orbiter_3d_optimize.log"), "w") do io
    println(io, "z_kick,retcode,t_final,drift_pct")
    for r in results; println(io, join(r, ",")); end
end
println("done")
