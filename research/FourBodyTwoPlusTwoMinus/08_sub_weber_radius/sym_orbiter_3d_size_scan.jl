# Diagonal scan R=rpp at the optimum orb and z_kick. Does larger "square size"
# give longer horizons?
using Printf
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

for L in (3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0)
    q0, p0, m, ch = ic_sym_orbiter_3d(; r_pp=L, R=L, orb=1.30, z_kick=0.13)
    sol = SharedSurvey.run(q0, p0, m, ch; tmax=2000.0, dt=1e-3, c=1.0, dims=3)
    s = SharedSurvey.summarize(sol)
    @printf("L=%.2f  rc=%s  tf=%.2f  drift=%.2e %%\n", L, s.retcode, s.t_final, s.E_drift_pct)
end
