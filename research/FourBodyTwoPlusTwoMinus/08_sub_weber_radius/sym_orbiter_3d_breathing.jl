# Test whether the "(+,+) at rest" condition is essential to the symmetric
# double-orbiter bound state. Add a symmetric radial breathing velocity to
# the positive pair (preserves Klein-four, σ_x preserved since v_1 = -v_2).
using Printf
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

# Extend IC to include breathing velocity on (+,+) pair.
function ic_breathing_orbiter(; r_pp=4.0, R=4.0, orb=1.30, z_kick=0.13, vbr=0.0, m=1.0)
    q0, p0, masses, charges = ic_sym_orbiter_3d(; r_pp=r_pp, R=R, orb=orb, z_kick=z_kick, m=m)
    # Particles 1,2 are on ±x. Add radial breathing: p_1x = -m*vbr, p_2x = +m*vbr (contracting)
    p0[1] += -m * vbr   # particle 1 x-momentum (index 1 in 3D layout x,y,z per particle)
    p0[4] +=  m * vbr   # particle 2 x-momentum (index 4 = 3+1)
    return (q0, p0, masses, charges)
end

for vbr in (-0.3, -0.15, -0.05, 0.0, 0.05, 0.15, 0.30)
    q0, p0, m, ch = ic_breathing_orbiter(; vbr=vbr)
    sol = SharedSurvey.run(q0, p0, m, ch; tmax=1500.0, dt=1e-3, c=1.0, dims=3)
    s = SharedSurvey.summarize(sol)
    @printf("vbr=%+.2f  rc=%s  tf=%.2f  drift=%.2e %%\n",
            vbr, s.retcode, s.t_final, s.E_drift_pct)
end
