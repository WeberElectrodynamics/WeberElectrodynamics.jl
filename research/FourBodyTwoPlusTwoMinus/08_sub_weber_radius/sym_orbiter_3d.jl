# Autonomous follow-up to REPORT §16 item 7: promote the symmetric double-orbiter
# family (the strongest 2D bound-state candidate) to 3D and check survival.
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
using .SharedSurvey

# Symmetric double-orbiter in 3D, motion confined to xy plane initially.
# + at (±r_pp/2, 0, 0); - at (0, ±R, 0) with ±v_orb x̂ tangential velocity.
# With purely planar IC the dynamics stay planar → sanity baseline.
# Then add a small out-of-plane z-kick to test transverse stability.
function ic_sym_orbiter_3d(; r_pp=3.0, R=3.0, orb=1.0, z_kick=0.0, m=1.0)
    X = [
         r_pp/2  0.0  0.0;
        -r_pp/2  0.0  0.0;
         0.0      R    0.0;
         0.0     -R    0.0;
    ]
    v_orb = orb * sqrt(2/R)
    P = zeros(4, 3)
    P[3, 1] =  m * v_orb
    P[4, 1] = -m * v_orb
    # transverse z-kick on + pair, symmetric (keeps p_total = 0)
    P[1, 3] =  m * z_kick
    P[2, 3] = -m * z_kick
    masses = fill(m, 4)
    charges = [1.0, 1.0, -1.0, -1.0]
    q0 = vec(collect(X'))
    p0 = vec(collect(P'))
    return (q0, p0, masses, charges)
end

results = []
for z in (0.0, 0.05, 0.10, 0.20)
    for (rpp, R, orb) in [(3.0,3.0,1.3), (4.0,4.0,1.3)]
        q0, p0, m, ch = ic_sym_orbiter_3d(; r_pp=rpp, R=R, orb=orb, z_kick=z)
        sol = try
            SharedSurvey.run(q0, p0, m, ch; tmax=1000.0, dt=1e-3, c=1.0, dims=3)
        catch e
            nothing
        end
        if sol === nothing
            @printf("z=%.2f rpp=%.1f R=%.1f orb=%.1f  ERR\n", z, rpp, R, orb)
            continue
        end
        s = SharedSurvey.summarize(sol)
        @printf("z=%.2f rpp=%.1f R=%.1f orb=%.1f  rc=%s  tf=%.2f  drift=%.2e %%\n",
                z, rpp, R, orb, s.retcode, s.t_final, s.E_drift_pct)
        push!(results, (z, rpp, R, orb, string(s.retcode), s.t_final, s.E_drift_pct))
    end
end

open(joinpath(@__DIR__, "sym_orbiter_3d.log"), "w") do io
    println(io, "z_kick,r_pp,R,orb,retcode,t_final,drift_pct")
    for r in results; println(io, join(r, ",")); end
end
println("done")
