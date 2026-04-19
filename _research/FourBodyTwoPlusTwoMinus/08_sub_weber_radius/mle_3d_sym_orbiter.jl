# Measure maximum Lyapunov exponent on the 3D symmetric double-orbiter,
# the strongest bound-state candidate in the 14-agent study.
using LinearAlgebra, Random, Printf
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

function mle(q0, p0, masses, charges; dims, tmax=200.0, dt_renorm=0.5,
             dt=1e-3, eps=1e-8, rng=MersenneTwister(17))
    n_q = length(q0)
    D = 2 * n_q
    N = Int(round(tmax / dt_renorm))
    u = randn(rng, D); u ./= norm(u)
    x_ref = vcat(q0, p0)
    x_sh  = x_ref .+ eps .* u
    log_sum = 0.0; completed = 0
    for k in 1:N
        qr, pr = x_ref[1:n_q], x_ref[n_q+1:end]
        qs, ps = x_sh[1:n_q],  x_sh[n_q+1:end]
        sr = SharedSurvey.run(qr, pr, masses, charges; tmax=dt_renorm, dt=dt, c=1.0, dims=dims)
        ss = SharedSurvey.run(qs, ps, masses, charges; tmax=dt_renorm, dt=dt, c=1.0, dims=dims)
        (sr.retcode != :Success || ss.retcode != :Success) && break
        x_ref = vcat(sr.q[end], sr.p[end])
        x_sh  = vcat(ss.q[end], ss.p[end])
        δ = norm(x_sh - x_ref)
        (!isfinite(δ) || δ <= 0) && break
        log_sum += log(δ/eps)
        x_sh = x_ref .+ eps .* (x_sh .- x_ref) ./ δ
        completed = k
    end
    return (lambda = log_sum / (completed * dt_renorm), completed = completed)
end

# Flagship 3D configurations
cases = [
    ("3D planar (3,3,1.0) z=0.00",  ic_sym_orbiter_3d(; r_pp=3.0, R=3.0, orb=1.0, z_kick=0.0)),
    ("3D z-kicked (3,3,1.0) z=0.10", ic_sym_orbiter_3d(; r_pp=3.0, R=3.0, orb=1.0, z_kick=0.10)),
    ("3D optimum (4,4,1.3) z=0.13",  ic_sym_orbiter_3d(; r_pp=4.0, R=4.0, orb=1.3, z_kick=0.13)),
]

for (name, ic) in cases
    q0, p0, m, ch = ic
    res = mle(q0, p0, m, ch; dims=3, tmax=200.0, dt_renorm=0.5)
    @printf("%-32s  λ_max = %.4f   (completed %d intervals)\n",
            name, res.lambda, res.completed)
end
