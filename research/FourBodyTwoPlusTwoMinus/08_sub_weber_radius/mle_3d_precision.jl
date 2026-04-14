# High-precision MLE test at the 3D optimum to discriminate genuine KAM torus
# from thin Arnold layer. Vary epsilon by 4 decades and increase N.
using LinearAlgebra, Random, Printf
include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
include(joinpath(@__DIR__, "sym_orbiter_3d.jl"))
using .SharedSurvey

function mle(q0, p0, masses, charges; dims, tmax, dt_renorm, dt, eps, seed)
    n_q = length(q0); D = 2*n_q
    N = Int(round(tmax / dt_renorm))
    rng = MersenneTwister(seed)
    u = randn(rng, D); u ./= norm(u)
    x_ref = vcat(q0, p0); x_sh = x_ref .+ eps .* u
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
    return (lambda = log_sum / (completed * dt_renorm), completed=completed, N=N)
end

ic = ic_sym_orbiter_3d(; r_pp=4.0, R=4.0, orb=1.3, z_kick=0.13)
q0, p0, m, ch = ic

for eps in (1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
    for seed in (17, 31, 47)
        res = mle(q0, p0, m, ch; dims=3, tmax=300.0, dt_renorm=0.3, dt=5e-4, eps=eps, seed=seed)
        @printf("ε=%.0e seed=%d  λ=%.5f  N=%d/%d\n", eps, seed, res.lambda, res.completed, res.N)
    end
end
