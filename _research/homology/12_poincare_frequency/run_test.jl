using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using LinearAlgebra

sys = HamiltonianSystem(4, 2)
a = 1.5; b = 1.45; eta = 0.75
X = [a 0.0; -a 0.0; 0.0 b; 0.0 -b]
masses = [1.0, 1.0, 1.0, 1.0]
charges = [1.0, 1.0, -1.0, -1.0]
U = 0.0
for i in 1:4, j in (i+1):4
    r = norm(X[i, :] - X[j, :])
    global U += charges[i] * charges[j] / r
end
T_target = eta * abs(U)
speed = sqrt(2 * T_target / 4)
P = zeros(4, 2)
for i in 1:4
    r = X[i, :]
    t = [-r[2], r[1]] / max(norm(r), eps())
    P[i, :] = speed * t
end
X .-= sum(X, dims = 1) / 4
P .-= sum(P, dims = 1) / 4
q0 = vec(transpose(X))
p0 = vec(transpose(P))

# Quick timing test
t1 = time()
prob = HamiltonianProblem(sys, (0.0, 50.0), q0, p0;
    masses = masses, charges = charges, c = 1.0, dt = 5e-4)
sol = solve(prob, SymmetricProjectionIntegrator(); callbacks = CollisionBounce(0.02))
t2 = time()
println("tmax=50, dt=5e-4: $(round(t2-t1, digits=2))s, retcode=$(sol.retcode), t_final=$(sol.t[end]), nsteps=$(length(sol.t))")

# Test tmax=200
t1 = time()
prob2 = HamiltonianProblem(sys, (0.0, 200.0), q0, p0;
    masses = masses, charges = charges, c = 1.0, dt = 5e-4)
sol2 = solve(prob2, SymmetricProjectionIntegrator(); callbacks = CollisionBounce(0.02))
t2 = time()
println("tmax=200, dt=5e-4: $(round(t2-t1, digits=2))s, retcode=$(sol2.retcode), t_final=$(sol2.t[end]), nsteps=$(length(sol2.t))")
