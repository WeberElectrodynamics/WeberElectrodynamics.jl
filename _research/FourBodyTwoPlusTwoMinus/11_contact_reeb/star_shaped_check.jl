"""
Star-shapedness check for the 4-body 2+/2− Weber Hill region {V <= E_sq}
relative to the alternating-square configuration q*.

For 50 random rays q* + t*v (||v||=1), we step t outward and inward and
test whether V(q* + t v) is monotone (i.e. V crosses E_sq at most once
between the star-point and infinity in each direction).
"""

using Random, LinearAlgebra, Printf
using WeberElectrodynamics

include(joinpath(@__DIR__, "..", "shared", "ic_generators.jl"))

const N = 4
const D = 2
sys = HamiltonianSystem(N, D)

q0, p0, masses, charges = alternating_square(side = 1.0, q = 1.0, m = 1.0,
                                              energy_fraction = 0.0,
                                              velocity_mode = :static)
c = 1.0
kappas = ones(N*(N-1) ÷ 2)
params = vcat(masses, charges, c, kappas)

# V(q) = H(q, p=0)
zero_p = zeros(length(q0))
V(q) = sys.hamiltonian_compiled(q, zero_p, params)

E_sq = V(q0)
@info "Star point V(q*) = $E_sq"

Random.seed!(20260414)
n_rays = 50
t_max = 50.0
n_samples = 400
ray_results = Tuple{Bool,Bool,Float64}[]  # (outward_monotone, inward_finite, max_V_along_ray)

for k in 1:n_rays
    v = randn(length(q0))
    v ./= norm(v)
    ts_out = range(1e-3, t_max; length = n_samples)
    Vs_out = [V(q0 .+ t .* v) for t in ts_out]
    # Outward: starting from V(q*) = E_sq, V should rise monotonically
    # (Hill region {V <= E} contains q* on its boundary -- we are AT level E).
    # Star-shaped test: for star-shapedness of {V <= E} w.r.t. q*, every
    # outward ray must satisfy V(q*+tv) >= E for all t >= 0 (q* is a
    # boundary star-point) OR equivalently, ∂V/∂t·t > 0 along the ray.
    diffs = diff(Vs_out)
    outward_monotone = all(diffs .>= -1e-10) || all(diffs .<= 1e-10)
    # We actually want: V is non-decreasing past q* in every direction
    # (since q* is a critical point of V at this energy).
    nondec = all(diffs .>= -1e-10)
    push!(ray_results, (nondec, true, maximum(Vs_out)))
end

frac = count(r -> r[1], ray_results) / n_rays
@printf "Fraction of monotonic-outward rays: %.3f (%d / %d)\n" frac count(r->r[1], ray_results) n_rays

# Additionally count how many rays at least cross the level set only once.
global single_cross = 0
for k in 1:n_rays
    global single_cross
    v = randn(length(q0)); v ./= norm(v)
    ts = range(1e-4, t_max; length = n_samples)
    Vs = [V(q0 .+ t .* v) for t in ts]
    crossings = 0
    for i in 1:(length(Vs)-1)
        if (Vs[i] - E_sq) * (Vs[i+1] - E_sq) < 0
            crossings += 1
        end
    end
    if crossings <= 1
        single_cross += 1
    end
end
@printf "Fraction of rays with <=1 level crossing: %.3f\n" (single_cross / n_rays)
