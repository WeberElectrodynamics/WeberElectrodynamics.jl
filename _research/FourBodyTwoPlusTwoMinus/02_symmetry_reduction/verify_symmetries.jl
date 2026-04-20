#!/usr/bin/env julia
# Numerical verification that the listed discrete symmetries leave the
# 2D 4-body Weber Hamiltonian invariant when m1=m2=m3=m4=m, q1=q2=+q, q3=q4=-q.
#
# Strategy: use the compiled Hamiltonian. Sample random (q,p) configurations,
# apply the symmetry transformation, evaluate H on both, and check equality
# to machine precision. For polynomial / rational H this is a definitive test.
#
# Run:  julia --project=. research/FourBodyTwoPlusTwoMinus/02_symmetry_reduction/verify_symmetries.jl

using WeberElectrodynamics
using Random

const N    = 4
const dims = 2
sys = HamiltonianSystem(N, dims)

# Equal-mass 2+/2- params, kappas all 1.
const params = Float64[1, 1, 1, 1,        # m1..m4
                       1, 1, -1, -1,       # q1..q4
                       10.0,               # c
                       1, 1, 1, 1, 1, 1]   # kappas

H(q,p) = sys.hamiltonian_compiled(q, p, params)

"""Apply permutation σ of particles plus optional sign flip on positions
and/or momenta. Layout is [x1,y1,x2,y2,x3,y3,x4,y4]."""
function transform(q, p, σ; flip_q=false, flip_p=false)
    qn = similar(q); pn = similar(p)
    for i in 1:4
        j = σ[i]
        qn[2(i-1)+1] = (flip_q ? -1 : 1) * q[2(j-1)+1]
        qn[2(i-1)+2] = (flip_q ? -1 : 1) * q[2(j-1)+2]
        pn[2(i-1)+1] = (flip_p ? -1 : 1) * p[2(j-1)+1]
        pn[2(i-1)+2] = (flip_p ? -1 : 1) * p[2(j-1)+2]
    end
    return qn, pn
end

function check(name, σ; flip_q=false, flip_p=false, ntests=20)
    Random.seed!(0xC0FFEE)
    maxerr = 0.0
    for _ in 1:ntests
        q = randn(2N); p = randn(2N)
        # Spread positions so r_ij isn't tiny.
        q .*= 3.0
        H0 = H(q,p)
        qn, pn = transform(q, p, σ; flip_q=flip_q, flip_p=flip_p)
        H1 = H(qn,pn)
        err = abs(H1 - H0) / max(1.0, abs(H0))
        maxerr = max(maxerr, err)
    end
    ok = maxerr < 1e-10
    println(rpad(name, 34), "  max rel err = ", maxerr, "   ", ok ? "OK" : "FAIL")
    return ok
end

println("=== 2+/2- 4-body 2D Weber Hamiltonian: discrete symmetry checks ===")

check("S12  (1<->2)",                    [2,1,3,4])
check("S34  (3<->4)",                    [1,2,4,3])
check("S12*S34",                          [2,1,4,3])
check("P    (parity: q,p -> -q,-p)",      [1,2,3,4]; flip_q=true, flip_p=true)
check("T    (time rev: p -> -p)",         [1,2,3,4]; flip_p=true)
check("PT",                               [1,2,3,4]; flip_q=true)        # q->-q only
check("C    (1<->3, 2<->4)",              [3,4,1,2])
check("CP   (1<->3,2<->4 + parity)",      [3,4,1,2]; flip_q=true, flip_p=true)
check("CT   (1<->3,2<->4 + time rev)",    [3,4,1,2]; flip_p=true)
check("CPT",                              [3,4,1,2]; flip_q=true)
# Cross swap S13 alone should NOT be a symmetry (charges differ); sanity check.
check("S13  (1<->3) [should FAIL]",       [3,2,1,4])
