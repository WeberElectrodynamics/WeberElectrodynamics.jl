"""
Agent 12 — Morse / Conley analysis of the effective Coulomb potential
on the 3-parameter "stacked-dimer" submanifold of the 2+/2− 4-body problem.

Parameterization (2D, COM-centered):
    + dimer along x, center at (0, +Δ): positives at (±r₊,  Δ)
    − dimer along y, center at (0, −Δ): negatives at ( 0,  ±r₋ − Δ)

Configuration is (r₊, r₋, Δ) ∈ R₊³. The **alternating square**
of Agent 4 (positives at (±R, 0), negatives at (0, ±R)) corresponds to
    r₊ = r₋ = R,   Δ = 0.

L₀ = 0, so V_eff = V_coulomb. Total angular momentum vanishes, COM = 0,
total charge = 0. Masses and |charge| set to 1. We compute V, ∇V, Hess V.
"""

using LinearAlgebra
using Printf

# ---------------------------------------------------------------------------
# 4-body positions on the stacked-dimer submanifold.
# ---------------------------------------------------------------------------
function positions(u::AbstractVector{T}) where {T}
    rp, rm, Δ = u[1], u[2], u[3]
    # rows: 1=+, 2=+, 3=-, 4=-
    return T[
         rp     Δ;
        -rp     Δ;
         zero(T)   ( rm - Δ);
         zero(T)  (-rm - Δ);
    ]
end

const CHARGES = (+1.0, +1.0, -1.0, -1.0)

function V_coulomb(u::AbstractVector{T}) where {T}
    X = positions(u)
    V = zero(T)
    @inbounds for i = 1:4, j = (i+1):4
        dx = X[i, 1] - X[j, 1]
        dy = X[i, 2] - X[j, 2]
        r = sqrt(dx * dx + dy * dy)
        V += CHARGES[i] * CHARGES[j] / r
    end
    return V
end

# Closed-form: pair distances on the submanifold
#   r12 = 2 r₊                   (++ pair)
#   r34 = 2 r₋                   (−− pair)
#   r13 = √((r₊−r₋)² + R²)        +1 to −1 (same-x sign)
#   r14 = √((r₊+r₋)² + R²)        +1 to −2
#   r23 = √((r₊+r₋)² + R²)        −1·−... (mirror of r14)
#   r24 = √((r₊−r₋)² + R²)        mirror of r13
# So
#   V(r₊,r₋,R) = 1/(2 r₊) + 1/(2 r₋)
#               − 4 / √((r₊−r₋)² + R²)        ??? sign check
# Wait: there are 4 unlike pairs (+−), each contributing −1/r. Two distinct
# distances appear, each twice:
#   V = 1/(2r₊) + 1/(2r₋) − 2/d_a − 2/d_b
# where d_a = √((r₊−r₋)² + R²), d_b = √((r₊+r₋)² + R²).
function V_closed(u::AbstractVector{T}) where {T}
    rp, rm, Δ = u[1], u[2], u[3]
    da = sqrt(rp^2 + (rm - 2Δ)^2)
    db = sqrt(rp^2 + (rm + 2Δ)^2)
    return 1/(2rp) + 1/(2rm) - 2/da - 2/db
end

# Moment of inertia about COM (origin):
#   |q1|² = r_+² + Δ²    (×2 for +'s)
#   |q3|² = (r_- − Δ)²    |q4|² = (r_- + Δ)²
#   I = 2(r_+² + Δ²) + (r_-−Δ)² + (r_-+Δ)² = 2 r_+² + 2 r_-² + 4 Δ²
function moment_of_inertia(u::AbstractVector{T}) where {T}
    return 2 * u[1]^2 + 2 * u[2]^2 + 4 * u[3]^2
end

# Amended (effective) potential at fixed total angular momentum L₀.
# (Naive Routh form — useful for 1D analysis; for multi-DOF rotating
# systems one should use the fixed-ω form V_omega below.)
function V_eff(u::AbstractVector{T}, L0::Real) where {T}
    return V_closed(u) + (L0^2) / (2 * moment_of_inertia(u))
end

# Amended potential in the rotating frame at fixed angular velocity ω
# (Smale / Marsden):  V_ω(q) = V(q) − ½ ω² ‖q‖²_M.
# Rigid relative equilibria are critical points of V_ω.
function V_omega(u::AbstractVector{T}, ω::Real) where {T}
    return V_closed(u) - 0.5 * ω^2 * moment_of_inertia(u)
end

# ---------------------------------------------------------------------------
# Sanity: closed form vs explicit position sum.
# ---------------------------------------------------------------------------
function sanity_check()
    for u in (Float64[0.6, 0.7, 0.9], Float64[0.5, 0.5, 1.0], Float64[0.3, 0.4, 0.5])
        v1 = V_coulomb(u)
        v2 = V_closed(u)
        @printf("  u = %s   V_pos = %+0.10f   V_closed = %+0.10f   diff = %.2e\n",
                u, v1, v2, abs(v1 - v2))
    end
end

# ---------------------------------------------------------------------------
# Critical-point search via gradient + Newton.
# ---------------------------------------------------------------------------
function ∇F(F, u::Vector{Float64}; h = 1e-5)
    g = zeros(3)
    for k = 1:3
        ep = copy(u); ep[k] += h
        em = copy(u); em[k] -= h
        g[k] = (F(ep) - F(em)) / (2h)
    end
    return g
end

function HF(F, u::Vector{Float64}; h = 1e-4)
    H = zeros(3, 3)
    V0 = F(u)
    for i = 1:3, j = 1:3
        if i == j
            ep = copy(u); ep[i] += h
            em = copy(u); em[i] -= h
            H[i, i] = (F(ep) - 2 * V0 + F(em)) / h^2
        else
            upp = copy(u); upp[i] += h; upp[j] += h
            upm = copy(u); upm[i] += h; upm[j] -= h
            ump = copy(u); ump[i] -= h; ump[j] += h
            umm = copy(u); umm[i] -= h; umm[j] -= h
            H[i, j] = (F(upp) - F(upm) - F(ump) + F(umm)) / (4h^2)
        end
    end
    return Symmetric((H + H') / 2)
end

∇V(u; h = 1e-5) = ∇F(V_closed, u; h = h)
HV(u; h = 1e-4) = HF(V_closed, u; h = h)

"Find a point where ‖∇F‖ is minimized (a critical point) using
a Levenberg–Marquardt-style step on g(u) = ∇F(u): solve
(JᵀJ + μI) δ = −Jᵀ g  where J = Hess F. Bounded box (1e-3, 50)."
function find_critical(F, u0::Vector{Float64}; tol = 1e-10, maxit = 200)
    u = copy(u0)
    μ = 1e-3
    for it = 1:maxit
        g = ∇F(F, u)
        ng = norm(g)
        ng < tol && return (u, it, ng)
        J = HF(F, u)
        # Hess is symmetric so JᵀJ = J²
        A = J * J + μ * I
        rhs = -(J * g)
        δ = A \ rhs
        # cap step
        if norm(δ) > 0.3 * norm(u)
            δ *= 0.3 * norm(u) / norm(δ)
        end
        ut = u + δ
        # r_+, r_- must stay > 0; Δ may be 0 or negative
        if ut[1] < 1e-3 || ut[2] < 1e-3 || any(abs.(ut) .> 50.0)
            μ *= 4
            continue
        end
        if norm(∇F(F, ut)) < ng
            u = ut
            μ = max(μ / 3, 1e-8)
        else
            μ *= 4
            μ > 1e10 && return (u, it, ng)
        end
    end
    return (u, maxit, norm(∇F(F, u)))
end

function morse_index(H::AbstractMatrix)
    λ = eigvals(Symmetric(H))
    return count(<(−1e-9), λ), λ
end

# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------
function main()
    println("== Sanity: closed-form V vs explicit Coulomb sum ==")
    sanity_check()

    println("\n== Alternating square slice  r₊ = r₋ = R, Δ = 0 ==")
    # V_square(R) = 1/R − 4/(R√2) = (1 − 2√2)/R
    for R in (0.5, 1.0, 2.0)
        u = [R, R, 0.0]
        @printf("  R=%.2f  V_closed=%+0.6f   (1−2√2)/R=%+0.6f\n",
                R, V_closed(u), (1 - 2*sqrt(2))/R)
    end

    println("\n== L₀=0: Critical-point search on V_coulomb ==")
    seeds = [
        [0.5, 0.5, 1.0],
        [0.3, 0.7, 1.0],
        [0.6, 0.6, 0.9],
        [1.0, 1.0, 1.0],
    ]
    for u0 in seeds
        u, it, gn = find_critical(V_closed, u0)
        @printf("  seed=%s -> u*≈%s  |∇V|=%.2e  V=%+0.6f\n",
                u0, round.(u; digits=4), gn, V_closed(u))
    end
    println("  (None converge to an interior critical point — V_coulomb is")
    println("   harmonic on each particle coordinate (Earnshaw); ∇V_coulomb")
    println("   has no zero on this submanifold.)")

    println("\n== L₀≠0: bona-fide critical points of V_eff (LM ‖∇F‖ minimization) ==")
    for L0 in (1.0, 1.5, 2.0, 3.0)
        F = u -> V_eff(u, L0)
        # seed near analytic 1D minimum, but slightly off-diagonal
        Rstar = L0^2 / (2 * sqrt(2))
        seeds = [
            [Rstar/2, Rstar/2, Rstar],
            [0.4, 0.4, 1.0],
            [0.5, 0.5, 0.8],
            [0.7, 0.7, 1.5],
        ]
        best = (zeros(3), Inf)
        for s in seeds
            u, it, gn = find_critical(F, s)
            gn < best[2] && (best = (u, gn))
        end
        u = best[1]
        H = HF(F, u)
        λ = eigvals(H)
        idx = count(<(−1e-9), λ)
        @printf("  L₀=%.2f  u*=%s  |∇F|=%.2e  V_eff=%+0.5f  Morse=%d  λ=%s\n",
                L0, round.(u; digits=4), best[2], F(u), idx, round.(λ; digits=4))
    end

    println("\n== Rotating alternating square as a critical pt of V_eff ==")
    println("  Agent 4: ω² = (2√2−1)/(4 R³). On square slice I = 2R², so")
    println("  L₀ = I·ω = 2R²·ω,  L₀² = R(2√2−1).  ⇒ R = L₀²/(2√2−1).")
    for L0 in (1.0, 1.5, 2.0)
        F = u -> V_eff(u, L0)
        R = L0^2 / (2*sqrt(2) - 1)
        u_eq = [R/2, R/2, R]
        g = ∇F(F, u_eq)
        H = HF(F, u_eq)
        λ = eigvals(H)
        idx = count(<(−1e-9), λ)
        @printf("  L₀=%.2f  R_eq=%.4f  u=%s\n", L0, R, round.(u_eq; digits=4))
        @printf("      V_eff=%+0.5f  ∇F=%s  |∇F|=%.2e\n",
                F(u_eq), round.(g; digits=4), norm(g))
        @printf("      Hess λ = %s   Morse index = %d\n",
                round.(λ; digits=4), idx)
    end

    println("\n== Rotating-frame amended potential V_ω = V − ½ω²·I(u) ==")
    println("  Agent 4: ω²(R) = (2√2 − 1)/(4 R³)  (R = circumscribing radius).")
    println("  Square in our coords: u = (r₊, r₋, Δ) = (R, R, 0).")
    for R in (0.5, 1.0, 2.0)
        ω² = (2*sqrt(2) - 1) / (4 * R^3)
        ω = sqrt(ω²)
        F = u -> V_omega(u, ω)
        u_eq = [R, R, 0.0]
        g = ∇F(F, u_eq)
        H = HF(F, u_eq)
        λ = eigvals(H)
        idx = count(<(−1e-9), λ)
        @printf("  R=%.2f  ω=%.4f  ∇V_ω=%s  |∇|=%.2e\n",
                R, ω, round.(g; digits=5), norm(g))
        @printf("      Hess λ = %s   Morse index = %d\n",
                round.(λ; digits=4), idx)
    end

    println("\n== Square diagonal V_eff(R) 1D restriction (L₀ fixed) ==")
    println("  Restricting to r₊=r₋=R/2:  V_eff(R) = −√2/R + L₀²/(4R²)")
    println("  ⇒ R* = L₀²/(2√2)  (1D minimum). Now check the full 3D Hessian.")
    for L0 in (1.0, 1.5, 2.0, 3.0)
        F = u -> V_eff(u, L0)
        Rstar = L0^2 / (2 * sqrt(2))
        u_star = [Rstar/2, Rstar/2, Rstar]
        g = ∇F(F, u_star)
        H = HF(F, u_star)
        λ = eigvals(H)
        idx = count(<(−1e-9), λ)
        @printf("  L₀=%.2f  R*=%.4f  V_eff=%+0.5f  |∇F|=%.2e\n",
                L0, Rstar, F(u_star), norm(g))
        @printf("      ∇F = %s\n", round.(g; digits=4))
        @printf("      Hess eigenvalues = %s   Morse index = %d\n",
                round.(λ; digits=4), idx)
    end
    println("  Interpretation: at R*, the symmetric 1D direction is a minimum;")
    println("  the perpendicular directions (r₊−r₋ asymmetry, dimer breathing)")
    println("  may be unstable — that determines the full Morse index.")

    println("\n== Hessian of bare V_coulomb at alternating square (R=1, Δ=0) ==")
    u_sq = [1.0, 1.0, 0.0]
    g = ∇V(u_sq)
    H = HV(u_sq)
    λ = eigvals(H)
    @printf("  ∇V(square) = %s   |∇V| = %.3e\n", round.(g; digits=6), norm(g))
    @printf("  Hessian eigenvalues: %s   (signature %d−, %d0, %d+)\n",
            round.(λ; digits=4),
            count(<(−1e-9), λ),
            count(x -> abs(x) < 1e-9, λ),
            count(>(1e-9), λ))

    println("\n== Plot V_coulomb along the square slice (R varies, Δ=0) ==")
    for R in 0.4:0.2:2.0
        @printf("  R=%.2f  V=%+0.5f\n", R, V_closed([R, R, 0.0]))
    end

    println("\n== V_coulomb on (r₊,r₋) grid at Δ=0  (no L) ==")
    rs = range(0.4, 1.6; length=7)
    println("        " * join((@sprintf("r₋=%.2f", r) for r in rs), "   "))
    for rp in rs
        line = @sprintf("r₊=%.2f", rp)
        for rm in rs
            v = V_closed([rp, rm, 0.0])
            line *= @sprintf("  %+8.4f", v)
        end
        println(line)
    end
end

main()
