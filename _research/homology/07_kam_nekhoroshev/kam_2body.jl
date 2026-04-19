#= KAM twist condition verification for 2-body Weber electrodynamics.

   Computes action-angle variables, verifies the twist condition, and estimates
   the perturbation size in Delaunay coordinates for c = 4, 10, 100.

   Run:  julia --project=. research/homology/07_kam_nekhoroshev/kam_2body.jl
=#

using Printf

# ── Physical setup ───────────────────────────────────────────────────────────
# Two opposite charges: q1 = +1, q2 = -1, equal masses m1 = m2 = 1.
# Reduced mass μ = m1*m2/(m1+m2) = 0.5
# Coupling k = q1*q2 = -1  (attractive)
# |k| = 1

const μ = 0.5
const k = -1.0      # q1*q2
const abs_k = 1.0   # |k|

# ── Kepler action-angle variables ────────────────────────────────────────────
#
# For bound Kepler orbits (E < 0):
#   Semi-major axis:  a = |k|/(2|E|)
#   Total action:     I = |k|√μ / √(2|E|)  =  √(μ |k| a)
#   Angular momentum: L = I √(1 - e²)
#   Radial action:    I_r = I - L
#
# Hamiltonian in action variables:
#   H₀(I) = -μ k² / (2 I²)
#
# Frequencies:
#   ω = ∂H₀/∂I = μ k² / I³
#
# Twist:
#   ∂²H₀/∂I² = -3μ k² / I⁴

"""
    kepler_action(E, L)

Compute the Delaunay total action I from energy E and angular momentum L.
Requires E < 0 (bound orbit).
"""
function kepler_action(E::Float64)
    @assert E < 0 "Bound orbit requires E < 0"
    return abs_k * sqrt(μ) / sqrt(2 * abs(E))
end

"""
    kepler_semimajor(E)

Semi-major axis from energy.
"""
kepler_semimajor(E::Float64) = abs_k / (2 * abs(E))

"""
    kepler_eccentricity(E, L)

Eccentricity from energy and angular momentum.
"""
function kepler_eccentricity(E::Float64, L::Float64)
    a = kepler_semimajor(E)
    e2 = 1.0 - L^2 / (μ * abs_k * a)
    return sqrt(max(e2, 0.0))
end

"""
    kepler_frequency(I)

Kepler orbital frequency ω = μk²/I³.
"""
kepler_frequency(I::Float64) = μ * abs_k^2 / I^3

"""
    kepler_twist(I)

Twist coefficient ∂²H₀/∂I² = -3μk²/I⁴.
"""
kepler_twist(I::Float64) = -3.0 * μ * abs_k^2 / I^4

# ── Weber perturbation in action-angle variables ─────────────────────────────
#
# H₁ = (q₁q₂)/(2μ²) × p_r² / r  =  (-1)/(2×0.25) × p_r² / r  =  -2 p_r²/r
#
# On a Kepler ellipse parametrized by eccentric anomaly E:
#   r = a(1 - e cos E)
#   p_r = μ ṙ = μ × (|k| e sin E) / (I × (1 - e cos E)) × r
#        ... simplifying:  p_r = (μ |k| e sin E) / I
#
# So:  p_r² = μ² k² e² sin²E / I²
#
# And: H₁ = k/(2μ²) × μ² k² e² sin²E / (I² a (1 - e cos E))
#         = k³ e² sin²E / (2 I² a (1 - e cos E))
#
# Note k < 0, so k³ < 0.  And H₁ enters as ε H₁ with ε = 1/c².
# The physical perturbation to the energy is:  δH = H₁/c²
#
# Period average:  ⟨sin²E/(1-e cos E)⟩ = 1  (standard Kepler integral)
#
# So:  ⟨H₁⟩ = k³ e² / (2 I² a) = (k < 0 → k³ < 0)
#
# Since a = I²/(μ|k|):
#   ⟨H₁⟩ = k³ e² / (2 I² × I²/(μ|k|))  =  k³ e² μ |k| / (2 I⁴)
#         = -|k|⁴ μ e² / (2 I⁴)   (since k³ = -|k|³ and |k|³×|k| = |k|⁴)
#         = -μ k⁴ e² / (2 I⁴)   (since k⁴ = |k|⁴)

"""
    weber_perturbation_avg(I, e)

Period-averaged Weber perturbation ⟨H₁⟩ in Delaunay variables.
"""
function weber_perturbation_avg(I::Float64, e::Float64)
    return -μ * abs_k^4 * e^2 / (2.0 * I^4)
end

"""
    perturbation_ratio(I, e, c)

Ratio |ε⟨H₁⟩/H₀| — relative size of the Weber perturbation.
"""
function perturbation_ratio(I::Float64, e::Float64, c::Float64)
    H0 = -μ * abs_k^2 / (2.0 * I^2)
    H1_avg = weber_perturbation_avg(I, e)
    ε = 1.0 / c^2
    return abs(ε * H1_avg / H0)
end

# ── Circular orbit reference ─────────────────────────────────────────────────
# For a circular orbit at radius r₀:
#   E = k/(2r₀) = -1/(2r₀)
#   L = √(μ |k| r₀) = √(r₀/2)
#   v_circ = √(|k|/(μ r₀)) × L/r₀ ... actually v = L/(μ r₀) = √(|k|/(μ r₀))

function circular_orbit(r0::Float64)
    E = k / (2.0 * r0)       # k < 0 → E < 0
    L = sqrt(μ * abs_k * r0) # angular momentum
    v = L / (μ * r0)         # orbital speed
    I = kepler_action(E)
    e = 0.0                   # circular
    ω = kepler_frequency(I)
    T = 2π / ω               # period
    return (; r0, E, L, v, I, e, ω, T)
end

function eccentric_orbit(r0::Float64, e::Float64)
    a = r0                    # semi-major axis
    E_energy = -abs_k / (2.0 * a)
    L = sqrt(μ * abs_k * a * (1 - e^2))
    I = kepler_action(E_energy)
    ω = kepler_frequency(I)
    T = 2π / ω
    v_peri = L / (μ * a * (1 - e))   # speed at periapsis
    v_apo = L / (μ * a * (1 + e))    # speed at apoapsis
    return (; a, E=E_energy, L, I, e, ω, T, v_peri, v_apo)
end

# ── Nekhoroshev stability time ───────────────────────────────────────────────

function nekhoroshev_time(c::Float64, n_dof::Int; method=:lochak)
    ε = 1.0 / c^2
    if method == :lochak
        a_exp = 1.0 / (2 * n_dof)
    elseif method == :poschel
        a_exp = 1.0 / (2 * (n_dof - 1))
    else
        error("Unknown method: $method")
    end
    return exp((1.0 / ε)^a_exp)
end

# ── KAM measure estimate ────────────────────────────────────────────────────

"""
    kam_fraction(c)

Estimated fraction of phase space occupied by KAM tori.
Scaling: 1 - C/c for some O(1) constant C.
"""
kam_fraction(c::Float64) = max(0.0, 1.0 - 1.0 / c)

# ══════════════════════════════════════════════════════════════════════════════
#  Main computation
# ══════════════════════════════════════════════════════════════════════════════

function main()
    println("="^76)
    println("  KAM / Nekhoroshev analysis for 2-body Weber electrodynamics")
    println("="^76)
    println()

    # ── 1. Reference circular orbits at different radii ──────────────────
    println("─"^76)
    println("  1. Reference circular Kepler orbits (q₁=+1, q₂=-1, m₁=m₂=1)")
    println("─"^76)
    for r0 in [1.0, 2.0, 4.0, 8.0]
        orb = circular_orbit(r0)
        @printf("  r₀ = %4.1f:  E = %8.4f  L = %6.4f  v = %6.4f  I = %6.4f  ω = %8.5f  T = %8.4f\n",
                r0, orb.E, orb.L, orb.v, orb.I, orb.ω, orb.T)
    end
    println()

    # ── 2. Twist condition verification ──────────────────────────────────
    println("─"^76)
    println("  2. KAM twist condition: ∂²H₀/∂I² = -3μk²/I⁴")
    println("─"^76)
    for r0 in [1.0, 2.0, 4.0, 8.0]
        orb = circular_orbit(r0)
        tw = kepler_twist(orb.I)
        @printf("  r₀ = %4.1f:  I = %6.4f  ∂²H₀/∂I² = %10.6f  (nonzero ✓)\n",
                r0, orb.I, tw)
    end
    println()
    println("  RESULT: Twist condition ∂²H₀/∂I² ≠ 0 is satisfied for ALL bound orbits.")
    println("  The Kepler Hamiltonian H₀ = -μk²/(2I²) is convex in I.")
    println()

    # ── 3. Weber perturbation size at c = 4, 10, 100 ────────────────────
    println("─"^76)
    println("  3. Weber perturbation size |ε⟨H₁⟩/H₀| for various (r₀, e, c)")
    println("─"^76)
    c_values = [1.0, 4.0, 10.0, 100.0]
    e_values = [0.0, 0.1, 0.3, 0.5, 0.8]
    r0 = 4.0  # matching the double-orbiter scale

    orb_circ = circular_orbit(r0)
    @printf("  Reference orbit: r₀ = %.1f, I = %.4f, ω = %.5f, T = %.4f\n\n",
            r0, orb_circ.I, orb_circ.ω, orb_circ.T)

    # Header
    @printf("  %6s", "e \\ c")
    for c in c_values
        @printf("  %12.0f", c)
    end
    println()
    println("  " * "-"^60)

    for e in e_values
        orb = eccentric_orbit(r0, e)
        @printf("  %6.2f", e)
        for c in c_values
            ratio = perturbation_ratio(orb.I, e, c)
            @printf("  %12.6f", ratio)
        end
        println()
    end
    println()
    println("  NOTE: For e=0 (circular), the radial Weber perturbation vanishes identically.")
    println("  The Weber force only affects non-circular (ṙ≠0) motion.")
    println()

    # ── 4. Nekhoroshev stability times ───────────────────────────────────
    println("─"^76)
    println("  4. Nekhoroshev stability times T_N (units of T₀)")
    println("─"^76)
    n_dof_values = [2, 4, 6, 9]
    @printf("  %6s", "c \\ n")
    for n in n_dof_values
        @printf("  %12d", n)
    end
    println()
    println("  " * "-"^60)

    for c in [1.0, 2.0, 4.0, 10.0, 100.0, 1000.0]
        @printf("  %6.0f", c)
        for n in n_dof_values
            TN = nekhoroshev_time(c, n)
            if TN > 1e6
                @printf("  %12.2e", TN)
            else
                @printf("  %12.2f", TN)
            end
        end
        println()
    end
    println()
    println("  n=2: 2-body in 2D (after reduction)")
    println("  n=4: 4-body in 2D (rhombus, after COM + symmetry)")
    println("  n=6: 4-body in 3D (double-orbiter, after COM + Klein-four)")
    println("  n=9: 4-body in 3D (no symmetry reduction, after COM)")
    println()

    # ── 5. KAM fraction estimates ────────────────────────────────────────
    println("─"^76)
    println("  5. Estimated KAM torus fraction vs c (heuristic: 1 - 1/c)")
    println("─"^76)
    for c in [1.0, 2.0, 4.0, 10.0, 100.0]
        frac = kam_fraction(c)
        @printf("  c = %6.0f:  ε = 1/c² = %10.6f  √ε = %8.5f  KAM fraction ≈ %6.3f\n",
                c, 1.0/c^2, 1.0/c, frac)
    end
    println()

    # ── 6. Comparison with observed quasi-bound lifetimes ────────────────
    println("─"^76)
    println("  6. Comparison with observed quasi-bound candidates (c=1)")
    println("─"^76)
    println()

    candidates = [
        ("3D Double-orbiter", 0.0205, 566.0, 6),
        ("Rhombus basin",     0.18,   440.0, 4),
        ("Breathing square",  228.6,  NaN,   4),
    ]

    for (name, λ, t_star, n) in candidates
        TN = nekhoroshev_time(1.0, n)
        τ_e = 1.0 / λ
        n_efold = isnan(t_star) ? NaN : t_star * λ

        @printf("  %-25s  λ_max = %8.4f  τ_e = %8.2f  t* = %7.1f  T_N = %5.2f  t*/τ_e = %6.1f\n",
                name, λ, τ_e, t_star, TN, n_efold)
    end
    println()
    println("  All candidates: t* >> T_N (Nekhoroshev bound is trivially satisfied).")
    println("  Double-orbiter: t*/τ_e ≈ 11.6 e-foldings before escape.")
    println("  Rhombus: t*/τ_e ≈ 79 e-foldings — wider basin relative to diffusion rate.")
    println()

    # ── 7. Scaling predictions ───────────────────────────────────────────
    println("─"^76)
    println("  7. Scaling predictions: t*(c) assuming fast-diffusion (t* ∝ c²)")
    println("─"^76)
    t_star_ref = 566.0  # at c=1
    for c in [1.0, 2.0, 4.0, 10.0, 100.0]
        t_pred = t_star_ref * c^2
        @printf("  c = %6.0f:  t*(predicted) = %12.1f\n", c, t_pred)
    end
    println()
    println("  These predictions are testable via c-continuation experiments.")
    println("  If t* grows faster than c², Nekhoroshev exponential stability is active.")
    println()

    println("="^76)
    println("  SUMMARY")
    println("="^76)
    println()
    println("  1. The KAM twist condition is satisfied for all 2-body Kepler orbits.")
    println("  2. The Weber perturbation is O(e²/c²) — small for low eccentricity or large c.")
    println("  3. At c=1 (ε=1), the perturbation is O(1): KAM/Nekhoroshev bounds are trivial.")
    println("  4. All observed quasi-bound lifetimes (t*~400-566) vastly exceed Nekhoroshev")
    println("     lower bounds (T_N~3), confirming the system is in the fast Arnold diffusion")
    println("     regime, not the Nekhoroshev exponential stability regime.")
    println("  5. Testable prediction: t* ∝ c² in the fast-diffusion regime.")
    println()
end

main()
