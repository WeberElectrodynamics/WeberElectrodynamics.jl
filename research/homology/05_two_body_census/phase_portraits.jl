"""
Effective potential analysis for the 2-body Weber problem.

Computes V_eff(r) for unlike and like charge cases, identifying turning points,
critical radii, and bound-orbit energy windows.

Output: phase_portrait_data.csv with columns:
  case, q1, q2, m1, m2, c, L, r, V_eff
"""

# ---------------------------------------------------------------------------
# Effective potential for the reduced 2-body Weber problem
#
# Unlike charges (q1*q2 < 0):
#   V_eff(r) = L^2 / (2 mu r^2) + q1*q2 / r
#   Standard Coulomb-like; bound orbits for E < 0.
#
# Like charges (q1*q2 > 0):
#   On the Weber plane the Hamiltonian is
#     H = r p_r^2 / (2(r - rho)) + L^2 / (2 r^2) + q1*q2 / r
#   where rho = q1*q2 / (mu c^2).  The "effective potential" (at p_r = 0) is:
#     V_eff(r) = L^2 / (2 r^2) + q1*q2 / r
#   same functional form, but the dynamics differ because of the metric factor.
# ---------------------------------------------------------------------------

using Printf

function effective_potential(r, mu, k, L)
    # k = q1*q2 (signed)
    return L^2 / (2 * mu * r^2) + k / r
end

function critical_radius(q1, q2, m1, m2, c)
    mu = m1 * m2 / (m1 + m2)
    k = q1 * q2
    if k <= 0
        return -1.0  # no critical radius for unlike charges
    end
    return k / (mu * c^2)
end

function critical_energy(q1, q2, m1, m2, c, L)
    rho = critical_radius(q1, q2, m1, m2, c)
    if rho <= 0
        return Inf
    end
    mu = m1 * m2 / (m1 + m2)
    k = q1 * q2
    return L^2 / (2 * mu * rho^2) + k / rho
end

# --- Generate effective potential curves ---

function generate_phase_portrait_data(outfile::String)
    open(outfile, "w") do io
        println(io, "case,q1,q2,m1,m2,c,L,r,V_eff,rho")

        # Unlike charges: hydrogen-like
        for (m1, m2) in [(1.0, 1.0), (1.0, 2.0), (1.0, 10.0)]
            for c in [2.0, 4.0, 10.0, 100.0]
                for L in [0.0, 0.1, 0.5, 1.0, 2.0]
                    q1, q2 = 1.0, -1.0
                    k = q1 * q2
                    mu = m1 * m2 / (m1 + m2)
                    rho = critical_radius(q1, q2, m1, m2, c)
                    # Scan r from 0.05 to 20
                    for r in range(0.05, 20.0, length=200)
                        V = effective_potential(r, mu, k, L)
                        @printf(io, "unlike,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.6f,%.8f,%.6f\n",
                                q1, q2, m1, m2, c, L, r, V, rho)
                    end
                end
            end
        end

        # Like charges: exotic sub-critical
        for (m1, m2) in [(1.0, 1.0), (1.0, 2.0), (1.0, 10.0)]
            for c in [1.0, 2.0, 4.0, 10.0]
                for L in [0.0, 0.1, 0.5, 1.0]
                    q1, q2 = 1.0, 1.0
                    k = q1 * q2
                    mu = m1 * m2 / (m1 + m2)
                    rho = critical_radius(q1, q2, m1, m2, c)
                    h_c = critical_energy(q1, q2, m1, m2, c, L)
                    # Scan inside critical radius
                    for r in range(0.001 * rho, 0.99 * rho, length=100)
                        V = effective_potential(r, mu, k, L)
                        @printf(io, "like_sub,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.6f,%.8f,%.6f\n",
                                q1, q2, m1, m2, c, L, r, V, rho)
                    end
                    # Scan outside critical radius
                    for r in range(1.01 * rho, 20.0, length=100)
                        V = effective_potential(r, mu, k, L)
                        @printf(io, "like_super,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.6f,%.8f,%.6f\n",
                                q1, q2, m1, m2, c, L, r, V, rho)
                    end
                end
            end
        end
    end
    println("Phase portrait data written to $outfile")
end

# --- Summary table of critical parameters ---

function generate_critical_summary(outfile::String)
    open(outfile, "w") do io
        println(io, "q1,q2,m1,m2,c,mu,rho,E_min_subcrit,L_values_and_hc")

        for (m1, m2) in [(1.0, 1.0), (1.0, 2.0), (1.0, 10.0), (1.0, 100.0)]
            for c in [1.0, 2.0, 4.0, 10.0, 100.0]
                # Like charges
                q1, q2 = 1.0, 1.0
                mu = m1 * m2 / (m1 + m2)
                rho = critical_radius(q1, q2, m1, m2, c)
                E_min = mu * c^2  # minimum energy for sub-critical binding
                L_hc_pairs = String[]
                for L in [0.0, 0.1, 0.5, 1.0, 2.0]
                    hc = critical_energy(q1, q2, m1, m2, c, L)
                    push!(L_hc_pairs, @sprintf("L=%.2f:hc=%.4f", L, hc))
                end
                @printf(io, "%.1f,%.1f,%.1f,%.1f,%.1f,%.6f,%.6f,%.6f,%s\n",
                        q1, q2, m1, m2, c, mu, rho, E_min, join(L_hc_pairs, ";"))
            end
        end
    end
    println("Critical summary written to $outfile")
end

# Run
basedir = @__DIR__
generate_phase_portrait_data(joinpath(basedir, "phase_portrait_data.csv"))
generate_critical_summary(joinpath(basedir, "critical_summary.csv"))
println("Phase portrait analysis complete.")
