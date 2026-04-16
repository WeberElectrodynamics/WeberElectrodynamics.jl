"""
Contact-type energy grid for 2-body Weber Hamiltonian.

Scans energy E over a range and for each E reports:
  (a) Is the Hill region star-shaped?
  (b) From what center?
  (c) % monotone rays

Writes results to contact_type_grid.csv.
"""

using Random, LinearAlgebra, Printf
using WeberElectrodynamics

include(joinpath(@__DIR__, "star_center_search.jl"))

# ============================================================================
# Unlike charges: q1=+1, q2=-1, m1=m2=1, c=1
# ============================================================================

function scan_unlike(;n_energies=50, n_rays=200)
    charges = [1.0, -1.0]
    masses = [1.0, 1.0]
    c = 1.0
    μ = 0.5

    V = make_2body_V(charges, masses, c)

    results = []

    energies = range(-10.0, -0.05; length=n_energies)

    for E in energies
        # Hill region for unlike charges: V(r) = -1/r <= E (E<0)
        # => 1/r >= |E| => r <= 1/|E|. Region is compact: (0, 1/|E|].
        # Place center at r = 1/(2|E|) where V = -2|E| < E (deep inside).
        r_boundary = 1.0 / abs(E)
        r_center = r_boundary / 2.0
        q_star = [r_center/2, 0.0, -r_center/2, 0.0]

        V_star = V(q_star)

        if V_star > E
            push!(results, (E=E, charge_type="unlike", star_shaped=false,
                           center_r=NaN, frac_monotone=0.0, frac_single=0.0))
            continue
        end

        t_max_ray = max(10.0, 3 * r_boundary)
        res = shoot_rays(V, q_star, E, n_rays; t_max=t_max_ray, n_samples=300, seed=42)

        # Try optimization if not already good
        if res.frac_monotone < 0.95
            opt = optimize_star_center(V, E, q_star, min(n_rays, 80);
                                       n_iters=80, step_size=r_center*0.3,
                                       t_max=t_max_ray, n_samples=200)
            if opt.frac_monotone > res.frac_monotone
                r_opt = norm([opt.q_star[1] - opt.q_star[3],
                              opt.q_star[2] - opt.q_star[4]])
                push!(results, (E=E, charge_type="unlike",
                               star_shaped=(opt.frac_monotone > 0.95),
                               center_r=r_opt, frac_monotone=opt.frac_monotone,
                               frac_single=opt.frac_single_cross))
                continue
            end
        end

        push!(results, (E=E, charge_type="unlike",
                       star_shaped=(res.frac_monotone > 0.95),
                       center_r=r_center, frac_monotone=res.frac_monotone,
                       frac_single=res.frac_single_cross))
    end

    return results
end

# ============================================================================
# Like charges: q1=+1, q2=+1, supercritical region r > ρ = 2.0
# ============================================================================

function scan_like(;n_energies=50, n_rays=200)
    charges = [1.0, 1.0]
    masses = [1.0, 1.0]
    c = 1.0
    μ = 0.5

    V = make_2body_V(charges, masses, c)

    # Critical radius: ρ = q1*q2/(μ*c²) = 1/(0.5*1) = 2.0
    ρ = 2.0

    results = []

    # V(r) = +1/r for like charges (at p=0). In supercritical region r > ρ = 2:
    # V ranges from 1/ρ = 0.5 down to 0 as r→∞.
    # Hill region {V <= E, r > ρ} exists only for E >= V(ρ) = 1/ρ = 0.5.
    # Wait — V(ρ) = 1/ρ = 0.5, V(∞) = 0. So {V <= E} with E >= 0.5 is all r > ρ.
    # For E < 0.5: {V <= E, r > ρ} = {1/r <= E, r > ρ} = {r >= 1/E, r > ρ} = {r >= 1/E}
    # (since 1/E > ρ = 2 when E < 0.5).

    energies = range(0.05, 0.5; length=n_energies)

    for E in energies
        r_inner = max(1.0 / E, ρ + 0.01)  # inner boundary of Hill region in supercritical zone
        r_center = r_inner + 1.0
        q_star = [r_center/2, 0.0, -r_center/2, 0.0]

        V_star = V(q_star)
        if V_star > E
            r_center = r_inner + 5.0
            q_star = [r_center/2, 0.0, -r_center/2, 0.0]
            V_star = V(q_star)
        end

        if V_star > E
            push!(results, (E=E, charge_type="like_supercrit", star_shaped=false,
                           center_r=NaN, frac_monotone=0.0, frac_single=0.0))
            continue
        end

        res = shoot_rays(V, q_star, E, n_rays; t_max=30.0, n_samples=300, seed=42)

        if res.frac_monotone < 0.95
            opt = optimize_star_center(V, E, q_star, min(n_rays, 100);
                                       n_iters=100, step_size=0.5,
                                       t_max=30.0, n_samples=200)
            if opt.frac_monotone > res.frac_monotone
                r_center = norm([opt.q_star[1] - opt.q_star[3],
                                 opt.q_star[2] - opt.q_star[4]])
                push!(results, (E=E, charge_type="like_supercrit",
                               star_shaped=(opt.frac_monotone > 0.95),
                               center_r=r_center, frac_monotone=opt.frac_monotone,
                               frac_single=opt.frac_single_cross))
                continue
            end
        end

        r_sep = norm([q_star[1] - q_star[3], q_star[2] - q_star[4]])
        push!(results, (E=E, charge_type="like_supercrit",
                       star_shaped=(res.frac_monotone > 0.95),
                       center_r=r_sep, frac_monotone=res.frac_monotone,
                       frac_single=res.frac_single_cross))
    end

    return results
end

# ============================================================================
# Write CSV
# ============================================================================

function write_csv(results, filename)
    open(filename, "w") do f
        println(f, "E,charge_type,star_shaped,center_r,frac_monotone,frac_single_cross")
        for r in results
            @printf(f, "%.6f,%s,%s,%.6f,%.4f,%.4f\n",
                    r.E, r.charge_type, r.star_shaped, r.center_r,
                    r.frac_monotone, r.frac_single)
        end
    end
end

function main()
    println("Scanning unlike charges (+1,-1)...")
    res_unlike = scan_unlike(n_energies=50, n_rays=200)

    println("Scanning like charges (+1,+1) supercritical...")
    res_like = scan_like(n_energies=50, n_rays=200)

    all_results = vcat(res_unlike, res_like)
    outfile = joinpath(@__DIR__, "contact_type_grid.csv")
    write_csv(all_results, outfile)
    println("Wrote $(length(all_results)) rows to $outfile")

    # Print summary
    println("\n--- UNLIKE CHARGES SUMMARY ---")
    n_star = count(r -> r.star_shaped, res_unlike)
    @printf "Star-shaped: %d / %d (%.0f%%)\n" n_star length(res_unlike) 100*n_star/length(res_unlike)
    @printf "Mean monotone fraction: %.1f%%\n" 100*mean(r.frac_monotone for r in res_unlike)

    println("\n--- LIKE CHARGES (SUPERCRITICAL) SUMMARY ---")
    n_star_like = count(r -> r.star_shaped, res_like)
    @printf "Star-shaped: %d / %d (%.0f%%)\n" n_star_like length(res_like) 100*n_star_like/length(res_like)
    @printf "Mean monotone fraction: %.1f%%\n" 100*mean(r.frac_monotone for r in res_like)
end

# Need Statistics for mean
using Statistics

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
