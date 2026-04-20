"""
Agent 10: Extended 4-body Weber bound orbit search.
Covers 3+/1-, unequal masses, unequal charges, c-variation, double-orbiter fine-tuning.
"""

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using LinearAlgebra

include(joinpath(@__DIR__, "ic_generators_extended.jl"))

# --------------------------------------------------------------------------
# System cache and runner
# --------------------------------------------------------------------------
const SYSTEMS = Dict{Tuple{Int,Int},Any}()
function get_system(n, d)
    key = (n, d)
    haskey(SYSTEMS, key) && return SYSTEMS[key]
    sys = HamiltonianSystem(n, d)
    SYSTEMS[key] = sys
    return sys
end

function run_4body(q0, p0, masses, charges; tmax=50.0, dt=1e-3, c=1.0, dims=2, bounce_r=0.0)
    local sys = get_system(4, dims)
    local prob = HamiltonianProblem(sys, (0.0, tmax), q0, p0;
        masses=masses, charges=charges, c=c, dt=dt)
    local sol = bounce_r > 0 ?
        solve(prob, SymmetricProjectionIntegrator(); callbacks=CollisionBounce(bounce_r)) :
        solve(prob, SymmetricProjectionIntegrator())
    return sol
end

# --------------------------------------------------------------------------
# Boundedness check
# --------------------------------------------------------------------------
function check_bound(sol; n=4, dims=2)
    local retcode = sol.retcode
    local t_final = sol.t[end]

    local en = compute_energy_timeseries(sol)
    local drift_pct = en.statistics.global_error_percent_max

    local q0 = sol.q[1]
    local d0_max = 0.0
    for i in 1:n, j in (i+1):n
        local dx = [q0[(i-1)*dims+k] - q0[(j-1)*dims+k] for k in 1:dims]
        d0_max = max(d0_max, norm(dx))
    end

    local d_max_ratio = 1.0
    local stride = max(1, length(sol.t) ÷ 200)
    for idx in 1:stride:length(sol.t)
        local q = sol.q[idx]
        for i in 1:n, j in (i+1):n
            local dx = [q[(i-1)*dims+k] - q[(j-1)*dims+k] for k in 1:dims]
            local ratio = norm(dx) / max(d0_max, eps())
            d_max_ratio = max(d_max_ratio, ratio)
        end
    end

    local is_bound = (retcode == :Success) && (d_max_ratio < 10.0) && (drift_pct < 1.0)
    # Handle NaN drift
    if isnan(drift_pct)
        is_bound = false
    end
    return (retcode=retcode, t_final=t_final, drift_pct=drift_pct,
            d_max_ratio=d_max_ratio, is_bound=is_bound)
end

# --------------------------------------------------------------------------
# CSV writer (incremental)
# --------------------------------------------------------------------------
const CSV_PATH = joinpath(@__DIR__, "survey_results.csv")
const RESULTS = Vector{NamedTuple}()

function init_csv()
    open(CSV_PATH, "w") do io
        println(io, "category,label,kwargs,retcode,t_final,drift_pct,d_max_ratio,is_bound,dims,c")
    end
end

function run_and_record(q0, p0, masses, charges; category, label, kwargs_str="",
                        tmax=100.0, dt=1e-3, c=1.0, dims=2, bounce_r=0.0)
    local sol = run_4body(q0, p0, masses, charges; tmax=tmax, dt=dt, c=c, dims=dims, bounce_r=bounce_r)
    local stats = check_bound(sol; n=4, dims=dims)
    local row = (category=category, label=label, kwargs_str=kwargs_str,
           retcode=string(stats.retcode), t_final=stats.t_final, drift_pct=stats.drift_pct,
           d_max_ratio=stats.d_max_ratio, is_bound=stats.is_bound, dims=dims, c=c)
    push!(RESULTS, row)

    # Append to CSV
    open(CSV_PATH, "a") do io
        println(io, "$(row.category),$(row.label),\"$(row.kwargs_str)\",$(row.retcode),$(round(row.t_final, digits=4)),$(round(row.drift_pct, sigdigits=6)),$(round(row.d_max_ratio, digits=4)),$(row.is_bound),$(row.dims),$(row.c)")
    end

    local bound_str = stats.is_bound ? "BOUND" : "unbound"
    println("  [$bound_str] $category / $label: t=$(round(stats.t_final, digits=2)), drift=$(round(stats.drift_pct, sigdigits=3))%, d_ratio=$(round(stats.d_max_ratio, digits=2))")
    return stats
end

init_csv()

# ======================================================================
# SECTION 1: 3+/1- configurations (~22 runs)
# ======================================================================
println("\n=== SECTION 1: 3+/1- configurations ===\n")

for eta in [0.1, 0.25, 0.5, 0.75]
    for mode in [:rotating, :breathing]
        local q0, p0, m, ch = triangular_trap_3p1m(R=2.0, energy_fraction=eta, velocity_mode=mode, dims=2)
        run_and_record(q0, p0, m, ch; category="3p1m", label="tri_R2_$(mode)_eta$(eta)",
                       kwargs_str="R=2,eta=$eta,mode=$mode", tmax=100.0, dims=2)
    end
end

for R in [1.5, 3.0, 4.0]
    local q0, p0, m, ch = triangular_trap_3p1m(R=R, energy_fraction=0.25, velocity_mode=:rotating, dims=2)
    run_and_record(q0, p0, m, ch; category="3p1m", label="tri_R$(R)_rot_eta0.25",
                   kwargs_str="R=$R,eta=0.25", tmax=100.0, dims=2)
end

# Triangular trap high-eta, larger R (promising from first run)
for R in [2.0, 3.0, 4.0]
    local q0, p0, m, ch = triangular_trap_3p1m(R=R, energy_fraction=0.75, velocity_mode=:rotating, dims=2)
    run_and_record(q0, p0, m, ch; category="3p1m", label="tri_R$(R)_rot_eta0.75",
                   kwargs_str="R=$R,eta=0.75", tmax=200.0, dims=2)
end

for kick in [0.05, 0.1, 0.2, 0.4]
    local q0, p0, m, ch = linear_3p1m(spacing=1.5, transverse_kick=kick, dims=2)
    run_and_record(q0, p0, m, ch; category="3p1m", label="linear_kick$(kick)",
                   kwargs_str="spacing=1.5,kick=$kick", tmax=100.0, dims=2)
end

# Tetrahedral 3D
for eta in [0.25, 0.5, 0.75]
    local q0, p0, m, ch = tetrahedral_3p1m(R=2.0, energy_fraction=eta, velocity_mode=:rotating)
    run_and_record(q0, p0, m, ch; category="3p1m", label="tetra_R2_rot_eta$(eta)",
                   kwargs_str="R=2,eta=$eta", tmax=100.0, dims=3)
end

# ======================================================================
# SECTION 2: Unequal masses (~16 runs)
# ======================================================================
println("\n=== SECTION 2: Unequal masses ===\n")

for eta in [0.25, 0.5, 0.75]
    local masses_hw = [10.0, 10.0, 1.0, 1.0]
    local charges_std = [1.0, 1.0, -1.0, -1.0]
    local q0, p0, m, ch = alternating_square_unequal(side=2.0, masses=masses_hw, charges=charges_std,
                                                energy_fraction=eta, velocity_mode=:rotating)
    run_and_record(q0, p0, m, ch; category="unequal_mass", label="heavy10_sq_rot_eta$(eta)",
                   kwargs_str="m=[10,10,1,1],eta=$eta", tmax=100.0, dims=2)
end

for eta in [0.25, 0.5]
    local masses_vhw = [100.0, 100.0, 1.0, 1.0]
    local charges_std = [1.0, 1.0, -1.0, -1.0]
    local q0, p0, m, ch = alternating_square_unequal(side=2.0, masses=masses_vhw, charges=charges_std,
                                                energy_fraction=eta, velocity_mode=:rotating)
    run_and_record(q0, p0, m, ch; category="unequal_mass", label="heavy100_sq_rot_eta$(eta)",
                   kwargs_str="m=[100,100,1,1],eta=$eta", tmax=100.0, dims=2)
end

for eta in [0.25, 0.5, 0.75]
    local masses_alt = [10.0, 1.0, 10.0, 1.0]
    local charges_std = [1.0, 1.0, -1.0, -1.0]
    local q0, p0, m, ch = alternating_square_unequal(side=2.0, masses=masses_alt, charges=charges_std,
                                                energy_fraction=eta, velocity_mode=:rotating)
    run_and_record(q0, p0, m, ch; category="unequal_mass", label="alt10_sq_rot_eta$(eta)",
                   kwargs_str="m=[10,1,10,1],eta=$eta", tmax=100.0, dims=2)
end

# Heavy backbone double-orbiter
for m_heavy in [10.0, 100.0]
    local masses_do = [m_heavy, m_heavy, 1.0, 1.0]
    local charges_do = [1.0, 1.0, -1.0, -1.0]
    local q0, p0, m, ch = symmetric_double_orbiter(r_pp=4.0, R=4.0, orb=1.3,
                                              masses=masses_do, charges=charges_do, dims=2)
    run_and_record(q0, p0, m, ch; category="unequal_mass",
                   label="doublorb_mh$(Int(m_heavy))",
                   kwargs_str="m=[$m_heavy,$m_heavy,1,1],rpp=4,R=4,orb=1.3",
                   tmax=200.0, dims=2)
end

# Heavy backbone rhombus
for m_heavy in [10.0, 100.0]
    local masses_rh = [m_heavy, m_heavy, 1.0, 1.0]
    local charges_rh = [1.0, 1.0, -1.0, -1.0]
    local q0, p0, m, ch = rhombus_extended(a=1.5, b=1.45, masses=masses_rh, charges=charges_rh,
                                      energy_fraction=0.75)
    run_and_record(q0, p0, m, ch; category="unequal_mass",
                   label="rhombus_mh$(Int(m_heavy))",
                   kwargs_str="m=[$m_heavy,$m_heavy,1,1],a=1.5,b=1.45,eta=0.75",
                   tmax=100.0, dims=2)
end

# ======================================================================
# SECTION 3: Unequal charges (~12 runs)
# ======================================================================
println("\n=== SECTION 3: Unequal charges ===\n")

for (qlabel, ch_vec) in [("q2p", [2.0, 2.0, -1.0, -1.0]), ("q2n", [1.0, 1.0, -2.0, -2.0])]
    for eta in [0.25, 0.5, 0.75]
        local q0, p0, m, ch = alternating_square_unequal(side=2.0, charges=ch_vec, energy_fraction=eta)
        run_and_record(q0, p0, m, ch; category="unequal_charge", label="$(qlabel)_sq_eta$(eta)",
                       kwargs_str="q=$ch_vec,eta=$eta", tmax=100.0, dims=2)
    end
end

for (qlabel, ch_vec) in [("q2p", [2.0, 2.0, -1.0, -1.0]), ("q2n", [1.0, 1.0, -2.0, -2.0])]
    local q0, p0, m, _ = symmetric_double_orbiter(r_pp=4.0, R=4.0, orb=1.3, charges=ch_vec, dims=2)
    run_and_record(q0, p0, m, ch_vec; category="unequal_charge",
                   label="doublorb_$(qlabel)_rpp4",
                   kwargs_str="q=$ch_vec,rpp=4,R=4,orb=1.3", tmax=200.0, dims=2)
end

for (qlabel, ch_vec) in [("q2p", [2.0, 2.0, -1.0, -1.0]), ("q2n", [1.0, 1.0, -2.0, -2.0])]
    local q0, p0, m, _ = symmetric_double_orbiter(r_pp=6.0, R=6.0, orb=1.3, charges=ch_vec, dims=2)
    run_and_record(q0, p0, m, ch_vec; category="unequal_charge",
                   label="doublorb_$(qlabel)_rpp6",
                   kwargs_str="q=$ch_vec,rpp=6,R=6,orb=1.3", tmax=200.0, dims=2)
end

# ======================================================================
# SECTION 4: c-variation (~12 runs, shorter tmax)
# ======================================================================
println("\n=== SECTION 4: c-variation ===\n")

# Double-orbiter 2D at c=1,2,4,10
for c_val in [1.0, 2.0, 4.0, 10.0]
    local q0, p0, m, ch = symmetric_double_orbiter(r_pp=4.0, R=4.0, orb=1.3, dims=2)
    run_and_record(q0, p0, m, ch; category="c_variation",
                   label="doublorb_2D_c$(c_val)",
                   kwargs_str="rpp=4,R=4,orb=1.3,c=$c_val",
                   tmax=600.0, dt=1e-3, c=c_val, dims=2)
end

# Double-orbiter 3D z_kick=0.13 at c=1,2,4
for c_val in [1.0, 2.0, 4.0]
    local q0, p0, m, ch = symmetric_double_orbiter(r_pp=4.0, R=4.0, orb=1.3, z_kick=0.13, dims=3)
    run_and_record(q0, p0, m, ch; category="c_variation",
                   label="doublorb_3D_c$(c_val)",
                   kwargs_str="rpp=4,R=4,orb=1.3,zk=0.13,c=$c_val",
                   tmax=600.0, dt=1e-3, c=c_val, dims=3)
end

# Rhombus at c=1,2,4
for c_val in [1.0, 2.0, 4.0]
    local q0, p0, m, ch = rhombus_extended(a=1.5, b=1.45, energy_fraction=0.75)
    run_and_record(q0, p0, m, ch; category="c_variation",
                   label="rhombus_c$(c_val)",
                   kwargs_str="a=1.5,b=1.45,eta=0.75,c=$c_val",
                   tmax=600.0, dt=1e-3, c=c_val, dims=2)
end

# ======================================================================
# SECTION 5: Double-orbiter fine-tuning (3D, ~18 runs)
# ======================================================================
println("\n=== SECTION 5: Double-orbiter fine-tuning (3D) ===\n")

for r_pp in [3.5, 3.8, 4.0, 4.2, 4.5]
    for R_val in [3.5, 4.0, 4.5]
        # Only 1 z_kick and 1 orb to keep budget
        # Skip redundant combos
        if r_pp == 4.0 && R_val == 4.0
            # orb scan
            for orb_val in [1.1, 1.2, 1.3, 1.4, 1.5]
                local q0, p0, m, ch = symmetric_double_orbiter(r_pp=r_pp, R=R_val, orb=orb_val,
                                                                z_kick=0.13, dims=3)
                run_and_record(q0, p0, m, ch; category="finetune",
                               label="rpp$(r_pp)_R$(R_val)_orb$(orb_val)_z0.13",
                               kwargs_str="rpp=$r_pp,R=$R_val,orb=$orb_val,zk=0.13",
                               tmax=600.0, dt=1e-3, c=1.0, dims=3)
            end
        elseif R_val == 4.0
            # r_pp scan at orb=1.3
            local q0, p0, m, ch = symmetric_double_orbiter(r_pp=r_pp, R=R_val, orb=1.3,
                                                            z_kick=0.13, dims=3)
            run_and_record(q0, p0, m, ch; category="finetune",
                           label="rpp$(r_pp)_R$(R_val)_orb1.3_z0.13",
                           kwargs_str="rpp=$r_pp,R=$R_val,orb=1.3,zk=0.13",
                           tmax=600.0, dt=1e-3, c=1.0, dims=3)
        elseif r_pp == 4.0
            # R scan at orb=1.3
            local q0, p0, m, ch = symmetric_double_orbiter(r_pp=r_pp, R=R_val, orb=1.3,
                                                            z_kick=0.13, dims=3)
            run_and_record(q0, p0, m, ch; category="finetune",
                           label="rpp$(r_pp)_R$(R_val)_orb1.3_z0.13",
                           kwargs_str="rpp=$r_pp,R=$R_val,orb=1.3,zk=0.13",
                           tmax=600.0, dt=1e-3, c=1.0, dims=3)
        end
    end
end

# z_kick scan at optimal (4,4,1.3)
for zk in [0.10, 0.16]
    local q0, p0, m, ch = symmetric_double_orbiter(r_pp=4.0, R=4.0, orb=1.3, z_kick=zk, dims=3)
    run_and_record(q0, p0, m, ch; category="finetune",
                   label="rpp4_R4_orb1.3_z$(zk)",
                   kwargs_str="rpp=4,R=4,orb=1.3,zk=$zk",
                   tmax=600.0, dt=1e-3, c=1.0, dims=3)
end

# ======================================================================
# Summary
# ======================================================================
println("\n=== SUMMARY ===\n")
bound_results = filter(r -> r.is_bound, RESULTS)
println("Total runs: $(length(RESULTS))")
println("Bound:      $(length(bound_results))")
println()

if !isempty(bound_results)
    println("--- Bound results by category ---")
    for cat in unique(r.category for r in bound_results)
        local cat_results = filter(r -> r.category == cat, bound_results)
        println("\n  $cat ($(length(cat_results)) bound):")
        for r in sort(collect(cat_results), by=x -> -x.t_final)
            println("    $(r.label): t*=$(round(r.t_final, digits=2)), drift=$(round(r.drift_pct, sigdigits=3))%, d_ratio=$(round(r.d_max_ratio, digits=2)), c=$(r.c)")
        end
    end
end

# c-variation scaling analysis
println("\n--- c-variation scaling analysis ---")
c_results = filter(r -> r.category == "c_variation", RESULTS)
for base in ["doublorb_2D", "doublorb_3D", "rhombus"]
    println("\n  $base:")
    local matching = sort(collect(filter(r -> startswith(r.label, base), c_results)), by=x -> x.c)
    for r in matching
        println("    c=$(r.c): t*=$(round(r.t_final, digits=2)), drift=$(round(r.drift_pct, sigdigits=3))%, bound=$(r.is_bound)")
    end
    if length(matching) >= 3
        local c_vals = [r.c for r in matching]
        local t_vals = [r.t_final for r in matching]
        local lc = log.(c_vals)
        local lt = log.(t_vals)
        local n = length(lc)
        if any(x -> x != 0.0, lc)  # avoid log(0)
            local valid = lc .!= -Inf
            local lc2 = lc[valid]
            local lt2 = lt[valid]
            if length(lc2) >= 2
                local n2 = length(lc2)
                local b = (n2 * sum(lc2 .* lt2) - sum(lc2) * sum(lt2)) / (n2 * sum(lc2.^2) - sum(lc2)^2)
                println("    Power-law fit: t* ~ c^$(round(b, digits=2))")
            end
        end
    end
end

# Fine-tuning summary
println("\n--- Fine-tuning results ---")
ft_results = sort(collect(filter(r -> r.category == "finetune", RESULTS)), by=x -> -x.t_final)
for r in ft_results
    println("  $(r.label): t*=$(round(r.t_final, digits=2)), drift=$(round(r.drift_pct, sigdigits=3))%, bound=$(r.is_bound)")
end

println("\nTotal runs written to CSV: $(length(RESULTS))")
println("Done.")
