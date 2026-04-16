"""
Analyze census_results.csv and print summary statistics.
"""

function analyze()
    lines = readlines(joinpath(@__DIR__, "census_results.csv"))
    header = split(lines[1], ",")
    data = [split(l, ",") for l in lines[2:end]]
    N = length(data)
    println("Total runs: $N")

    # Parse helper
    pf(s) = (v = tryparse(Float64, s); v === nothing ? NaN : v)

    # ----- Overall counts -----
    println("\n========== OVERALL ==========")
    labels = Dict{String,Int}()
    retcodes = Dict{String,Int}()
    types = Dict{String,Int}()
    bound_count = 0
    for row in data
        lab = row[16]
        labels[lab] = get(labels, lab, 0) + 1
        rc = row[15]
        retcodes[rc] = get(retcodes, rc, 0) + 1
        ot = row[11]
        types[ot] = get(types, ot, 0) + 1
        row[9] == "true" && (bound_count += 1)
    end
    println("By case label:")
    for (k, v) in sort(collect(labels))
        bnd = count(r -> r[16] == k && r[9] == "true", data)
        println("  $k: $v runs, $bnd bound")
    end
    println("\nReturn codes:")
    for (k, v) in sort(collect(retcodes))
        println("  $k: $v")
    end
    println("\nOrbit types:")
    for (k, v) in sort(collect(types))
        println("  $k: $v")
    end
    println("\nTotal bound: $bound_count / $N")

    # ----- Unlike charges detail -----
    println("\n========== UNLIKE CHARGES (q1*q2 < 0) ==========")

    unlike = filter(r -> r[16] == "unlike", data)
    println("Total unlike runs: $(length(unlike))")

    # By c
    println("\nBound by c:")
    for c_val in ["1.00", "2.00", "4.00", "10.00", "100.00"]
        sub = filter(r -> r[5] == c_val, unlike)
        bnd = count(r -> r[9] == "true", sub)
        println("  c=$c_val: $bnd/$(length(sub)) bound")
    end

    # By L
    println("\nBound by L:")
    for L_val in ["0.0000", "0.1000", "0.2500", "0.5000", "1.0000", "2.0000"]
        sub = filter(r -> r[7] == L_val, unlike)
        bnd = count(r -> r[9] == "true", sub)
        println("  L=$L_val: $bnd/$(length(sub)) bound")
    end

    # By mass ratio
    println("\nBound by mass ratio:")
    for (m1s, m2s) in [("1.00","1.00"), ("1.00","2.00"), ("1.00","10.00"), ("1.00","100.00")]
        sub = filter(r -> r[3] == m1s && r[4] == m2s, unlike)
        bnd = count(r -> r[9] == "true", sub)
        println("  m1=$m1s, m2=$m2s: $bnd/$(length(sub)) bound")
    end

    # By energy
    println("\nBound by energy:")
    E_vals = sort(unique([r[8] for r in unlike]))
    for e in E_vals
        sub = filter(r -> r[8] == e, unlike)
        bnd = count(r -> r[9] == "true", sub)
        println("  E=$e: $bnd/$(length(sub)) bound")
    end

    # Bound orbit details
    println("\nBound unlike orbits (sample):")
    bound_unlike = filter(r -> r[9] == "true", unlike)
    for (i, row) in enumerate(bound_unlike)
        i > 30 && break
        println("  m=($(row[3]),$(row[4])) c=$(row[5]) L=$(row[7]) E=$(row[8]) r0=$(row[6]) " *
                "type=$(row[11]) T=$(row[10]) drift=$(row[12])%")
    end

    # Energy drift statistics for bound orbits
    drifts_unlike = Float64[]
    for row in bound_unlike
        d = pf(row[12])
        !isnan(d) && push!(drifts_unlike, d)
    end
    if !isempty(drifts_unlike)
        sort!(drifts_unlike)
        println("\nEnergy drift (bound unlike):")
        println("  N=$(length(drifts_unlike))")
        println("  min=$(minimum(drifts_unlike))%")
        println("  median=$(drifts_unlike[length(drifts_unlike)÷2+1])%")
        println("  max=$(maximum(drifts_unlike))%")
        println("  <0.01%: $(count(d -> d < 0.01, drifts_unlike))")
        println("  <0.1%: $(count(d -> d < 0.1, drifts_unlike))")
        println("  <1.0%: $(count(d -> d < 1.0, drifts_unlike))")
    end

    # Period detection
    periods_found = count(r -> r[9] == "true" && pf(r[10]) > 0 && !isnan(pf(r[10])), unlike)
    println("\nPeriods detected: $periods_found / $(length(bound_unlike)) bound orbits")

    # ----- Like charges detail -----
    println("\n========== LIKE CHARGES (q1*q2 > 0) ==========")

    like = filter(r -> r[16] == "like", data)
    println("Total like runs: $(length(like))")

    println("\nBy retcode:")
    for rc in ["Success", "Failure", "Error"]
        cnt = count(r -> r[15] == rc, like)
        cnt > 0 && println("  $rc: $cnt")
    end

    println("\nBound by c:")
    for c_val in ["1.00", "2.00", "4.00", "10.00", "100.00"]
        sub = filter(r -> r[5] == c_val, like)
        bnd = count(r -> r[9] == "true", sub)
        succ = count(r -> r[15] == "Success", sub)
        println("  c=$c_val: $bnd bound, $succ success / $(length(sub)) total")
    end

    println("\nBound by L:")
    for L_val in ["0.0000", "0.1000", "0.5000"]
        sub = filter(r -> r[7] == L_val, like)
        bnd = count(r -> r[9] == "true", sub)
        succ = count(r -> r[15] == "Success", sub)
        println("  L=$L_val: $bnd bound, $succ success / $(length(sub)) total")
    end

    # Success cases detail
    println("\nSuccessful like-charge runs:")
    succ_like = filter(r -> r[15] == "Success", like)
    for row in succ_like
        println("  m=($(row[3]),$(row[4])) c=$(row[5]) L=$(row[7]) E=$(row[8]) r0=$(row[6]) " *
                "bound=$(row[9]) type=$(row[11]) drift=$(row[12])%")
    end

    # Bound like-charge runs
    bound_like = filter(r -> r[9] == "true", like)
    if !isempty(bound_like)
        println("\nBound like-charge orbits:")
        for row in bound_like
            println("  m=($(row[3]),$(row[4])) c=$(row[5]) L=$(row[7]) E=$(row[8]) r0=$(row[6]) " *
                    "type=$(row[11]) T=$(row[10]) drift=$(row[12])%")
        end
    else
        println("\nNo bound like-charge orbits found.")
        println("This is expected: the symplectic integrator's fixed-point projection")
        println("fails for most sub-critical configurations due to the metric singularity")
        println("at r=rho. The Weber-metric signature change makes standard symplectic")
        println("methods unreliable in this regime.")
    end

    # ----- Asymmetric charges detail -----
    println("\n========== ASYMMETRIC CHARGES ==========")

    for label in ["asym_1_-2", "asym_2_-1", "asym_1_-0.5", "asym_0.5_-1"]
        sub = filter(r -> r[16] == label, data)
        bnd = count(r -> r[9] == "true", sub)
        succ = count(r -> r[15] == "Success", sub)
        println("\n$label: $bnd bound, $succ success / $(length(sub)) total")
        bound_sub = filter(r -> r[9] == "true", sub)
        for row in bound_sub
            println("  c=$(row[5]) m=($(row[3]),$(row[4])) L=$(row[7]) E=$(row[8]) " *
                    "type=$(row[11]) T=$(row[10]) drift=$(row[12])%")
        end
    end

    # ----- Weber effects: precession analysis -----
    println("\n========== WEBER EFFECTS ==========")
    println("\nCompare unlike orbits at different c (same E, L, masses):")
    # Group by (m1, m2, L, E) and compare across c
    groups = Dict{String,Vector}()
    for row in bound_unlike
        key = "m=($(row[3]),$(row[4])) L=$(row[7]) E=$(row[8])"
        push!(get!(groups, key, []), row)
    end
    for (key, rows) in sort(collect(groups))
        length(rows) > 1 || continue
        println("  $key:")
        for row in sort(rows, by=r->pf(r[5]))
            T = pf(row[10])
            T_str = isnan(T) ? "NaN" : string(round(T, digits=4))
            println("    c=$(row[5]) type=$(row[11]) T=$T_str r0=$(row[6]) drift=$(row[12])%")
        end
    end
end

analyze()
