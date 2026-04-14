#!/usr/bin/env julia
# Reproducibility harness for the 4-body 2+/2- numerical atlas.
#
# Usage:
#     julia --project=. research/FourBodyTwoPlusTwoMinus/14_numerical_atlas/reproduce.jl
#     julia --project=. .../reproduce.jl --subset F1b
#
# Rebuilds every row of atlas.csv that has non-empty q0_json/p0_json,
# re-runs the integration via SharedSurvey.run, and prints PASS/FAIL
# against the stored retcode / E drift.
#
# NOTE: this script is NOT run by CI. It is a manual reproducibility tool.

const HERE   = @__DIR__
const SHARED = normpath(joinpath(HERE, "..", "shared"))
include(joinpath(SHARED, "run_survey.jl"))
using .SharedSurvey
using Printf

# Tiny CSV parser handling quoted fields with embedded commas. Returns
# Vector{Dict{String,String}}.
function parse_csv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return Dict{String,String}[]
    header = _split_csv_line(lines[1])
    rows = Dict{String,String}[]
    for l in lines[2:end]
        isempty(strip(l)) && continue
        fields = _split_csv_line(l)
        if length(fields) != length(header)
            @warn "skipping malformed row (nfields=$(length(fields)))" line=l
            continue
        end
        push!(rows, Dict(header[i] => fields[i] for i in eachindex(header)))
    end
    return rows
end

function _split_csv_line(l::AbstractString)
    out = String[]
    buf = IOBuffer()
    inq = false
    i = 1
    chars = collect(l)
    while i <= length(chars)
        c = chars[i]
        if inq
            if c == '"'
                if i+1 <= length(chars) && chars[i+1] == '"'
                    write(buf, '"'); i += 2; continue
                else
                    inq = false; i += 1; continue
                end
            else
                write(buf, c); i += 1
            end
        else
            if c == '"'; inq = true; i += 1
            elseif c == ','; push!(out, String(take!(buf))); i += 1
            else; write(buf, c); i += 1
            end
        end
    end
    push!(out, String(take!(buf)))
    return out
end

# Parse a JSON-ish flat numeric array "[1.0,2.0,...]" without a JSON dep.
function parse_jarray(s::AbstractString)
    s = strip(s)
    (isempty(s) || s == "[]") && return Float64[]
    s = replace(s, "[" => "", "]" => "")
    return [parse(Float64, strip(t)) for t in split(s, ",") if !isempty(strip(t))]
end

function parse_jarray_any(s::AbstractString)
    # same as above for masses/charges
    parse_jarray(s)
end

function main(args)
    subset = nothing
    for a in args
        if startswith(a, "--subset")
            subset = length(a) > 9 ? a[10:end] : (isempty(args) ? nothing : popfirst!(args))
        end
    end
    # crude: also accept "--subset PAT" as two args
    for (i, a) in enumerate(args)
        if a == "--subset" && i < length(args)
            subset = args[i+1]
        end
    end

    atlas = joinpath(HERE, "atlas.csv")
    rows = parse_csv(atlas)
    @printf("Loaded %d rows from %s\n", length(rows), atlas)

    npass = 0; nfail = 0; nskip = 0
    for r in rows
        id = r["id"]
        if subset !== nothing && !occursin(subset, id)
            continue
        end
        if isempty(r["q0_json"]) || isempty(r["p0_json"])
            @printf("[SKIP] %-40s  (no q0/p0 stored)\n", id)
            nskip += 1; continue
        end
        try
            q0 = parse_jarray(r["q0_json"])
            p0 = parse_jarray(r["p0_json"])
            masses  = parse_jarray_any(r["masses"])
            charges = parse_jarray_any(r["charges"])
            c    = parse(Float64, r["c"])
            dt   = parse(Float64, r["dt"])
            tmax = parse(Float64, r["tmax"])
            za   = parse(Float64, r["zollner_a"])
            br   = parse(Float64, r["bounce_r"])
            sol  = SharedSurvey.run(q0, p0, masses, charges;
                                    tmax=tmax, dt=dt, c=c,
                                    zollner_enabled = za > 0,
                                    zollner_a = za, bounce_r = br)
            s = SharedSurvey.summarize(sol)
            stored_rc    = r["retcode"]
            stored_drift = tryparse(Float64, r["E_drift_pct"])
            rc_ok    = occursin(string(s.retcode), stored_rc) || occursin(stored_rc, string(s.retcode))
            drift_ok = stored_drift === nothing || isnan(stored_drift) ||
                       abs(s.E_drift_pct - stored_drift) < max(0.5, 0.1 * abs(stored_drift))
            status = (rc_ok && drift_ok) ? "PASS" : "FAIL"
            @printf("[%s] %-40s  rc=%s stored=%s drift=%.3g stored=%s\n",
                    status, id, s.retcode, stored_rc, s.E_drift_pct, r["E_drift_pct"])
            status == "PASS" ? (npass += 1) : (nfail += 1)
        catch e
            @printf("[FAIL] %-40s  EXCEPTION: %s\n", id, sprint(showerror, e))
            nfail += 1
        end
    end
    @printf("\nSummary: %d pass, %d fail, %d skipped (no IC)\n", npass, nfail, nskip)
    return (npass, nfail, nskip)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(copy(ARGS))
end
