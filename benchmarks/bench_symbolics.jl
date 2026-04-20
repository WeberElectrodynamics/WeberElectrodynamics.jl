#!/usr/bin/env julia
# Symbolics compile-cost benchmark for Phase 0 of the refactor.
#
# Purpose: measure symbolic-build + build_function compile times for the
# current Weber Hamiltonian path across a sweep of (n_particles, dims), and
# compare with `cse=true` to decide whether CSE should be a core default in
# the upcoming generic HamiltonianSystem builder.
#
# Run from the project root with:
#   julia --project=. benchmarks/bench_symbolics.jl
#
# Output: markdown table to stdout + benchmarks/results_phase0.md.

using Symbolics
using Printf

# Replicate the current Weber build path without pulling in the rest of the
# package — we want to measure the symbolic/compile cost in isolation.

function generate_phase_space_symbols(n::Int, d::Int)
    coord_names = [:x, :y, :z]
    p_names = [:px, :py, :pz]
    qs = [Symbolics.variable(Symbol(string(coord_names[k]) * string(i)))
          for i = 1:n for k = 1:d]
    ps = [Symbolics.variable(Symbol(string(p_names[k]) * string(i)))
          for i = 1:n for k = 1:d]
    return qs, ps
end

function generate_param_symbols(n::Int)
    ms = [Symbolics.variable(Symbol("m$i")) for i = 1:n]
    chs = [Symbolics.variable(Symbol("q$i")) for i = 1:n]
    c_var = Symbolics.variable(:c)
    kappas = [Symbolics.variable(Symbol("kappa_$(i)_$(j)"))
              for i = 1:n for j = (i+1):n]
    return ms, chs, c_var, kappas, vcat(ms, chs, [c_var], kappas)
end

@inline pair_index(i, j, n) = (i - 1) * (2n - i) ÷ 2 + (j - i)

function build_weber_hamiltonian(qs, ps, ms, chs, c_var, kappas, n, d)
    H = zero(eltype(qs))
    @inbounds for i = 1:n
        pstart = (i - 1) * d + 1
        pend = i * d
        pkin = sum(ps[pstart:pend] .^ 2)
        H = H + pkin / (2 * ms[i])
    end
    c2 = c_var^2
    @inbounds for i = 1:n
        for j = (i+1):n
            qis = (i - 1) * d + 1
            qjs = (j - 1) * d + 1
            pis = (i - 1) * d + 1
            pjs = (j - 1) * d + 1
            r2 = zero(eltype(qs))
            rv = zero(eltype(qs))
            for k = 1:d
                dq = qs[qis+k-1] - qs[qjs+k-1]
                dv = ps[pis+k-1] / ms[i] - ps[pjs+k-1] / ms[j]
                r2 = r2 + dq^2
                rv = rv + dq * dv
            end
            r = sqrt(r2)
            rdot = rv / r
            κ = kappas[pair_index(i, j, n)]
            k_scalar = κ * chs[i] * chs[j]
            H = H + k_scalar / r * (1 - rdot^2 / (2 * c2))
        end
    end
    return H
end

struct VariantResult
    n::Int
    dims::Int
    dof::Int
    n_params::Int
    variant::String
    t_build_H::Float64
    t_grad_q::Float64
    t_grad_p::Float64
    t_compile_dq::Float64
    t_compile_dp::Float64
    t_compile_H::Float64
    t_first_call::Float64
    t_per_call_ns::Float64
    expr_bytes_dq::Int
    expr_bytes_dp::Int
end

function run_variant(n::Int, d::Int; cse::Bool)
    qs, ps = generate_phase_space_symbols(n, d)
    ms, chs, c_var, kappas, params = generate_param_symbols(n)

    t_H = @elapsed H = build_weber_hamiltonian(qs, ps, ms, chs, c_var, kappas, n, d)

    t_dq_sym = @elapsed dq = [Symbolics.derivative(H, p) for p in ps]
    t_dp_sym = @elapsed dp = [-Symbolics.derivative(H, q) for q in qs]

    # build_function returns (out-of-place, in-place). We compile the in-place
    # variant which is what the hot path uses.
    build_kwargs = cse ? (:cse => true,) : ()

    t_dq_c = @elapsed dq_f = Symbolics.build_function(
        dq, qs, ps, params; expression = Val{false}, build_kwargs...
    )[2]

    t_dp_c = @elapsed dp_f = Symbolics.build_function(
        dp, qs, ps, params; expression = Val{false}, build_kwargs...
    )[2]

    t_H_c = @elapsed H_f = Symbolics.build_function(
        H, qs, ps, params; expression = Val{false}, build_kwargs...
    )

    # Expression size: stringify the *expression* (not compiled) and measure.
    # Uses a fresh build_function with expression=Val{true} to keep it textual.
    dq_expr_str = string(Symbolics.build_function(dq, qs, ps, params; expression = Val{true}, build_kwargs...)[2])
    dp_expr_str = string(Symbolics.build_function(dp, qs, ps, params; expression = Val{true}, build_kwargs...)[2])
    sz_dq = sizeof(dq_expr_str)
    sz_dp = sizeof(dp_expr_str)

    dof = n * d
    q0 = rand(dof)
    p0 = rand(dof)
    p_vals = vcat(ones(n), ones(n), [10.0], ones(n * (n - 1) ÷ 2))
    out = zeros(dof)

    # First-call latency (includes JIT specialization).
    t_first = @elapsed dq_f(out, q0, p0, p_vals)

    # Steady-state throughput: run enough calls to get a stable ns/call.
    n_iter = 10_000
    for _ = 1:10
        dq_f(out, q0, p0, p_vals)
    end
    t_loop = @elapsed for _ = 1:n_iter
        dq_f(out, q0, p0, p_vals)
    end
    ns_per = (t_loop / n_iter) * 1e9

    return VariantResult(
        n, d, dof, length(p_vals),
        cse ? "cse=true" : "baseline",
        t_H, t_dq_sym, t_dp_sym,
        t_dq_c, t_dp_c, t_H_c,
        t_first, ns_per,
        sz_dq, sz_dp,
    )
end

function format_row(r::VariantResult)
    total_compile = r.t_compile_dq + r.t_compile_dp + r.t_compile_H
    @sprintf(
        "| n=%-2d | d=%d | dof=%-3d | %-9s | H %5.2fs | ∂q %5.2fs | ∂p %5.2fs | bf_dq %5.2fs | bf_dp %5.2fs | bf_H %5.2fs | Σ %6.2fs | first %5.3fs | %8.1f ns | dq %7d B | dp %7d B |",
        r.n, r.dims, r.dof, r.variant,
        r.t_build_H, r.t_grad_q, r.t_grad_p,
        r.t_compile_dq, r.t_compile_dp, r.t_compile_H, total_compile,
        r.t_first_call, r.t_per_call_ns,
        r.expr_bytes_dq, r.expr_bytes_dp,
    )
end

function main()
    sweep = [(3, 2), (3, 3), (5, 2), (5, 3), (8, 2), (8, 3), (12, 2), (12, 3)]

    results = VariantResult[]

    println("# Phase 0 — Symbolics compile benchmark")
    println()
    println("Julia: ", VERSION, "  |  Symbolics: v", pkgversion(Symbolics))
    println()

    # Warm up Symbolics so the first real measurement isn't polluted by
    # one-off package compilation.
    print("Warming up Symbolics... ")
    flush(stdout)
    _ = run_variant(2, 2; cse = false)
    println("done.")
    println()

    for (n, d) in sweep, cse in (false, true)
        print("Running n=$n, d=$d, ", cse ? "cse=true" : "baseline", " ... ")
        flush(stdout)
        r = run_variant(n, d; cse = cse)
        push!(results, r)
        println("Σ compile = ",
                round(r.t_compile_dq + r.t_compile_dp + r.t_compile_H; digits = 2), "s")
    end

    println()
    println("## Results")
    println()
    for r in results
        println(format_row(r))
    end

    return results
end

results = main()

# Also write machine-readable results for the report script.
open(joinpath(@__DIR__, "results_phase0.raw.txt"), "w") do io
    for r in results
        println(io, join([
            r.n, r.dims, r.dof, r.n_params, r.variant,
            r.t_build_H, r.t_grad_q, r.t_grad_p,
            r.t_compile_dq, r.t_compile_dp, r.t_compile_H,
            r.t_first_call, r.t_per_call_ns,
            r.expr_bytes_dq, r.expr_bytes_dp,
        ], "\t"))
    end
end

println()
println("Raw results written to benchmarks/results_phase0.raw.txt")
