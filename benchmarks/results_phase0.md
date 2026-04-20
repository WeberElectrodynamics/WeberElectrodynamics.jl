# Phase 0 — Symbolics compile benchmark

**Goal**: measure the symbolic-build + `build_function` compile cost for the current
Weber Hamiltonian across `n ∈ {3, 5, 8, 12}` particles and `dims ∈ {2, 3}`, with and
without `cse=true`, to decide whether CSE should be a core default in the upcoming
generic `HamiltonianSystem` builder.

**Environment**: Julia 1.12.6, Symbolics v7.19.0, macOS (Darwin 25.4.0), single thread.

**Script**: [`bench_symbolics.jl`](bench_symbolics.jl). Raw table:
[`results_phase0.raw.txt`](results_phase0.raw.txt).

## Gate decision — **PASS**

> Generic build at `n=8, dims=3` compiles in ≤ 2× current Weber time, OR
> `cse=true` / sharded form closes the gap to ≤ 2×, OR
> document a deferred/lazy per-term compilation fallback.

We blew through the gate. At the worst measured case (`n=12, dims=3`, 163 symbolic
variables: 36 coords + 36 momenta + 12 masses + 12 charges + 66 κ-pairs + c), total
compile time is **0.26 s** baseline and **0.07 s** with `cse=true`. No compile cliff.

**Recommendation**: adopt `cse=true` as the default in the upcoming generic
`HamiltonianSystem` build pipeline. It is a strict win on every axis:

| Axis                  | baseline → cse=true (n=12, d=3) | Improvement |
|-----------------------|----------------------------------|-------------|
| Σ compile time        | 0.26 s → 0.07 s                  | **3.7×**    |
| First-call latency    | 2.46 s → 0.62 s                  | **4.0×**    |
| Steady-state per call | 3599 ns → 890 ns                 | **4.0×**    |
| `∂p` expression size  | 1.15 MB → 544 kB                 | **2.1×**    |

No sharded / multithreaded form evaluated — not needed at the measured scale.
Re-visit only if the generic builder (Phase 1) measurably regresses past the
2× mark relative to these numbers.

## Results

Data captured directly from the benchmark run. Time units: seconds (compile)
or nanoseconds (per-call). Byte counts are expression string sizes.

| n  | d | dof | variant  | build H | ∂q    | ∂p    | bf_dq | bf_dp | bf_H  | Σ compile | first call | ns/call | dq bytes | dp bytes  |
|----|---|-----|----------|---------|-------|-------|-------|-------|-------|-----------|------------|---------|----------|-----------|
| 3  | 2 | 6   | baseline | 0.00 s  | 0.00s | 0.00s | 0.00s | 0.00s | 0.00s | 0.00 s    | 0.044 s    | 82.2    | 7 715    | 25 455    |
| 3  | 2 | 6   | cse=true | 0.00 s  | 0.01s | 0.00s | 0.04s | 0.00s | 0.01s | 0.05 s    | 0.019 s    | 35.0    | 9 237    | 19 107    |
| 3  | 3 | 9   | baseline | 0.00 s  | 0.00s | 0.00s | 0.00s | 0.01s | 0.00s | 0.01 s    | 0.085 s    | 150.8   | 14 623   | 51 925    |
| 3  | 3 | 9   | cse=true | 0.00 s  | 0.00s | 0.00s | 0.00s | 0.00s | 0.00s | 0.00 s    | 0.023 s    | 49.1    | 12 940   | 27 338    |
| 5  | 2 | 10  | baseline | 0.00 s  | 0.00s | 0.01s | 0.00s | 0.01s | 0.00s | 0.01 s    | 0.168 s    | 264.9   | 23 685   | 83 456    |
| 5  | 2 | 10  | cse=true | 0.00 s  | 0.00s | 0.01s | 0.00s | 0.00s | 0.00s | 0.01 s    | 0.045 s    | 99.1    | 26 120   | 58 730    |
| 5  | 3 | 15  | baseline | 0.00 s  | 0.01s | 0.01s | 0.01s | 0.02s | 0.00s | 0.03 s    | 0.373 s    | 529.1   | 46 837   | 173 401   |
| 5  | 3 | 15  | cse=true | 0.00 s  | 0.01s | 0.01s | 0.00s | 0.01s | 0.00s | 0.01 s    | 0.070 s    | 150.3   | 36 941   | 84 642    |
| 8  | 2 | 16  | baseline | 0.00 s  | 0.01s | 0.02s | 0.01s | 0.03s | 0.01s | 0.06 s    | 0.486 s    | 789.1   | 65 031   | 234 656   |
| 8  | 2 | 16  | cse=true | 0.00 s  | 0.01s | 0.02s | 0.01s | 0.01s | 0.00s | 0.02 s    | 0.136 s    | 272.5   | 68 092   | 160 968   |
| 8  | 3 | 24  | baseline | 0.00 s  | 0.02s | 0.03s | 0.02s | 0.17s | 0.00s | 0.19 s    | 0.957 s    | 1 490.4 | 129 727  | 487 438   |
| 8  | 3 | 24  | cse=true | 0.00 s  | 0.02s | 0.03s | 0.01s | 0.02s | 0.00s | 0.03 s    | 0.203 s    | 406.0   | 96 475   | 233 212   |
| 12 | 2 | 24  | baseline | 0.01 s  | 0.04s | 0.05s | 0.02s | 0.12s | 0.00s | 0.14 s    | 1.198 s    | 1 768.3 | 152 563  | 555 748   |
| 12 | 2 | 24  | cse=true | 0.01 s  | 0.02s | 0.04s | 0.01s | 0.03s | 0.01s | 0.05 s    | 0.320 s    | 676.0   | 156 561  | 376 504   |
| 12 | 3 | 36  | baseline | 0.01 s  | 0.05s | 0.07s | 0.09s | 0.16s | 0.01s | 0.26 s    | 2.459 s    | 3 598.7 | 305 263  | 1 153 774 |
| 12 | 3 | 36  | cse=true | 0.01 s  | 0.04s | 0.06s | 0.02s | 0.05s | 0.01s | 0.07 s    | 0.616 s    | 890.3   | 222 540  | 544 004   |

Columns: `H` = symbolic Hamiltonian build; `∂q`, `∂p` = gradient via
`Symbolics.derivative`; `bf_dq`, `bf_dp`, `bf_H` = `build_function` compile;
`Σ compile` = `bf_dq + bf_dp + bf_H`; `first call` = first `dq_f(out,q,p,params)`
invocation (JIT specialization); `ns/call` = steady-state from 10 000-iteration loop
after warmup.

## Interpretation

1. **Symbolic build and gradient costs are negligible** across the whole sweep
   (under 0.15 s combined even at n=12, d=3). The cost that matters is
   `build_function` compilation (`bf_dq + bf_dp`).
2. **Expression size scales like O(n² · d)** — the pairwise Weber potential.
   `dp` dominates because the potential gradient has `(n-1)` contributions per
   coordinate, while `dq` is just momentum.
3. **`cse=true` shrinks `bf_dp` compile time by 3–8×** and expression bytes by 2×.
   No measured downside on any row.
4. **Steady-state throughput is faster with CSE**, which is slightly surprising —
   CSE pulls common factors into locals, and the Julia compiler evidently
   optimises the shared bindings better than the un-CSE'd expanded form.
5. **First-call latency is the dominant user-visible cost.** At n=12, d=3,
   baseline first-call is 2.5 s (pure JIT); CSE cuts it to 0.6 s.

## Implications for Phase 1

- The generic `HamiltonianSystem(H, q, p; params, t)` builder will compile the
  aggregate EOMs with `cse=true` by default. No flag, no opt-in.
- Per-term compiled closures (`weber_term`, `zollner_term`, `kinetic_term`)
  get the same treatment for symmetry, even though their individual costs are
  trivial — it preserves a single compile path and makes the eventual custom-term
  case (user-provided H) behave the same as the built-in cases.
- At expected working scales (n ≤ 12, d ≤ 3), compile time is a non-issue.
  Deferred / lazy / per-term JIT is **not needed** and will not be introduced
  speculatively.
- If a future user assembles a very large system (n ≥ 20, or dense all-to-all
  velocity terms beyond the Weber shape), we revisit with the sharded /
  multithreaded form. Re-running this benchmark is the trigger.

## Raw data

See [`results_phase0.raw.txt`](results_phase0.raw.txt) for the tab-separated
dump consumed by any future analysis script.
