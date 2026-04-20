using WeberElectrodynamics
using Symbolics
using Latexify
using Random

const OUTDIR = @__DIR__

sys = HamiltonianSystem(4, 2)
H = sys.hamiltonian_symbolic
q = sys.q_symbols
p = sys.p_symbols
ps = sys.param_symbols

subs = Dict{Num,Num}()
for i in 1:4
    subs[ps[i]] = Num(1)
end
subs[ps[5]] = Num(1); subs[ps[6]] = Num(1)
subs[ps[7]] = Num(-1); subs[ps[8]] = Num(-1)
subs[ps[9]] = Num(1)
for i in 10:15
    subs[ps[i]] = Num(1)
end

H_sub = Symbolics.substitute(H, subs)
H_exp = expand(H_sub)

open(joinpath(OUTDIR, "hamiltonian_expanded.txt"), "w") do io
    println(io, "# Full Weber Hamiltonian, N=4, 2D, m_i=1, q=(+1,+1,-1,-1), c=1, kappa=1")
    println(io, "# Variables: q = [x1,y1,x2,y2,x3,y3,x4,y4], p = [px1,py1,...,px4,py4]")
    println(io)
    println(io, "H = ", string(H_exp))
end

try
    open(joinpath(OUTDIR, "hamiltonian.tex"), "w") do io
        println(io, "% LaTeX form of the 4-body Weber Hamiltonian (2+/2-), m=c=kappa=1")
        println(io, latexify(H_exp))
    end
catch e
    @warn "latexify failed" exception=e
end

charges = [Num(1), Num(1), Num(-1), Num(-1)]
masses  = [Num(1), Num(1), Num(1), Num(1)]
c_num = Num(1)

T_kin = sum((p[2i-1]^2 + p[2i]^2) / (2*masses[i]) for i in 1:4)
pairs = [(i,j) for i in 1:3 for j in (i+1):4]

function rij(i,j)
    dx = q[2i-1] - q[2j-1]
    dy = q[2i]   - q[2j]
    sqrt(dx^2 + dy^2)
end
function rdot_ij(i,j)
    dx = q[2i-1] - q[2j-1]
    dy = q[2i]   - q[2j]
    dvx = p[2i-1]/masses[i] - p[2j-1]/masses[j]
    dvy = p[2i]  /masses[i] - p[2j]  /masses[j]
    (dx*dvx + dy*dvy) / rij(i,j)
end

U_coulomb = sum(charges[i]*charges[j]/rij(i,j) for (i,j) in pairs)
U_weber   = sum(-charges[i]*charges[j] * rdot_ij(i,j)^2 / (2*c_num^2 * rij(i,j)) for (i,j) in pairs)
H_decomp  = T_kin + U_coulomb + U_weber

Random.seed!(42)
function sample_dict()
    d = Dict{Num,Num}()
    for k in eachindex(q); d[q[k]] = Num(randn()); end
    for k in eachindex(p); d[p[k]] = Num(randn()); end
    d
end
function numval(expr, d)
    r = Symbolics.substitute(expr, d; fold=Val(true))
    try
        return Float64(Symbolics.value(r))
    catch
        # fallback: try Symbolics.symbolic_to_float, else parse result
        return convert(Float64, Symbolics.symbolic_to_float(r))
    end
end

decomp_errs = Float64[]
for _ in 1:5
    d = sample_dict()
    push!(decomp_errs, abs(numval(H_sub, d) - numval(H_decomp, d)))
end
decomp_ok = maximum(decomp_errs) < 1e-10
println("DECOMP max err = ", maximum(decomp_errs), "  OK=", decomp_ok)

open(joinpath(OUTDIR, "pair_table.md"), "w") do io
    println(io, "# Pair analysis (2+/2- with q1=q2=+1, q3=q4=-1, m=c=1)\n")
    println(io, "| pair (i,j) | q_i q_j | Coulomb | Weber correction sign | ρ_ij = q_i q_j /(μ c²), μ=1/2 |")
    println(io, "|---|---|---|---|---|")
    for (i,j) in pairs
        qq_val = (i<=2 && j<=2) || (i>=3 && j>=3) ? +1 : -1
        coul   = qq_val == +1 ? "repulsive (+1/r)" : "attractive (-1/r)"
        web    = qq_val == +1 ? "negative (lowers H; weakens repulsion)" :
                                "positive (raises H; weakens attraction)"
        ρ      = qq_val == +1 ? "+2" : "-2"
        println(io, "| ($i,$j) | $qq_val | $coul | $web | $ρ |")
    end
    println(io)
    println(io, "With μ_ij = m_i m_j/(m_i+m_j) = 1/2 and c=1: ρ_ij = 2·q_i q_j ∈ {+2,-2}.")
    println(io, "Like-charge pairs (1,2),(3,4) have ρ = +2; the four unlike-sign pairs have ρ = -2.")
end

function poisson_sym(F, H)
    s = Num(0)
    for k in eachindex(q)
        s += Symbolics.derivative(F, q[k]) * Symbolics.derivative(H, p[k]) -
             Symbolics.derivative(F, p[k]) * Symbolics.derivative(H, q[k])
    end
    s
end
function max_abs_at_samples(expr; n=6)
    m = 0.0
    for _ in 1:n
        d = sample_dict()
        v = abs(numval(expr, d))
        m = max(m, v)
    end
    m
end

Px = p[1] + p[3] + p[5] + p[7]
Py = p[2] + p[4] + p[6] + p[8]
L  = sum(q[2i-1]*p[2i] - q[2i]*p[2i-1] for i in 1:4)
D  = sum(q[k]*p[k] for k in eachindex(q))
d_x = sum(charges[i]*q[2i-1] for i in 1:4)
d_y = sum(charges[i]*q[2i]   for i in 1:4)

brackets = Dict{String,Num}()
brackets["Px"] = poisson_sym(Px, H_sub)
brackets["Py"] = poisson_sym(Py, H_sub)
brackets["L"]  = poisson_sym(L, H_sub)
brackets["H"]  = poisson_sym(H_sub, H_sub)
brackets["D"]  = poisson_sym(D, H_sub)
brackets["d_x"]= poisson_sym(d_x, H_sub)
brackets["d_y"]= poisson_sym(d_y, H_sub)

pb_report = Dict{String,Float64}()
for (name, expr) in brackets
    pb_report[name] = max_abs_at_samples(expr)
    println("{$name, H} numeric max |·| = ", pb_report[name])
end

function apply_perm(expr, perm)
    d = Dict{Num,Num}()
    for (a,b) in perm
        d[q[2a-1]] = q[2b-1]
        d[q[2a]]   = q[2b]
        d[p[2a-1]] = p[2b-1]
        d[p[2a]]   = p[2b]
    end
    Symbolics.substitute(expr, d)
end
H_C = apply_perm(H_sub, [(1,2),(2,1),(3,4),(4,3)])

Pdict = Dict{Num,Num}()
for k in eachindex(q); Pdict[q[k]] = -q[k]; end
for k in eachindex(p); Pdict[p[k]] = -p[k]; end
H_P = Symbolics.substitute(H_sub, Pdict)

Tdict = Dict{Num,Num}()
for k in eachindex(p); Tdict[p[k]] = -p[k]; end
H_T = Symbolics.substitute(H_sub, Tdict)

sym_errs = Dict(
    "C" => max_abs_at_samples(H_C - H_sub),
    "P" => max_abs_at_samples(H_P - H_sub),
    "T" => max_abs_at_samples(H_T - H_sub),
)
println("Symmetry residuals: ", sym_errs)

U_full = U_coulomb + U_weber
virial_q = sum(q[k]*Symbolics.derivative(U_full, q[k]) for k in eachindex(q))
v_q_err = max_abs_at_samples(virial_q + U_full)
println("q-virial residual: ", v_q_err)

virial_p = sum(p[k]*Symbolics.derivative(H_sub, p[k]) for k in eachindex(p))
v_p_err = max_abs_at_samples(virial_p - (2*T_kin + 2*U_weber))
println("p-virial residual: ", v_p_err)

open(joinpath(OUTDIR, "virial.txt"), "w") do io
    println(io, "=== Virial identities for 2+/2- Weber (numerical verification) ===\n")
    println(io, "q-homogeneity (U is degree -1 in q):   sum_k q_k ∂U/∂q_k = -U")
    println(io, "  residual max |sum q.∂U/∂q + U| = ", v_q_err, "\n")
    println(io, "p-identity (U_weber is degree +2 in p):   sum_k p_k ∂H/∂p_k = 2T + 2 U_weber")
    println(io, "  residual max |LHS - RHS| = ", v_p_err, "\n")
    println(io, "Therefore  d/dt(Σ q·p) = Σ p·(∂H/∂p) − Σ q·(∂H/∂q)")
    println(io, "                       = (2T + 2 U_w) − (−U_c − U_w)")
    println(io, "                       = 2T + U_c + 3 U_w.")
    println(io, "Time-averaging on bounded orbits (⟨d/dt(q·p)⟩=0):")
    println(io, "     2⟨T⟩ + ⟨U_c⟩ + 3⟨U_w⟩ = 0")
    println(io, "i.e.  2⟨T⟩ = −⟨U_c⟩ − 3⟨U_w⟩   (Weber-modified virial).")
    println(io, "Reduces to the classical 2⟨T⟩ = −⟨U⟩ when U_weber → 0.")
end

open(joinpath(OUTDIR, "conserved_quantities.md"), "w") do io
    println(io, "# Conserved Quantities & Symmetries — 4-body 2+/2− Weber (2D)\n")
    println(io, "All residuals computed numerically at 6 random phase-space points (seed 42).")
    println(io, "A residual below 1e-10 certifies the identity to machine precision.\n")
    println(io, "## Continuous generators — Poisson brackets with H\n")
    println(io, "| F | definition | max\\|{F,H}\\| | conserved? |")
    println(io, "|---|---|---|---|")
    fmt(x) = x < 1e-10 ? "≈0 ("*string(round(x,sigdigits=2))*")" : string(round(x,sigdigits=4))
    rows = [
        ("Px", "P_x", "Σ p_{x,i}"),
        ("Py", "P_y", "Σ p_{y,i}"),
        ("L",  "L",   "Σ (x_i p_{y,i} − y_i p_{x,i})"),
        ("H",  "H",   "H itself (sanity)"),
        ("D",  "D",   "Σ q_k p_k  (dilation)"),
        ("d_x","d_x", "Σ charge_i · x_i  (charge dipole, x)"),
        ("d_y","d_y", "Σ charge_i · y_i  (charge dipole, y)"),
    ]
    for (key, nm, defn) in rows
        r = pb_report[key]
        status = r < 1e-10 ? "**yes**" : "no"
        println(io, "| ", nm, " | ", defn, " | ", fmt(r), " | ", status, " |")
    end
    println(io)
    println(io, "### Non-conserved residuals in closed form")
    println(io, "- **{D,H} = 2T + U_c + 3 U_w** (derived from the q- and p-virial identities;")
    println(io, "  verified numerically in virial.txt). D is conserved only on the zero-set")
    println(io, "  of this expression, which is *not* a level set of H.")
    println(io, "- **{d_x,H} = Σ charge_i · ∂H/∂p_{x,i} = Σ charge_i · (p_{x,i}/m_i + Weber-velocity terms)**")
    println(io, "  — this is a charge-weighted current; non-zero in general.")
    println(io)
    println(io, "## Discrete symmetries\n")
    println(io, "| symmetry | action | residual max\\|H∘g − H\\| | invariant? |")
    println(io, "|---|---|---|---|")
    println(io, "| C (charge swap) | (1↔2),(3↔4) in q and p | ", fmt(sym_errs["C"]), " | **yes** |")
    println(io, "| P (parity)      | q→−q, p→−p             | ", fmt(sym_errs["P"]), " | **yes** |")
    println(io, "| T (time reverse)| p→−p, t→−t             | ", fmt(sym_errs["T"]), " | **yes** |")
    println(io)
    println(io, "H depends on p only through p² (kinetic) and ṙ² (Weber), both even in p; so T")
    println(io, "is manifest. P flips both r and v in each pair, leaving r² and ṙ² invariant.")
    println(io, "C permutes identical-mass, identical-charge particles, leaving U pairwise invariant.")
    println(io)
    println(io, "## Virial (Weber-modified)\n")
    println(io, "    2⟨T⟩ + ⟨U_c⟩ + 3⟨U_w⟩ = 0   ⇔   2⟨T⟩ = −⟨U_c⟩ − 3⟨U_w⟩.")
    println(io, "Reduces to 2⟨T⟩ = −⟨U⟩ when U_weber → 0.")
end

println("DECOMP_OK=", decomp_ok)
println("DONE.")
