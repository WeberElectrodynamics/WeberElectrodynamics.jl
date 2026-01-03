using Symbolics

struct SymbolicHamiltonian
    H_sym::Num
    qdot_sym::Vector{Num}
    pdot_sym::Vector{Num}
    param_vars::Vector{Num}
    q_vars::Vector{Num}
    p_vars::Vector{Num}
end

struct HamiltonianVectorField
    qdot_func::Function
    pdot_func::Function
end

function create_phase_space_variables(n_particles::UInt8, dims::UInt8)::Tuple{Vector{Symbol}, Vector{Symbol}}
    coord_names = [:x, :y, :z]
    momentum_names = [:px, :py, :pz]

    q_syms = Symbol[]
    p_syms = Symbol[]

    for i in 1:n_particles
        for d in 1:dims
            push!(q_syms, Symbol(string(coord_names[d]) * string(i)))
            push!(p_syms, Symbol(string(momentum_names[d]) * string(i)))
        end
    end

    return (q_syms, p_syms)
end

function symbolize(H_func::Function, q_syms::Vector{Symbol}, p_syms::Vector{Symbol}, param_names::Vector{Symbol})::SymbolicHamiltonian
    q_vars = [Symbolics.variable(q_sym) for q_sym in q_syms]
    p_vars = [Symbolics.variable(p_sym) for p_sym in p_syms]
    param_vars = [Symbolics.variable(param_name) for param_name in param_names]
    t = Symbolics.variable(:t)

    H_sym = H_func(q_vars, p_vars, param_vars, t)

    qdot_sym = [Symbolics.derivative(H_sym, p_vars[i]) for i in eachindex(p_vars)]
    pdot_sym = [-Symbolics.derivative(H_sym, q_vars[i]) for i in eachindex(q_vars)]

    return SymbolicHamiltonian(
        H_sym,
        qdot_sym,
        pdot_sym,
        param_vars,
        q_vars,
        p_vars,
    )
end

function compile(sh::SymbolicHamiltonian)::HamiltonianVectorField
    qdot_func = build_function(sh.qdot_sym, sh.q_vars, sh.p_vars, sh.param_vars, expression=Val{false})[2]
    pdot_func = build_function(sh.pdot_sym, sh.q_vars, sh.p_vars, sh.param_vars, expression=Val{false})[2]

    return HamiltonianVectorField(
        qdot_func,
        pdot_func,
    )
end
