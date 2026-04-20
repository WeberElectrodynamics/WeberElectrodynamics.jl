@testset "NamedTerm queries" begin
    @testset "convenience ctor decomposes into weber + zollner" begin
        sys = HamiltonianSystem(3, 2)
        @test term_names(sys) == [:weber, :zollner]
        @test has_term(sys, :weber)
        @test has_term(sys, :zollner)
        weber = get_term(sys, :weber)
        zol = get_term(sys, :zollner)
        @test weber.name === :weber
        @test zol.name === :zollner
        # Sum of the two symbolic term Hamiltonians reconstructs the full system H.
        @test isequal(weber.H_symbolic + zol.H_symbolic, sys.hamiltonian_symbolic)
        @test weber.pair_decomposition === nothing
        @test zol.pair_decomposition === nothing
    end

    @testset "generic ctor defaults to :hamiltonian" begin
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 tt
        q = [x1, y1, x2, y2]
        p = [px1, py1, px2, py2]
        H = sum(p .^ 2) / 2
        sys = HamiltonianSystem(
            H, q, p;
            param_symbols = [m1, m2], t = tt,
            n_particles = 2, dims = 2,
        )
        @test term_names(sys) == [:hamiltonian]
        @test_throws KeyError get_term(sys, :missing)
    end

    @testset "generic ctor accepts an explicit terms vector" begin
        @variables x1 y1 x2 y2 x3 y3 px1 py1 px2 py2 px3 py3
        @variables m1 m2 m3 q1 q2 q3 cc k12 k13 k23 tt
        q_syms = [x1, y1, x2, y2, x3, y3]
        p_syms = [px1, py1, px2, py2, px3, py3]
        param_syms = [m1, m2, m3, q1, q2, q3, cc, k12, k13, k23]

        H_pure = weber_term(
            q_syms, p_syms;
            masses = [m1, m2, m3], charges = [q1, q2, q3], c = cc,
            kappas = [1.0, 1.0, 1.0], n_particles = 3, dims = 2,
        )
        H_corr = zollner_term(
            q_syms, p_syms;
            masses = [m1, m2, m3], charges = [q1, q2, q3], c = cc,
            kappas = [k12, k13, k23], n_particles = 3, dims = 2,
        )
        sys = HamiltonianSystem(
            H_pure + H_corr, q_syms, p_syms;
            param_symbols = param_syms, t = tt,
            n_particles = 3, dims = 2,
            terms = [NamedTerm(:weber, H_pure), NamedTerm(:zollner, H_corr)],
        )

        @test term_names(sys) == [:weber, :zollner]
        @test has_term(sys, :weber) && has_term(sys, :zollner)
        @test isequal(get_term(sys, :weber).H_symbolic, H_pure)
        @test isequal(get_term(sys, :zollner).H_symbolic, H_corr)
    end
end
