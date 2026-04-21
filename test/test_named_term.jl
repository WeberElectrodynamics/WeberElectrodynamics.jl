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
        @test weber.pair_decomposition !== nothing
        @test zol.pair_decomposition !== nothing
    end

    @testset "pair_decomposition closures match compute_pair_weber_components" begin
        sys = HamiltonianSystem(2, 2)
        weber = get_term(sys, :weber)
        zol = get_term(sys, :zollner)

        q = [1.0, 0.2, -0.5, 0.3]
        p = [0.1, 0.2, -0.15, -0.1]
        # params = [m1, m2, q1, q2, c, κ12]
        params = [1.0, 2.0, 1.0, -1.0, 5.0, 1.3]

        masses = params[1:2]
        charges = params[3:4]
        c_val = params[5]
        κ = params[6]

        coulomb, velocity, rdot, zollner_extra =
            WeberElectrodynamics.compute_pair_weber_components(
                q,
                p,
                1,
                2,
                masses,
                charges,
                c_val,
                2,
                κ,
            )

        wr = weber.pair_decomposition(1, 2, q, p, params)
        zr = zol.pair_decomposition(1, 2, q, p, params)

        @test wr.coulomb ≈ coulomb
        @test wr.velocity ≈ velocity
        @test wr.rdot ≈ rdot
        @test wr.r ≈ sqrt((q[1] - q[3])^2 + (q[2] - q[4])^2)

        @test zr.zollner_extra ≈ zollner_extra
        @test zr.r ≈ wr.r
        @test zr.rdot ≈ rdot

        # κ = 1 → Zöllner contribution vanishes identically.
        params_k1 = [1.0, 2.0, 1.0, -1.0, 5.0, 1.0]
        zr1 = zol.pair_decomposition(1, 2, q, p, params_k1)
        @test zr1.zollner_extra == 0.0
    end

    @testset "pair_decomposition closures in 3D" begin
        sys = HamiltonianSystem(3, 3)
        weber = get_term(sys, :weber)
        zol = get_term(sys, :zollner)

        q = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.5, 0.2]
        p = [0.1, 0.0, 0.0, -0.1, 0.05, 0.0, 0.0, -0.05, 0.02]
        # params = [m1, m2, m3, q1, q2, q3, c, κ12, κ13, κ23]
        params = [1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 8.0, 1.1, 1.0, 1.2]

        masses = params[1:3]
        charges = params[4:6]
        c_val = params[7]

        for (i, j, κ) in ((1, 2, params[8]), (1, 3, params[9]), (2, 3, params[10]))
            coulomb, velocity, rdot, zollner_extra =
                WeberElectrodynamics.compute_pair_weber_components(
                    q,
                    p,
                    i,
                    j,
                    masses,
                    charges,
                    c_val,
                    3,
                    κ,
                )

            wr = weber.pair_decomposition(i, j, q, p, params)
            zr = zol.pair_decomposition(i, j, q, p, params)

            @test wr.coulomb ≈ coulomb
            @test wr.velocity ≈ velocity
            @test wr.rdot ≈ rdot
            @test zr.zollner_extra ≈ zollner_extra
        end
    end

    @testset "generic ctor defaults to :hamiltonian" begin
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 tt
        q = [x1, y1, x2, y2]
        p = [px1, py1, px2, py2]
        H = sum(p .^ 2) / 2
        sys = HamiltonianSystem(
            H,
            q,
            p;
            param_symbols = [m1, m2],
            t = tt,
            n_particles = 2,
            dims = 2,
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
            q_syms,
            p_syms;
            masses = [m1, m2, m3],
            charges = [q1, q2, q3],
            c = cc,
            kappas = [1.0, 1.0, 1.0],
            n_particles = 3,
            dims = 2,
        )
        H_corr = zollner_term(
            q_syms,
            p_syms;
            masses = [m1, m2, m3],
            charges = [q1, q2, q3],
            c = cc,
            kappas = [k12, k13, k23],
            n_particles = 3,
            dims = 2,
        )
        sys = HamiltonianSystem(
            H_pure + H_corr,
            q_syms,
            p_syms;
            param_symbols = param_syms,
            t = tt,
            n_particles = 3,
            dims = 2,
            terms = [NamedTerm(:weber, H_pure), NamedTerm(:zollner, H_corr)],
        )

        @test term_names(sys) == [:weber, :zollner]
        @test has_term(sys, :weber) && has_term(sys, :zollner)
        @test isequal(get_term(sys, :weber).H_symbolic, H_pure)
        @test isequal(get_term(sys, :zollner).H_symbolic, H_corr)
    end
end
