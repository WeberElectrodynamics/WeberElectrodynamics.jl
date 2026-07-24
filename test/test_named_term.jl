@testset "NamedTerm queries" begin
    @testset "convenience ctor exposes the Weber term" begin
        sys = HamiltonianSystem(3, 2)
        @test term_names(sys) == [:weber]
        @test has_term(sys, :weber)
        weber = get_term(sys, :weber)
        @test weber.name === :weber
        @test isequal(weber.H_symbolic, sys.hamiltonian_symbolic)
        @test weber.pair_decomposition !== nothing
    end

    @testset "pair_decomposition matches Weber components" begin
        sys = HamiltonianSystem(2, 2)
        weber = get_term(sys, :weber)
        q = [1.0, 0.2, -0.5, 0.3]
        p = [0.1, 0.2, -0.15, -0.1]
        params = [1.0, 2.0, 1.0, -1.0, 5.0]

        coulomb, velocity, rdot =
            WeberElectrodynamics.compute_pair_weber_components(
                q,
                p,
                1,
                2,
                params[1:2],
                params[3:4],
                params[5],
                2,
            )
        result = weber.pair_decomposition(1, 2, q, p, params)
        @test result.coulomb ≈ coulomb
        @test result.velocity ≈ velocity
        @test result.rdot ≈ rdot
        @test result.r ≈ sqrt((q[1] - q[3])^2 + (q[2] - q[4])^2)
    end

    @testset "pair_decomposition closures in 3D" begin
        sys = HamiltonianSystem(3, 3)
        weber = get_term(sys, :weber)
        q = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.5, 0.2]
        p = [0.1, 0.0, 0.0, -0.1, 0.05, 0.0, 0.0, -0.05, 0.02]
        params = [1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 8.0]
        for (i, j) in ((1, 2), (1, 3), (2, 3))
            expected = WeberElectrodynamics.compute_pair_weber_components(
                q,
                p,
                i,
                j,
                params[1:3],
                params[4:6],
                params[7],
                3,
            )
            result = weber.pair_decomposition(i, j, q, p, params)
            @test result.coulomb ≈ expected[1]
            @test result.velocity ≈ expected[2]
            @test result.rdot ≈ expected[3]
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
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 q1 q2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]
        kinetic =
            kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2)
        coulomb = coulomb_term(
            q_syms;
            charges = [q1, q2],
            n_particles = 2,
            dims = 2,
        )
        sys = HamiltonianSystem(
            kinetic + coulomb,
            q_syms,
            p_syms;
            param_symbols = [m1, m2, q1, q2, cc],
            t = tt,
            n_particles = 2,
            dims = 2,
            terms = [NamedTerm(:kinetic, kinetic), NamedTerm(:coulomb, coulomb)],
        )

        @test term_names(sys) == [:kinetic, :coulomb]
        @test has_term(sys, :kinetic) && has_term(sys, :coulomb)
        @test isequal(get_term(sys, :kinetic).H_symbolic, kinetic)
        @test isequal(get_term(sys, :coulomb).H_symbolic, coulomb)
    end
end
