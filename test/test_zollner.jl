@testset "Zöllner Electrogravitation" begin

    @testset "ZollnerOptions construction" begin
        # Default: disabled
        z = ZollnerOptions()
        @test z.enabled == false
        @test z.a == 0.0

        # Enabled with valid a
        z2 = ZollnerOptions(enabled = true, a = 0.01)
        @test z2.enabled == true
        @test z2.a == 0.01

        # Validation: a must be positive when enabled
        @test_throws AssertionError ZollnerOptions(enabled = true, a = 0.0)
        @test_throws AssertionError ZollnerOptions(enabled = true, a = -0.1)

        # Disabled with a=0 is valid
        @test ZollnerOptions(enabled = false, a = 0.0).enabled == false

        # Disabled with nonzero a is allowed (a is simply unused)
        @test ZollnerOptions(enabled = false, a = 0.5).a == 0.5
    end

    @testset "Kappa computation — unlike charges" begin
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, 0.0];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.01,
            zollner = ZollnerOptions(enabled = true, a = 0.1),
        )
        @test length(kappas(prob)) == 1
        @test kappas(prob)[1] ≈ 1.1  # 1 + a for unlike charges
        @test prob.zollner.enabled == true
        @test prob.zollner.a == 0.1
    end

    @testset "Kappa computation — like charges" begin
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, 0.0];
            masses = [1.0, 1.0], charges = [1.0, 1.0], c = 10.0, dt = 0.01,
            zollner = ZollnerOptions(enabled = true, a = 0.1),
        )
        @test length(kappas(prob)) == 1
        @test kappas(prob)[1] ≈ 1.0  # like charges: κ = 1 regardless of a
    end

    @testset "Kappa computation — 3-particle mixed" begin
        # Particles: +1, -1, +1
        # Pair (1,2): unlike → κ = 1+a
        # Pair (1,3): like   → κ = 1
        # Pair (2,3): unlike → κ = 1+a
        sys = HamiltonianSystem(3, 2)
        a = 0.05
        prob = HamiltonianProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, 0.0, 0.0, -1.0, 0.0], zeros(6);
            masses = [1.0, 1.0, 1.0], charges = [1.0, -1.0, 1.0],
            c = 10.0, dt = 0.01, zollner = ZollnerOptions(enabled = true, a = a),
        )
        @test length(kappas(prob)) == 3
        @test kappas(prob)[1] ≈ 1.0 + a   # pair (1,2): unlike
        @test kappas(prob)[2] ≈ 1.0        # pair (1,3): like
        @test kappas(prob)[3] ≈ 1.0 + a   # pair (2,3): unlike
    end

    @testset "Params vector length" begin
        # For N particles: params length = 2N + 1 + N*(N-1)/2
        for N in [2, 3, 4]
            sys = HamiltonianSystem(N, 2)
            q0 = zeros(N * 2)
            q0[1] = 1.0
            p0 = zeros(N * 2)
            charges = [(-1)^i * 1.0 for i in 1:N]
            prob = HamiltonianProblem(
                sys, (0.0, 1.0), q0, p0;
                masses = ones(N), charges = charges, c = 10.0, dt = 0.01,
            )
            expected_len = 2N + 1 + N * (N - 1) ÷ 2
            @test length(params(prob)) == expected_len
            @test length(kappas(prob)) == N * (N - 1) ÷ 2
            # params must literally end with the kappas values
            n_pairs = N * (N - 1) ÷ 2
            @test params(prob)[end-n_pairs+1:end] == kappas(prob)
        end
    end

    @testset "Zöllner disabled ≡ κ=1 physics" begin
        # A disabled Zöllner with any a should give κ=1 for all pairs
        sys = HamiltonianSystem(2, 2)
        prob_off = HamiltonianProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.01,
            zollner = ZollnerOptions(enabled = false, a = 0.5),
        )
        @test all(k -> k ≈ 1.0, kappas(prob_off))
        # params tail should be all 1.0
        @test params(prob_off)[end] ≈ 1.0
    end

    @testset "Integration runs with Zöllner enabled" begin
        # Smoke test: ensure the solver accepts and runs a Zöllner-modified problem
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner = ZollnerOptions(enabled = true, a = 0.01),
        )
        sol = solve(prob, SymmetricProjectionIntegrator())
        @test sol.retcode == :Success
        @test length(sol.t) > 1
        # Zöllner is conservative: energy should be well-conserved
        energy = compute_energy_timeseries(sol)
        E0 = energy.total_energy[1]
        if abs(E0) > 1e-10
            max_drift = maximum(abs.((energy.total_energy .- E0) ./ E0))
            @test max_drift < 1e-3
        end
    end

    @testset "Zöllner energy statistics" begin
        sys = HamiltonianSystem(2, 2)

        # Disabled: Zöllner residual should be zero everywhere
        prob_off = HamiltonianProblem(
            sys, (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
        )
        sol_off = solve(prob_off, SymmetricProjectionIntegrator())
        en_off = compute_energy_timeseries(sol_off)
        @test all(x -> x ≈ 0.0, en_off.total_zollner_residual)
        pair_off = en_off.pair_energies[(1, 2)]
        @test pair_off.kappa ≈ 1.0
        @test all(x -> x ≈ 0.0, pair_off.zollner_extra_potential)

        # Enabled with unlike charges: Zöllner residual should be nonzero
        prob_on = HamiltonianProblem(
            sys, (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner = ZollnerOptions(enabled = true, a = 0.05),
        )
        sol_on = solve(prob_on, SymmetricProjectionIntegrator())
        en_on = compute_energy_timeseries(sol_on)
        pair_on = en_on.pair_energies[(1, 2)]
        @test pair_on.kappa ≈ 1.05
        # For unlike charges, charge product is negative, so extra potential is negative (attractive)
        @test all(x -> x <= 0.0 || abs(x) < 1e-15, pair_on.zollner_extra_potential)
        @test !all(x -> x ≈ 0.0, en_on.total_zollner_residual)
    end

    @testset "Regularization + Zöllner compatibility" begin
        # Both features can be enabled simultaneously without errors.
        # Orbit is stable (no close encounter) so regularization stays idle,
        # but the infrastructure must accept both flags at once.
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.5),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.5, 0.0, -0.5];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner = ZollnerOptions(enabled = true, a = 0.01),
        )
        alg = RegularizedIntegrator(SymmetricProjectionIntegrator())
        sol = solve(prob, alg)
        @test sol.retcode == :Success
        @test length(sol.t) > 1
    end

    @testset "Zöllner kappas respected during regularized substeps" begin
        # Verify that κ ≠ 1 is actually used inside regularization substeps,
        # not just during unregularized steps.  Use a circular orbit at r0 < r_on
        # so regularization fires on every step.  Energy conservation confirms
        # the params_pair kappas path is exercised correctly.
        a = 0.05
        m1 = m2 = 1.0
        q1 = 0.1
        q2 = -0.1
        c = 10.0
        r0 = 0.1
        M = m1 + m2
        v_circ = sqrt(abs(q1 * q2) * M / (m1 * m2 * r0))

        sys = HamiltonianSystem(2, 2)
        q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
        p0 = [0.0, m1 * (-m2 / M) * v_circ, 0.0, m2 * (m1 / M) * v_circ]
        prob = HamiltonianProblem(
            sys, (0.0, 1.0), q0, p0;
            masses = [m1, m2], charges = [q1, q2], c = c, dt = 0.001,
            zollner = ZollnerOptions(enabled = true, a = a),
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            r_on = 0.15,
            r_off = 0.25,
        )
        sol = solve(prob, alg)
        @test sol.retcode == :Success
        # r0 < r_on so regularization fires on every step.
        @test sol.regularization.pair_steps > 0
        # Energy conservation validates the params_pair kappas path.
        energy = compute_energy_timeseries(sol)
        E0 = energy.total_energy[1]
        if abs(E0) > 1e-10
            max_drift = maximum(abs.((energy.total_energy .- E0) ./ E0))
            @test max_drift < 1e-3
        end
        # κ for this unlike-charge pair should be 1+a.
        @test kappas(sol.prob)[1] ≈ 1.0 + a
    end

end
