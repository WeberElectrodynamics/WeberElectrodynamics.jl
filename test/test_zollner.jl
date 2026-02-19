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
        sys = WeberSystem(2, 2)
        prob = WeberProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, 0.0];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.01,
            zollner_enabled = true, zollner_a = 0.1,
        )
        @test length(prob.kappas) == 1
        @test prob.kappas[1] ≈ 1.1  # 1 + a for unlike charges
        @test prob.zollner.enabled == true
        @test prob.zollner.a == 0.1
    end

    @testset "Kappa computation — like charges" begin
        sys = WeberSystem(2, 2)
        prob = WeberProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, 0.0];
            masses = [1.0, 1.0], charges = [1.0, 1.0], c = 10.0, dt = 0.01,
            zollner_enabled = true, zollner_a = 0.1,
        )
        @test length(prob.kappas) == 1
        @test prob.kappas[1] ≈ 1.0  # like charges: κ = 1 regardless of a
    end

    @testset "Kappa computation — 3-particle mixed" begin
        # Particles: +1, -1, +1
        # Pair (1,2): unlike → κ = 1+a
        # Pair (1,3): like   → κ = 1
        # Pair (2,3): unlike → κ = 1+a
        sys = WeberSystem(3, 2)
        a = 0.05
        prob = WeberProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, 0.0, 0.0, -1.0, 0.0], zeros(6);
            masses = [1.0, 1.0, 1.0], charges = [1.0, -1.0, 1.0],
            c = 10.0, dt = 0.01, zollner_enabled = true, zollner_a = a,
        )
        @test length(prob.kappas) == 3
        @test prob.kappas[1] ≈ 1.0 + a   # pair (1,2): unlike
        @test prob.kappas[2] ≈ 1.0        # pair (1,3): like
        @test prob.kappas[3] ≈ 1.0 + a   # pair (2,3): unlike
    end

    @testset "Params vector length" begin
        # For N particles: params length = 2N + 1 + N*(N-1)/2
        for N in [2, 3, 4]
            sys = WeberSystem(N, 2)
            q0 = zeros(N * 2)
            q0[1] = 1.0
            p0 = zeros(N * 2)
            charges = [(-1)^i * 1.0 for i in 1:N]
            prob = WeberProblem(
                sys, (0.0, 1.0), q0, p0;
                masses = ones(N), charges = charges, c = 10.0, dt = 0.01,
            )
            expected_len = 2N + 1 + N * (N - 1) ÷ 2
            @test length(prob.params) == expected_len
            @test length(prob.kappas) == N * (N - 1) ÷ 2
        end
    end

    @testset "Zöllner disabled ≡ κ=1 physics" begin
        # A disabled Zöllner with any a should give κ=1 for all pairs
        sys = WeberSystem(2, 2)
        prob_off = WeberProblem(
            sys, (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.01,
            zollner_enabled = false, zollner_a = 0.5,
        )
        @test all(k -> k ≈ 1.0, prob_off.kappas)
        # params tail should be all 1.0
        @test prob_off.params[end] ≈ 1.0
    end

    @testset "Integration runs with Zöllner enabled" begin
        # Smoke test: ensure the solver accepts and runs a Zöllner-modified problem
        sys = WeberSystem(2, 2)
        prob = WeberProblem(
            sys, (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner_enabled = true, zollner_a = 0.01,
            regularization_enabled = false,
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
        sys = WeberSystem(2, 2)

        # Disabled: Zöllner residual should be zero everywhere
        prob_off = WeberProblem(
            sys, (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner_enabled = false, regularization_enabled = false,
        )
        sol_off = solve(prob_off, SymmetricProjectionIntegrator())
        en_off = compute_energy_timeseries(sol_off)
        @test all(x -> x ≈ 0.0, en_off.total_zollner_residual)
        pair_off = en_off.pair_energies[(1, 2)]
        @test pair_off.kappa ≈ 1.0
        @test all(x -> x ≈ 0.0, pair_off.zollner_extra_potential)

        # Enabled with unlike charges: Zöllner residual should be nonzero
        prob_on = WeberProblem(
            sys, (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner_enabled = true, zollner_a = 0.05, regularization_enabled = false,
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
        sys = WeberSystem(2, 2)
        prob = WeberProblem(
            sys, (0.0, 0.5),
            [1.0, 0.0, -1.0, 0.0], [0.0, 0.5, 0.0, -0.5];
            masses = [1.0, 1.0], charges = [1.0, -1.0], c = 10.0, dt = 0.005,
            zollner_enabled = true, zollner_a = 0.01,
            regularization_enabled = true,
        )
        sol = solve(prob, SymmetricProjectionIntegrator())
        @test sol.retcode == :Success
        @test length(sol.t) > 1
    end

end
