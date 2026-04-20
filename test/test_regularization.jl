@testset "Regularization" begin
    function make_orbit_problem(dims::Int;
        dt::Float64,
        t_end::Float64,
        v_scale::Float64,
        regularization_enabled::Bool,
        backend::Symbol,
        warn_on_fallback::Bool = false,
        r_on::Float64 = 0.6,
        r_off::Float64 = 0.9,
        max_substeps::Int = 512,
    )
        m1, m2 = 1.0, 0.1
        q1, q2 = sqrt(0.1), -sqrt(0.1)
        c = 4.0

        r0 = 2.0
        M = m1 + m2
        k = q1 * q2
        v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))

        if dims == 1
            q0 = [-0.2, 0.2]
            p0 = [0.0, 0.0]
            masses = [1.0, 1.0]
            charges = [0.2, -0.2]
        elseif dims == 2
            q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
            p0 = [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, m2 * (m1 / M * v_circ * v_scale)]
            masses = [m1, m2]
            charges = [q1, q2]
        else
            q0 = [-m2 / M * r0, 0.0, 0.0, m1 / M * r0, 0.0, 0.0]
            p0 = [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, 0.0, m2 * (m1 / M * v_circ * v_scale), 0.0]
            masses = [m1, m2]
            charges = [q1, q2]
        end

        system = HamiltonianSystem(2, dims)
        prob = HamiltonianProblem(
            system,
            (0.0, t_end),
            q0,
            p0;
            masses = masses,
            charges = charges,
            c = c,
            dt = dt,
        )
        alg = if regularization_enabled
            RegularizedIntegrator(
                SymmetricProjectionIntegrator();
                backend = backend,
                warn_on_fallback = warn_on_fallback,
                r_on = r_on,
                r_off = r_off,
                max_substeps = max_substeps,
            )
        else
            SymmetricProjectionIntegrator()
        end
        return prob, alg
    end

    state_error(sol_a::HamiltonianSolution, sol_b::HamiltonianSolution) =
        norm(sol_a.q[end] - sol_b.q[end]) + norm(sol_a.p[end] - sol_b.p[end])

    @testset "API and fallback" begin
        sys2 = HamiltonianSystem(2, 2)
        q0_2d = [-1.0, 0.0, 1.0, 0.0]
        p0_2d = [0.0, -0.05, 0.0, 0.05]

        alg_valid = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
        )
        @test alg_valid.options.backend == WeberElectrodynamics.REG_BACKEND_ADAPTIVE

        @test_throws AssertionError RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :invalid_backend,
        )

        sys1 = HamiltonianSystem(2, 1)
        q0_1d = [-0.08, 0.08]
        p0_1d = [0.0, 0.0]

        prob_fb = HamiltonianProblem(
            sys1,
            (0.0, 0.01),
            q0_1d,
            p0_1d;
            masses = [1.0, 1.0],
            charges = [0.2, -0.2],
            c = 4.0,
            dt = 0.001,
        )
        alg_fb = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = false,
            r_on = 0.2,
            r_off = 0.26,
        )
        int_fb = init(prob_fb, alg_fb)
        rb_fb = int_fb.buffers.regularization_buffers
        @test rb_fb.effective_backend == WeberElectrodynamics.REG_BACKEND_ADAPTIVE
        @test rb_fb.backend_fallback

        step!(int_fb)
        @test int_fb.diagnostics.adaptive_pair_steps == 1
        @test int_fb.diagnostics.lifted_pair_steps == 0
        @test int_fb.diagnostics.backend_fallback_steps == 1

        alg_warn = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = true,
            r_on = 0.2,
            r_off = 0.26,
        )
        @test_logs (:warn, r"falling back to :adaptive_cartesian") init(prob_fb, alg_warn)
    end

    @testset "3D lifted pair fallback to adaptive_cartesian" begin
        sys3 = HamiltonianSystem(2, 3)
        # Start within r_on (separation 0.16 < r_on 0.2) so regularization fires
        # on the first step and backend_fallback_steps is immediately exercised.
        q0_3d = [-0.08, 0.0, 0.0, 0.08, 0.0, 0.0]
        p0_3d = [0.0, -0.05, 0.0, 0.0, 0.05, 0.0]

        prob_3d = HamiltonianProblem(
            sys3,
            (0.0, 0.01),
            q0_3d,
            p0_3d;
            masses = [1.0, 0.5],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
        )
        alg_3d = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = false,
            r_on = 0.2,
            r_off = 0.3,
        )
        int_3d = init(prob_3d, alg_3d)
        rb_3d = int_3d.buffers.regularization_buffers
        @test rb_3d.effective_backend == WeberElectrodynamics.REG_BACKEND_ADAPTIVE
        @test rb_3d.backend_fallback

        step!(int_3d)
        @test int_3d.diagnostics.backend_fallback_steps == 1
        @test int_3d.diagnostics.lifted_pair_steps == 0

        alg_3d_warn = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = true,
            r_on = 0.2,
            r_off = 0.3,
        )
        @test_logs (:warn, r"falling back to :adaptive_cartesian") init(prob_3d, alg_3d_warn)
    end

    @testset "Transform identities" begin
        @testset "Levi-Civita" begin
            Random.seed!(42)
            rb = WeberElectrodynamics.RegularizationBuffers(2, 2, 4, 0.1, 0.2, :lifted_pair, false, RegularizationOptions())
            q_rel = zeros(2)
            p_rel = zeros(2)

            for _ = 1:32
                u1 = randn()
                u2 = randn()
                x = u1^2 - u2^2
                y = 2u1 * u2
                r = sqrt(x^2 + y^2)
                rho2 = u1^2 + u2^2
                @test r ≈ rho2 rtol = 1e-12 atol = 1e-12

                J = [2u1 -2u2; 2u2 2u1]
                @test J * transpose(J) ≈ 4rho2 * I(2) rtol = 1e-12 atol = 1e-12

                rb.rel_q[1] = x
                rb.rel_q[2] = y
                rb.rel_p[1] = randn()
                rb.rel_p[2] = randn()

                WeberElectrodynamics._lc_lift!(rb)
                WeberElectrodynamics._lc_project!(q_rel, p_rel, rb.lc_u, rb.lc_U)

                @test q_rel[1] ≈ x rtol = 1e-10 atol = 1e-10
                @test q_rel[2] ≈ y rtol = 1e-10 atol = 1e-10
            end
        end

        @testset "KS" begin
            rb = WeberElectrodynamics.RegularizationBuffers(2, 3, 6, 0.1, 0.2, :adaptive_cartesian, false, RegularizationOptions())
            J = rb.ks_J

            for _ = 1:32
                u = randn(4)
                x1 = u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2
                x2 = 2 * (u[1] * u[2] - u[3] * u[4])
                x3 = 2 * (u[1] * u[3] + u[2] * u[4])
                r = sqrt(x1^2 + x2^2 + x3^2)
                rho2 = sum(abs2, u)
                @test r ≈ rho2 rtol = 1e-12 atol = 1e-12

                rb.ks_u .= u
                WeberElectrodynamics._ks_jacobian!(J, rb.ks_u)
                @test J * transpose(J) ≈ 4rho2 * I(3) rtol = 1e-12 atol = 1e-12

                p = randn(3)
                rb.ks_U[1] = J[1, 1] * p[1] + J[2, 1] * p[2] + J[3, 1] * p[3]
                rb.ks_U[2] = J[1, 2] * p[1] + J[2, 2] * p[2] + J[3, 2] * p[3]
                rb.ks_U[3] = J[1, 3] * p[1] + J[2, 3] * p[2] + J[3, 3] * p[3]
                rb.ks_U[4] = J[1, 4] * p[1] + J[2, 4] * p[2] + J[3, 4] * p[3]

                psi = WeberElectrodynamics._ks_constraint(rb.ks_u, rb.ks_U)
                @test abs(psi) < 1e-10

                rb.ks_U .+= randn(4) * 1e-4
                c_err = WeberElectrodynamics._ks_project_constraint!(rb.ks_U, rb.ks_u, rb.ks_n)
                @test c_err < 1e-10
            end
        end
    end

    @testset "Switching and hysteresis" begin
        rb = WeberElectrodynamics.RegularizationBuffers(3, 2, 6, 0.2, 0.3, :adaptive_cartesian, false, RegularizationOptions())

        q_activate = [-0.09, 0.0, 0.09, 0.0, 0.8, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_activate, true)
        @test active
        @test mode == WeberElectrodynamics.REG_MODE_PAIR

        q_between = [-0.13, 0.0, 0.13, 0.0, 0.8, 0.0]
        for _ = 1:8
            active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_between, true)
            @test active
            @test mode == WeberElectrodynamics.REG_MODE_PAIR
        end

        q_off = [-0.17, 0.0, 0.17, 0.0, 0.8, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_off, true)
        @test !active
        @test mode == WeberElectrodynamics.REG_MODE_NONE

        q_pair = [-0.09, 0.0, 0.09, 0.0, 1.0, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_pair, true)
        @test active
        @test mode == WeberElectrodynamics.REG_MODE_PAIR

        q_chain = [-0.09, 0.0, 0.09, 0.0, 0.2, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_chain, true)
        @test active
        @test mode == WeberElectrodynamics.REG_MODE_CHAIN

        q_pair_again = [-0.09, 0.0, 0.09, 0.0, 0.6, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_pair_again, true)
        @test active
        @test mode == WeberElectrodynamics.REG_MODE_PAIR
    end

    @testset "Chain-disabled overlap fallback" begin
        rb = WeberElectrodynamics.RegularizationBuffers(3, 2, 6, 0.2, 0.3, :adaptive_cartesian, false, RegularizationOptions())
        q_overlap = [-0.09, 0.0, 0.09, 0.0, 0.2, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_overlap, false)
        @test active
        @test rb.active_count > 2
        @test mode == WeberElectrodynamics.REG_MODE_NONE

        sys = HamiltonianSystem(3, 2)
        q0 = [-0.09, 0.0, 0.09, 0.0, 0.2, 0.0]
        p0 = [0.0, 0.02, 0.0, 0.0, 0.0, -0.02]
        prob = HamiltonianProblem(
            sys,
            (0.0, 0.1),
            q0,
            p0;
            masses = [1.0, 1.0, 1.0],
            charges = [0.1, -0.2, 0.1],
            c = 4.0,
            dt = 0.001,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            r_on = 0.2,
            r_off = 0.3,
            chain_enabled = false,
        )
        sol = solve(prob, alg)
        @test sol.retcode == :Success
        @test sol.regularization.pair_steps == 0
        @test sol.regularization.chain_steps == 0
        @test sol.regularization.unregularized_steps == length(sol) - 1
    end

    @testset "Pair mode correctness (2D lifted)" begin
        prob_cart, alg_cart = make_orbit_problem(2;
            dt = 0.004,
            t_end = 3.0,
            v_scale = 0.2,
            regularization_enabled = false,
            backend = :lifted_pair,
            r_on = 0.6,
            r_off = 0.9,
        )
        prob_adaptive, alg_adaptive = make_orbit_problem(2;
            dt = 0.004,
            t_end = 3.0,
            v_scale = 0.2,
            regularization_enabled = true,
            backend = :adaptive_cartesian,
            r_on = 0.6,
            r_off = 0.9,
        )
        prob_lifted, alg_lifted = make_orbit_problem(2;
            dt = 0.004,
            t_end = 3.0,
            v_scale = 0.2,
            regularization_enabled = true,
            backend = :lifted_pair,
            r_on = 0.6,
            r_off = 0.9,
        )
        prob_ref, alg_ref = make_orbit_problem(2;
            dt = 0.001,
            t_end = 3.0,
            v_scale = 0.2,
            regularization_enabled = true,
            backend = :lifted_pair,
            r_on = 0.6,
            r_off = 0.9,
        )

        sol_cart = solve(prob_cart, alg_cart)
        sol_adaptive = solve(prob_adaptive, alg_adaptive)
        sol_lifted = solve(prob_lifted, alg_lifted)
        sol_ref = solve(prob_ref, alg_ref)

        @test sol_cart.retcode == :Success
        @test sol_adaptive.retcode == :Success
        @test sol_lifted.retcode == :Success
        @test sol_ref.retcode == :Success

        @test all(isfinite, sol_lifted.q[end])
        @test all(isfinite, sol_lifted.p[end])
        @test sol_lifted.regularization.pair_steps > 0
        @test sol_lifted.regularization.lifted_pair_steps > 0
        @test sol_adaptive.regularization.adaptive_pair_steps > 0

        err_cart = state_error(sol_cart, sol_ref)
        err_lifted = state_error(sol_lifted, sol_ref)
        err_adaptive = state_error(sol_adaptive, sol_ref)

        @test err_lifted <= 0.5 * err_cart
        @test err_adaptive <= err_cart
    end

    @testset "Sub-critical like-charge oscillation (collision bounce)" begin
        # Two equal positive charges inside the critical radius ρ.
        # With ℓ=0 (head-on), they attract and oscillate between r₀ and r≈0.
        # The collision at r=0 is C⁰-continuable (Frauenfelder & Weber 2024).
        # ρ = q₁q₂/(μc²) = 1/(0.5·16) = 0.125
        m1 = m2 = 1.0
        q1 = q2 = 1.0
        c = 4.0
        r0 = 0.05     # < ρ = 0.125
        mu = m1 * m2 / (m1 + m2)
        rho = q1 * q2 / (mu * c^2)
        T_est = (2π / c) * rho

        sys = HamiltonianSystem(2, 2)
        q_init = [r0 / 2, 0.0, -r0 / 2, 0.0]
        p_init = [0.0, 0.0, 0.0, 0.0]

        prob = HamiltonianProblem(
            sys,
            (0.0, 10 * T_est),
            q_init,
            p_init;
            masses = [m1, m2],
            charges = [q1, q2],
            c = c,
            dt = T_est / 2000,
        )

        sol = solve(
            prob,
            SymmetricProjectionIntegrator();
            callbacks = CollisionBounce(0.01),
        )
        @test sol.retcode == :Success

        # All states finite.
        @test all(isfinite, sol.q[end])
        @test all(isfinite, sol.p[end])

        # Particles stay bounded: r ≤ r₀ (they never exceed the starting separation).
        rs = [sqrt((sol.q[k][1] - sol.q[k][3])^2 + (sol.q[k][2] - sol.q[k][4])^2)
              for k in 1:length(sol.t)]
        @test maximum(rs) <= r0 + 1e-6

        # Energy conservation < 1%.
        en = compute_energy_timeseries(sol)
        @test en.statistics.global_error_percent_max < 1.0

        # Multiple oscillation cycles occur (at least 10 half-periods in 10T).
        n_minima = count(k -> rs[k] < rs[k-1] && rs[k] < rs[k+1], 2:length(rs)-1)
        @test n_minima >= 10
    end

    @testset "Collision bounce with regularization enabled" begin
        # Smoke test: collision bounce + adaptive_cartesian regularization simultaneously.
        # Documented as "works best without LC" but must not crash with adaptive Cartesian.
        m1 = m2 = 1.0
        q1 = q2 = 1.0
        c = 4.0
        r0 = 0.05
        mu = m1 * m2 / (m1 + m2)
        rho = q1 * q2 / (mu * c^2)
        T_est = (2π / c) * rho

        sys = HamiltonianSystem(2, 2)
        q_init = [r0 / 2, 0.0, -r0 / 2, 0.0]
        p_init = [0.0, 0.0, 0.0, 0.0]

        prob = HamiltonianProblem(
            sys,
            (0.0, 5 * T_est),
            q_init,
            p_init;
            masses = [m1, m2],
            charges = [q1, q2],
            c = c,
            dt = T_est / 1000,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            # Set r_on below bounce_r so regularization stays idle;
            # bounce handles the singularity and the two features coexist.
            r_on = 0.005,
            r_off = 0.015,
            collision_bounce_radius = 0.01,
        )

        sol = solve(prob, alg)
        @test sol.retcode == :Success
        @test all(isfinite, sol.q[end])
        @test all(isfinite, sol.p[end])

        en = compute_energy_timeseries(sol)
        @test en.statistics.global_error_percent_max < 1.0
    end

    @testset "Chain mode correctness" begin
        sys = HamiltonianSystem(3, 2)
        q0 = [-0.12, 0.0, 0.0, 0.0, 0.12, 0.0]
        p0 = [0.0, 0.02, 0.0, 0.0, 0.0, -0.02]

        prob = HamiltonianProblem(
            sys,
            (0.0, 0.2),
            q0,
            p0;
            masses = [1.0, 1.0, 1.0],
            charges = [0.1, -0.2, 0.1],
            c = 4.0,
            dt = 0.002,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = false,
            r_on = 0.28,
            r_off = 0.39,
            chain_enabled = true,
            max_substeps = 256,
        )

        sol = solve(prob, alg)
        @test sol.retcode == :Success
        @test sol.regularization.chain_steps > 0
        @test sol.regularization.lifted_pair_steps == 0
    end

    @testset "Backward mode parity" begin
        sys = HamiltonianSystem(2, 2)
        q0 = [-1.0, 0.0, 1.0, 0.0]
        p0 = [0.0, -0.05, 0.0, 0.05]

        kwargs = (
            masses = [1.0, 0.5],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
        )

        prob_shared = HamiltonianProblem(
            sys,
            (0.0, 0.2),
            q0,
            p0;
            kwargs...,
        )
        alg_on = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            r_on = 0.05,
            r_off = 0.08,
        )

        sol_off = solve(prob_shared, SymmetricProjectionIntegrator())
        sol_on = solve(prob_shared, alg_on)

        @test sol_on.regularization.pair_steps == 0
        @test sol_on.regularization.chain_steps == 0
        @test sol_on.regularization.unregularized_steps == length(sol_on) - 1
        @test sol_on.regularization.used_backend == WeberElectrodynamics.REG_BACKEND_DISABLED
        @test sol_off.t == sol_on.t
        @test all(sol_off.q .== sol_on.q)
        @test all(sol_off.p .== sol_on.p)
    end

    @testset "Single-particle regularization safety" begin
        sys = HamiltonianSystem(1, 2)
        q0 = [0.0, 0.0]
        p0 = [0.0, 0.0]
        prob = HamiltonianProblem(
            sys,
            (0.0, 0.1),
            q0,
            p0;
            masses = [1.0],
            charges = [0.0],
            c = 4.0,
            dt = 0.01,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
        )
        sol = solve(prob, alg)
        @test sol.retcode == :Success
        @test sol.regularization.pair_steps == 0
        @test sol.regularization.chain_steps == 0
        @test sol.regularization.unregularized_steps == length(sol) - 1
    end

    @testset "Diagnostics validity" begin
        prob, alg = make_orbit_problem(2;
            dt = 0.004,
            t_end = 3.0,
            v_scale = 0.2,
            regularization_enabled = true,
            backend = :lifted_pair,
            r_on = 0.6,
            r_off = 0.9,
        )
        sol = solve(prob, alg)
        d = sol.regularization

        @test d.pair_steps + d.chain_steps + d.unregularized_steps == length(sol) - 1
        @test d.pair_steps == d.adaptive_pair_steps + d.lifted_pair_steps
        @test d.total_substeps >= d.pair_steps + d.chain_steps + d.unregularized_steps
        @test d.max_substeps_used >= 1
        @test isfinite(d.min_encounter_distance)
        @test d.max_constraint_violation >= 0
        @test d.requested_backend == WeberElectrodynamics.REG_BACKEND_LIFTED
        @test d.used_backend == WeberElectrodynamics.REG_BACKEND_LIFTED
    end

    @testset "Allocation checks" begin
        sys = HamiltonianSystem(2, 2)

        q0 = [-1.0, 0.0, 1.0, 0.0]
        p0 = [0.0, -0.05, 0.0, 0.05]
        prob_unreg = HamiltonianProblem(
            sys,
            (0.0, 0.02),
            q0,
            p0;
            masses = [1.0, 0.5],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
        )
        int_unreg = init(prob_unreg)
        step!(int_unreg)
        alloc_unreg = @allocated step!(int_unreg)
        @test alloc_unreg <= 512

        q1 = [-0.08, 0.01, 0.08, -0.01]
        p1 = [0.0, -0.02, 0.0, 0.02]
        prob_lifted = HamiltonianProblem(
            sys,
            (0.0, 0.02),
            q1,
            p1;
            masses = [1.0, 1.0],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
        )
        alg_lifted = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = false,
            r_on = 0.2,
            r_off = 0.26,
            max_substeps = 256,
        )
        int_lifted = init(prob_lifted, alg_lifted)
        step!(int_lifted)
        @test int_lifted.diagnostics.mode_history[1] == WeberElectrodynamics.REG_MODE_PAIR
        alloc_lifted = @allocated step!(int_lifted)
        @test alloc_lifted <= 2048
    end
end
