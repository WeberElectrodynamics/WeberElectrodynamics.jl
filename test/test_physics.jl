@testset "Physics Validation" begin
    @testset "Energy conservation - Coulomb-like" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 5.0), dt = 0.001)
        sol = solve(prob)
        ms = masses(prob)
        qs = charges(prob)

        E(q, p) = coulomb_like_energy_2body_2d(q, p, ms, qs)
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i = 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-8
    end

    @testset "Energy conservation - Weber" begin
        prob = make_weber_problem(tspan = (0.0, 2.0), dt = 0.0005)
        sol = solve(prob)
        ms = masses(prob)
        qs = charges(prob)
        c_val = speed_of_light(prob)

        E(q, p) = weber_energy_2body_2d(q, p, ms, qs, c_val)
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i = 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-6
    end

    @testset "Weber force decomposition - Coulomb-like" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
        sol = solve(prob)

        forces = compute_pair_force_timeseries(
            sol,
            (1, 2),
            2,
            2,
            masses(prob),
            charges(prob),
            speed_of_light(prob),
        )

        # For large c (Coulomb limit), Weber correction terms should be small
        # Vector form: F ≈ Coulomb (v·v, r·a, rv² terms should be negligible)
        for t = 1:length(forces.t)
            coulomb_mag = norm(forces.coulomb[t])
            vv_mag = norm(forces.vector_term_vv[t])
            ra_mag = norm(forces.vector_term_ra[t])
            rv2_mag = norm(forces.vector_term_rv2[t])

            # Weber corrections should be tiny relative to Coulomb
            @test vv_mag / coulomb_mag < 1e-10
            @test ra_mag / coulomb_mag < 1e-10
            @test rv2_mag / coulomb_mag < 1e-10
        end
    end

    @testset "Weber force decomposition - consistency" begin
        prob = make_weber_problem(tspan = (0.0, 1.0), dt = 0.001)
        sol = solve(prob)

        forces = compute_pair_force_timeseries(
            sol,
            (1, 2),
            2,
            2,
            masses(prob),
            charges(prob),
            speed_of_light(prob),
        )

        # Vector form and radial form should give same total force
        for t = 1:length(forces.t)
            force_vector =
                forces.coulomb[t] .+ forces.vector_term_vv[t] .+ forces.vector_term_ra[t] .+
                forces.vector_term_rv2[t]
            force_radial =
                forces.coulomb[t] .+ forces.radial_term_rdot2[t] .+
                forces.radial_term_rddot[t]

            @test forces.force[t] ≈ force_vector rtol = 1e-12
            @test force_vector ≈ force_radial rtol = 1e-10
        end
    end

    @testset "Center of mass conservation" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 2.0), dt = 0.01)
        sol = solve(prob)

        m1, m2 = masses(prob)
        M = m1 + m2

        # Initial center of mass position and momentum
        x_cm_0 = (m1 * sol.q[1][1] + m2 * sol.q[1][3]) / M
        y_cm_0 = (m1 * sol.q[1][2] + m2 * sol.q[1][4]) / M
        px_cm_0 = sol.p[1][1] + sol.p[1][3]
        py_cm_0 = sol.p[1][2] + sol.p[1][4]

        # Check conservation at all times
        for i = 1:length(sol)
            x_cm = (m1 * sol.q[i][1] + m2 * sol.q[i][3]) / M
            y_cm = (m1 * sol.q[i][2] + m2 * sol.q[i][4]) / M
            px_cm = sol.p[i][1] + sol.p[i][3]
            py_cm = sol.p[i][2] + sol.p[i][4]

            # Position drifts with velocity
            expected_x = x_cm_0 + px_cm_0 / M * sol.t[i]
            expected_y = y_cm_0 + py_cm_0 / M * sol.t[i]

            @test x_cm ≈ expected_x rtol = 1e-8
            @test y_cm ≈ expected_y rtol = 1e-8
            @test px_cm ≈ px_cm_0 rtol = 1e-10
            @test py_cm ≈ py_cm_0 rtol = 1e-10
        end
    end

    @testset "Angular momentum conservation (Coulomb-like)" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 2.0), dt = 0.01)
        sol = solve(prob)

        function angular_momentum(q, p)
            # L = r1 × p1 + r2 × p2 (z-component in 2D)
            q[1] * p[2] - q[2] * p[1] + q[3] * p[4] - q[4] * p[3]
        end

        L0 = angular_momentum(sol.q[1], sol.p[1])

        for i = 1:length(sol)
            Li = angular_momentum(sol.q[i], sol.p[i])
            @test Li ≈ L0 rtol = 1e-8
        end
    end

    @testset "Symplectic integrator: no secular drift" begin
        # Long simulation to check for secular drift
        prob = make_coulomb_like_problem(tspan = (0.0, 50.0), dt = 0.01)
        sol = solve(prob)

        # Energy at beginning, middle, and end
        ms = masses(prob)
        qs = charges(prob)
        E(q, p) = coulomb_like_energy_2body_2d(q, p, ms, qs)
        E_start = E(sol.q[1], sol.p[1])
        E_mid = E(sol.q[length(sol)÷2], sol.p[length(sol)÷2])
        E_end = E(sol.q[end], sol.p[end])

        # No secular drift: errors should be bounded, not growing
        error_mid = abs(E_mid - E_start)
        error_end = abs(E_end - E_start)

        # End error should not be significantly larger than mid error
        # (would indicate secular drift)
        @test error_end < 2 * error_mid + 1e-10
    end

    @testset "Phase space orbit closure (Weber)" begin
        # For a near-circular bound orbit, phase space should approximately close
        # Use Weber system
        prob = make_weber_problem(tspan = (0.0, 20.0), dt = 0.001)
        sol = solve(prob)

        # For a bound orbit, the particle should stay bounded
        r_min = Inf
        r_max = 0.0
        for i = 1:length(sol)
            dx = sol.q[i][1] - sol.q[i][3]
            dy = sol.q[i][2] - sol.q[i][4]
            r = sqrt(dx^2 + dy^2)
            r_min = min(r_min, r)
            r_max = max(r_max, r)
        end

        # Orbit should be bounded (not escaping)
        # Weber orbits precess and may have larger excursions than Kepler
        @test r_max < 25.0
        @test r_min > 0.05  # Don't get too close (no collision)
    end

    @testset "Bound vs unbound orbits" begin
        # Negative total energy = bound orbit
        prob_bound = make_coulomb_like_problem(tspan = (0.0, 2.0), dt = 0.01)
        sol_bound = solve(prob_bound)

        ms = masses(prob_bound)
        qs = charges(prob_bound)
        E0 = coulomb_like_energy_2body_2d(sol_bound.q[1], sol_bound.p[1], ms, qs)
        @test E0 < 0  # Bound orbit has negative energy

        # Separation should remain bounded
        r_max = 0.0
        for i = 1:length(sol_bound)
            dx = sol_bound.q[i][1] - sol_bound.q[i][3]
            dy = sol_bound.q[i][2] - sol_bound.q[i][4]
            r = sqrt(dx^2 + dy^2)
            r_max = max(r_max, r)
        end
        @test r_max < 10.0  # Should stay bounded
    end

    @testset "Weber vs Coulomb limit" begin
        # As c → ∞, Weber should approach Coulomb
        m1, m2 = 1.0, 0.5
        q1, q2 = 1.0, -1.0

        # Single system, different params
        sys = HamiltonianSystem(2, 2)

        # Test at specific point
        q = [1.0, 0.0, -1.0, 0.0]
        p = [0.0, 0.5, 0.0, -0.5]

        # Params layout: [m1, m2, q1, q2, c]
        params_large_c = [m1, m2, q1, q2, 1e10]
        params_small_c = [m1, m2, q1, q2, 10.0]

        out_large_c = zeros(4)
        out_small_c = zeros(4)

        sys.dp_dt_compiled(out_large_c, q, p, 0.0, params_large_c)
        sys.dp_dt_compiled(out_small_c, q, p, 0.0, params_small_c)

        # Force magnitude should be similar but with small Weber correction for finite c
        @test norm(out_large_c) ≈ norm(out_small_c) rtol = 0.1
    end
end
