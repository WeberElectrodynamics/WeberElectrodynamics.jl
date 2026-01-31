@testset "Integration Tests" begin
    @testset "Full two-body workflow" begin
        # Setup Weber system
        m1, m2 = 1.0, 0.1
        q1_charge = sqrt(0.1)
        q2_charge = -sqrt(0.1)  # q1*q2 = -0.1 (attractive)
        c = 4.0

        system = WeberSystem(2, 2)

        r0 = 2.0
        M = m1 + m2
        k = q1_charge * q2_charge
        v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))
        q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
        p0 = [0.0, m1 * (-m2 / M * v_circ * 0.9), 0.0, m2 * (m1 / M * v_circ * 0.9)]

        prob = WeberProblem(system, (0.0, 1.0), q0, p0;
            masses = [m1, m2], charges = [q1_charge, q2_charge], c = c, dt = 0.001)

        # Solve
        sol = solve(prob)
        @test sol.retcode == :Success

        # Trajectories
        traj = compute_trajectory_data(sol, 2, 2; stride = 10)
        @test traj.n_particles == 2

        # Energy
        function total_energy(q, p, params, t)
            weber_energy_2body_2d(q, p, [m1, m2], [q1_charge, q2_charge], c)
        end
        energy = compute_energy_timeseries(
            sol,
            total_energy,
            nothing,
            nothing,
            nothing;
            stride = 10,
        )
        @test energy.relative_energy_range < 1e-6

        # Forces
        forces = compute_force_timeseries(
            sol,
            2,
            2,
            [m1, m2],
            [q1_charge, q2_charge],
            c;
            stride = 10,
        )
        n3 = check_newtons_third_law(forces)
        @test n3.global_max_violation < 1e-6

        # Phase space
        ps = compute_phase_space_data(sol, 2, 2, [m1, m2]; stride = 10)
        @test all(ps.separation_distance .> 0)  # Particles don't collide
    end

    @testset "Stepped integration matches full solve" begin
        prob = make_weber_problem(tspan = (0.0, 0.5))

        # Full solve
        sol_full = solve(prob)

        # Stepped solve
        integrator = init(prob)
        while step!(integrator)
        end
        sol_stepped = solve!(integrator)

        # Should match
        @test length(sol_full.t) == length(sol_stepped.t)
        @test sol_full.t ≈ sol_stepped.t
    end

    @testset "Reproducibility" begin
        prob = make_weber_problem(tspan = (0.0, 0.5))

        sol1 = solve(prob)
        sol2 = solve(prob)

        @test sol1.t == sol2.t
        @test all(sol1.q .== sol2.q)
        @test all(sol1.p .== sol2.p)
    end

    @testset "Different initial conditions" begin
        # Test with different orbital configurations
        m1, m2 = 1.0, 0.1
        q1_charge, q2_charge = 0.1, -0.1
        c = 4.0

        for scale in [0.5, 1.0, 2.0]
            system = WeberSystem(2, 2)
            r0 = 2.0 * scale
            M = m1 + m2
            k = q1_charge * q2_charge
            v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))
            q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
            p0 = [0.0, m1 * (-m2 / M * v_circ * 0.9), 0.0, m2 * (m1 / M * v_circ * 0.9)]

            prob = WeberProblem(system, (0.0, 1.0), q0, p0;
                masses = [m1, m2], charges = [q1_charge, q2_charge], c = c, dt = 0.001)
            sol = solve(prob)
            @test sol.retcode == :Success
        end
    end

    @testset "Long-time stability" begin
        # Run for many orbits with Coulomb-like (large c)
        prob = make_coulomb_like_problem(tspan = (0.0, 10.0), dt = 0.01)
        sol = solve(prob)

        @test sol.retcode == :Success

        # Energy should be conserved
        masses = prob.masses
        charges = prob.charges
        E(q, p) = coulomb_like_energy_2body_2d(q, p, masses, charges)
        E0 = E(sol.q[1], sol.p[1])
        E_final = E(sol.q[end], sol.p[end])
        @test abs(E_final - E0) / abs(E0) < 1e-6
    end

    @testset "API consistency" begin
        # Verify all exported functions work together
        prob = make_coulomb_like_problem(tspan = (0.0, 0.5), dt = 0.01)

        # Method 1: solve directly
        sol1 = solve(prob)

        # Method 2: init + solve!
        int2 = init(prob)
        sol2 = solve!(int2)

        # Method 3: init + step! loop + solve!
        int3 = init(prob)
        for _ = 1:5
            step!(int3)
        end
        sol3 = solve!(int3)

        # All should produce valid solutions
        @test sol1.retcode == :Success
        @test sol2.retcode == :Success
        @test sol3.retcode == :Success

        # sol1 and sol2 should be identical
        @test sol1.t == sol2.t
        @test all(sol1.q .== sol2.q)
    end

    @testset "Statistics on same solution" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
        sol = solve(prob)
        masses = prob.masses
        charges = prob.charges

        # All statistics should work on the same solution
        traj = compute_trajectory_data(sol, 2, 2)
        @test traj isa TrajectoryData

        energy_func(q, p, params, t) = coulomb_like_energy_2body_2d(q, p, masses, charges)
        energy = compute_energy_timeseries(sol, energy_func, nothing, nothing, nothing)
        @test energy isa EnergyData

        forces = compute_force_timeseries(sol, 2, 2, masses, charges, 1e10)
        @test forces isa ForceData

        n3 = check_newtons_third_law(forces)
        @test n3 isa NewtonsThirdLawData

        ps = compute_phase_space_data(sol, 2, 2, masses)
        @test ps isa PhaseSpaceData

        # All should have consistent time lengths
        @test length(traj.trajectories[1][:, 1]) == length(sol.t)
        @test length(energy.t) == length(sol.t)
        @test length(ps.t) == length(sol.t)
    end
end
