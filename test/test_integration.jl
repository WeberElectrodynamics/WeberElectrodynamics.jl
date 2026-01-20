@testset "Integration Tests" begin
    @testset "Full two-body workflow" begin
        # Build Hamiltonian
        H = build_hamiltonian(weber_H, 2, 2; param_names=[:m1, :m2, :k, :c])

        # Setup problem
        m1, m2, k, c = 1.0, 0.1, -0.1, 4.0
        r0 = 2.0
        M = m1 + m2
        v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))
        q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
        p0 = [0.0, m1 * (-m2 / M * v_circ * 0.9), 0.0, m2 * (m1 / M * v_circ * 0.9)]

        prob = WeberProblem(H, (0.0, 1.0), q0, p0; params=[m1, m2, k, c], dt=0.001)

        # Solve
        sol = solve(prob)
        @test sol.retcode == :Success

        # Trajectories
        traj = create_trajectory_data(sol, 2, 2; stride=10)
        @test traj.n_particles == 2

        # Energy
        total_energy(q, p, params, t) = weber_H(q, p, params)
        energy = compute_energy_timeseries(sol, total_energy, nothing, nothing, [m1, m2, k, c]; stride=10)
        @test energy.relative_energy_range < 1e-6

        # Forces
        forces = compute_force_timeseries(sol, 2, 2, [m1, m2], [1.0, -1.0], c; stride=10)
        n3 = check_newtons_third_law(forces)
        @test n3.global_max_violation < 1e-6

        # Phase space
        ps = compute_phase_space_data(sol, 2, 2, [m1, m2]; stride=10)
        @test all(ps.r .> 0)  # Particles don't collide
    end

    @testset "Stepped integration matches full solve" begin
        prob = make_weber_problem(tspan=(0.0, 0.5))

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
        prob = make_weber_problem(tspan=(0.0, 0.5))

        sol1 = solve(prob)
        sol2 = solve(prob)

        @test sol1.t == sol2.t
        @test all(sol1.q .== sol2.q)
        @test all(sol1.p .== sol2.p)
    end

    @testset "Different initial conditions" begin
        H = build_hamiltonian(harmonic_oscillator_H, 1, 1; param_names=[:m, :k])

        # Different amplitudes
        for amp in [0.1, 1.0, 10.0]
            prob = WeberProblem(H, (0.0, 1.0), [amp], [0.0]; params=[1.0, 1.0], dt=0.01)
            sol = solve(prob)
            @test sol.retcode == :Success

            # Check amplitude is preserved (approximately)
            max_q = maximum(abs.(getindex.(sol.q, 1)))
            @test max_q ≈ amp rtol = 0.05
        end
    end

    @testset "Long-time stability" begin
        # Run for many orbits
        prob = make_coulomb_problem(tspan=(0.0, 10.0), dt=0.01)
        sol = solve(prob)

        @test sol.retcode == :Success

        # Energy should be conserved
        params = [1.0, 0.5, 1.0]
        E(q, p) = coulomb_H(q, p, params)
        E0 = E(sol.q[1], sol.p[1])
        E_final = E(sol.q[end], sol.p[end])
        @test abs(E_final - E0) / abs(E0) < 1e-6
    end

    @testset "API consistency" begin
        # Verify all exported functions work together
        prob = make_coulomb_problem(tspan=(0.0, 0.5), dt=0.01)

        # Method 1: solve directly
        sol1 = solve(prob)

        # Method 2: init + solve!
        int2 = init(prob)
        sol2 = solve!(int2)

        # Method 3: init + step! loop + solve!
        int3 = init(prob)
        for _ in 1:5
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
        prob = make_coulomb_problem(tspan=(0.0, 1.0), dt=0.01)
        sol = solve(prob)
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        params = [1.0, 0.5, 1.0]

        # All statistics should work on the same solution
        traj = create_trajectory_data(sol, 2, 2)
        @test traj isa TrajectoryData

        energy_func(q, p, params, t) = coulomb_H(q, p, params)
        energy = compute_energy_timeseries(sol, energy_func, nothing, nothing, params)
        @test energy isa EnergyData

        forces = compute_force_timeseries(sol, 2, 2, masses, charges, 1e6)
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
