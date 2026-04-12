using Plots

@testset "Plots extension" begin
    @testset "2D smoke tests (Coulomb-like)" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
        sol = solve(prob)

        traj = compute_trajectory_data(sol, 2, 2)
        energy = compute_energy_timeseries(sol)
        momentum = compute_momentum_timeseries(sol)
        forces = compute_pair_force_timeseries(
            sol, (1, 2), 2, 2, [1.0, 0.5], [1.0, -1.0], 1e10,
        )

        @test plot_trajectories(traj) isa Plots.Plot
        @test plot_energy(energy) isa Plots.Plot
        @test plot_pair_energy(energy, (1, 2)) isa Plots.Plot
        @test plot_energy_errors(energy) isa Plots.Plot
        @test plot_pair_forces(forces) isa Plots.Plot
        @test plot_phase_space(forces) isa Plots.Plot
        @test plot_momentum(momentum) isa Plots.Plot

        # Validation: unknown pair should throw
        @test_throws ArgumentError plot_pair_energy(energy, (1, 5))
    end

    @testset "1D trajectories" begin
        sys1d = WeberSystem(2, 1)
        prob1d = WeberProblem(
            sys1d, (0.0, 0.5),
            [1.0, -1.0], [0.1, -0.1];
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 1e10, dt = 0.01,
        )
        sol1d = solve(prob1d)
        traj1d = compute_trajectory_data(sol1d, 2, 1)
        @test plot_trajectories(traj1d) isa Plots.Plot
    end

    @testset "3D trajectories" begin
        sys3d = WeberSystem(2, 3)
        prob3d = WeberProblem(
            sys3d, (0.0, 0.5),
            [1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
            [0.0, 0.1, 0.0, 0.0, -0.1, 0.0];
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 1e10, dt = 0.01,
        )
        sol3d = solve(prob3d)
        traj3d = compute_trajectory_data(sol3d, 2, 3)
        @test plot_trajectories(traj3d) isa Plots.Plot
    end

    @testset "Zöllner plots" begin
        sys = WeberSystem(2, 2)
        prob_z = WeberProblem(
            sys, (0.0, 1.0),
            [-1.0, 0.0, 1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 10.0, dt = 0.01,
            zollner = ZollnerOptions(enabled = true, a = 0.05),
        )
        prob_w = WeberProblem(
            sys, (0.0, 1.0),
            [-1.0, 0.0, 1.0, 0.0], [0.0, 0.3, 0.0, -0.3];
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 10.0, dt = 0.01,
        )

        sol_z = solve(prob_z)
        sol_w = solve(prob_w)

        energy_z = compute_energy_timeseries(sol_z)
        forces_z = compute_pair_force_timeseries(
            sol_z, (1, 2), 2, 2, [1.0, 1.0], [1.0, -1.0], 10.0,
        )
        forces_w = compute_pair_force_timeseries(
            sol_w, (1, 2), 2, 2, [1.0, 1.0], [1.0, -1.0], 10.0,
        )

        @test plot_zollner_energy(energy_z) isa Plots.Plot
        @test plot_zollner_force_residual(forces_z) isa Plots.Plot
        @test plot_weber_vs_zollner(sol_w, sol_z) isa Plots.Plot
        @test plot_weber_vs_zollner(sol_w, sol_z; labels = ["A", "B"]) isa Plots.Plot
        @test plot_zollner_phase_space(forces_w, forces_z) isa Plots.Plot
    end

    @testset "plot_weber_vs_zollner 3D branch" begin
        sys3d = WeberSystem(2, 3)
        q0 = [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        p0 = [0.0, 0.2, 0.0, 0.0, -0.2, 0.0]
        prob_a = WeberProblem(
            sys3d, (0.0, 0.5), q0, p0;
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 10.0, dt = 0.01,
        )
        prob_b = WeberProblem(
            sys3d, (0.0, 0.5), q0, p0;
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 10.0, dt = 0.01,
            zollner = ZollnerOptions(enabled = true, a = 0.02),
        )
        sol_a = solve(prob_a)
        sol_b = solve(prob_b)
        @test plot_weber_vs_zollner(sol_a, sol_b) isa Plots.Plot
    end
end
