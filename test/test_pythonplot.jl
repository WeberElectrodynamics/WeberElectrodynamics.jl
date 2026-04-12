@testset "PythonPlot extension (optional)" begin
    can_load = try
        @eval using PythonPlot
        true
    catch err
        @info "PythonPlot not available — skipping PythonPlot extension tests" exception = err
        false
    end
    if !can_load
        return
    end

    prob = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
    sol = solve(prob)

    traj = compute_trajectory_data(sol, 2, 2)
    energy = compute_energy_timeseries(sol)
    forces = compute_pair_force_timeseries(
        sol, (1, 2), 2, 2,
        sol.prob.masses, sol.prob.charges, sol.prob.c,
    )
    momentum = compute_momentum_timeseries(sol)

    # Non-interactive backend to avoid opening windows on headless CI/local.
    PythonPlot.matplotlib.use("Agg")

    @testset "figure creation smoke tests" begin
        @test plot_trajectories(traj)            !== nothing
        @test plot_energy(energy)                !== nothing
        @test plot_pair_energy(energy, (1, 2))   !== nothing
        @test plot_energy_errors(energy)         !== nothing
        @test plot_pair_forces(forces)           !== nothing
        @test plot_phase_space(forces)           !== nothing
        @test plot_momentum(momentum)            !== nothing

        # Zöllner plots also run on a standard (non-Zöllner) problem — residuals are zero.
        @test plot_zollner_energy(energy)                   !== nothing
        @test plot_zollner_force_residual(forces)           !== nothing
        @test plot_weber_vs_zollner(sol, sol)               !== nothing
        @test plot_zollner_phase_space(forces, forces)      !== nothing
    end

    PythonPlot.close("all")
end
