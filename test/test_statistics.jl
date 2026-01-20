@testset "Statistics" begin
    # Setup: run simulations for testing
    prob_coulomb = make_coulomb_problem(tspan=(0.0, 1.0), dt=0.01)
    sol_coulomb = solve(prob_coulomb)

    prob_harmonic = make_harmonic_problem(tspan=(0.0, 2.0), dt=0.01)
    sol_harmonic = solve(prob_harmonic)

    @testset "TrajectoryData" begin
        @testset "create_trajectory_data basic" begin
            traj = create_trajectory_data(sol_coulomb, 2, 2)

            @test traj isa TrajectoryData
            @test traj.n_particles == 2
            @test traj.dims == 2
            @test length(traj.trajectories) == 2
            @test size(traj.trajectories[1]) == (length(sol_coulomb.t), 2)
            @test length(traj.initial_positions) == 2
            @test length(traj.final_positions) == 2
            @test length(traj.initial_positions[1]) == 2
        end

        @testset "create_trajectory_data with stride" begin
            traj = create_trajectory_data(sol_coulomb, 2, 2; stride=10)

            expected_points = length(1:10:length(sol_coulomb.t))
            @test size(traj.trajectories[1], 1) == expected_points
        end

        @testset "Initial and final positions" begin
            traj = create_trajectory_data(sol_coulomb, 2, 2)

            # Initial positions should match solution
            for p in 1:2
                for d in 1:2
                    idx = (p - 1) * 2 + d
                    @test traj.initial_positions[p][d] == sol_coulomb.q[1][idx]
                    @test traj.final_positions[p][d] == sol_coulomb.q[end][idx]
                end
            end
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError create_trajectory_data(sol_coulomb, 2, 2; stride=0)
            @test_throws ArgumentError create_trajectory_data(sol_coulomb, 2, 2; stride=-1)
            @test_throws ArgumentError create_trajectory_data(sol_coulomb, 3, 2)  # Wrong n_particles
            @test_throws ArgumentError create_trajectory_data(sol_coulomb, 2, 3)  # Wrong dims
        end

        @testset "1D system" begin
            traj = create_trajectory_data(sol_harmonic, 1, 1)
            @test traj.n_particles == 1
            @test traj.dims == 1
            @test size(traj.trajectories[1], 2) == 1
        end
    end

    @testset "EnergyData" begin
        params = [1.0, 1.0]  # m, k

        total_energy(q, p, params, t) = harmonic_oscillator_H(q, p, params)
        kinetic_energy(q, p, params, t) = sum(p .^ 2) / (2 * params[1])
        potential_energy(q, p, params, t) = params[2] * sum(q .^ 2) / 2

        @testset "compute_energy_timeseries basic" begin
            energy = compute_energy_timeseries(sol_harmonic, total_energy, nothing, nothing, params)

            @test energy isa EnergyData
            @test length(energy.t) == length(sol_harmonic.t)
            @test length(energy.total) == length(sol_harmonic.t)
            @test isnothing(energy.kinetic)
            @test isnothing(energy.potential)
        end

        @testset "compute_energy_timeseries with components" begin
            energy = compute_energy_timeseries(sol_harmonic, total_energy, kinetic_energy, potential_energy, params)

            @test !isnothing(energy.kinetic)
            @test !isnothing(energy.potential)
            @test length(energy.kinetic) == length(sol_harmonic.t)
            @test length(energy.potential) == length(sol_harmonic.t)

            # KE + PE = Total (approximately)
            for i in 1:length(energy.t)
                @test energy.kinetic[i] + energy.potential[i] ≈ energy.total[i] rtol = 1e-10
            end
        end

        @testset "Energy statistics" begin
            energy = compute_energy_timeseries(sol_harmonic, total_energy, nothing, nothing, params)

            @test energy.max_local_error >= 0
            # For symplectic integrator, relative range should be small
            if !isnan(energy.relative_energy_range)
                @test energy.relative_energy_range < 1e-6
            end
        end

        @testset "compute_energy_timeseries with stride" begin
            energy = compute_energy_timeseries(sol_harmonic, total_energy, nothing, nothing, params; stride=10)

            expected_points = length(1:10:length(sol_harmonic.t))
            @test length(energy.t) == expected_points
            @test length(energy.total) == expected_points
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_energy_timeseries(sol_harmonic, total_energy, nothing, nothing, params; stride=0)
            @test_throws ArgumentError compute_energy_timeseries(sol_harmonic, total_energy, nothing, nothing, params; stride=-1)
        end
    end

    @testset "ForceData" begin
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        c = 1e6  # Large c for Coulomb limit

        @testset "compute_force_timeseries basic" begin
            forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c)

            @test forces isa ForceData
            @test forces.n_particles == 2
            @test forces.dims == 2
            @test length(forces.t) == length(sol_coulomb.t) - 1  # Forces need acceleration

            # Should have forces for both (1,2) and (2,1)
            @test haskey(forces.forces, (1, 2))
            @test haskey(forces.forces, (2, 1))
        end

        @testset "Force symmetry (Newton's third law)" begin
            forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c)

            F_12 = forces.forces[(1, 2)]
            F_21 = forces.forces[(2, 1)]

            # F_12 = -F_21
            for t in 1:length(forces.t)
                @test F_12[t] ≈ -F_21[t] atol = 1e-10
            end
        end

        @testset "compute_force_timeseries with stride" begin
            forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c; stride=5)

            # Fewer time points due to stride
            @test length(forces.t) < length(sol_coulomb.t) - 1
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c; stride=0)
            @test_throws ArgumentError compute_force_timeseries(sol_coulomb, 3, 2, masses, charges, c)  # Wrong n_particles
            @test_throws ArgumentError compute_force_timeseries(sol_coulomb, 2, 2, [1.0], charges, c)  # Wrong masses length
            @test_throws ArgumentError compute_force_timeseries(sol_coulomb, 2, 2, masses, [1.0], c)  # Wrong charges length
            @test_throws ArgumentError compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, -1.0)  # Negative c
        end
    end

    @testset "NewtonsThirdLawData" begin
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        c = 1e6

        forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c)

        @testset "check_newtons_third_law" begin
            n3 = check_newtons_third_law(forces)

            @test n3 isa NewtonsThirdLawData
            @test n3.n_pairs == 1  # One pair for 2 particles
            @test length(n3.t) == length(forces.t)
            @test haskey(n3.pair_violations, (1, 2))
            @test haskey(n3.max_violations, (1, 2))
            @test haskey(n3.mean_violations, (1, 2))
            @test haskey(n3.rms_violations, (1, 2))

            # Violations should be small for Coulomb
            @test n3.global_max_violation < 1e-8
        end
    end

    @testset "PhaseSpaceData" begin
        masses = [1.0, 0.5]

        @testset "compute_phase_space_data basic" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses)

            @test ps isa PhaseSpaceData
            @test length(ps.t) == length(sol_coulomb.t)
            @test length(ps.r) == length(sol_coulomb.t)
            @test length(ps.rdot) == length(sol_coulomb.t)
            @test !isnothing(ps.theta)
            @test !isnothing(ps.L)
        end

        @testset "compute_phase_space_data with stride" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses; stride=10)

            expected_points = length(1:10:length(sol_coulomb.t))
            @test length(ps.t) == expected_points
            @test length(ps.r) == expected_points
        end

        @testset "Different particle pairs" begin
            ps_12 = compute_phase_space_data(sol_coulomb, 2, 2, masses; particle_pair=(1, 2))
            ps_21 = compute_phase_space_data(sol_coulomb, 2, 2, masses; particle_pair=(2, 1))

            # r should be the same (distance is symmetric)
            @test ps_12.r ≈ ps_21.r

            # rdot = dot(r_vec, v_vec) / |r|
            # Swapping particles flips both r_vec and v_vec, so rdot is the same
            @test ps_12.rdot ≈ ps_21.rdot
        end

        @testset "Disable optional computations" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses;
                compute_angle=false,
                compute_angular_momentum=false)

            @test isnothing(ps.theta)
            @test isnothing(ps.L)
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_phase_space_data(sol_coulomb, 2, 2, masses; stride=0)
            @test_throws ArgumentError compute_phase_space_data(sol_coulomb, 2, 2, masses; particle_pair=(0, 1))
            @test_throws ArgumentError compute_phase_space_data(sol_coulomb, 2, 2, masses; particle_pair=(1, 3))
        end

        @testset "Positive separation distance" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses)
            @test all(ps.r .> 0)  # Particles should never collide
        end
    end
end
