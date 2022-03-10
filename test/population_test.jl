@testset "simulate_population" begin

    t_expected = [0.0, 1.0, 2.0, 3.0, 4.0]
    x_expected = [1.0 0.0; 2.0 2.0; 3.0 4.0; 4.0 6.0; 5.0 8.0]

    # Simulate population
    model = create_dummy("ode")
    solution = simulate_population(model, 2)
    @test solution.time == t_expected
    @test solution.mean == x_expected
    @test solution.trajectories[:, :, 1] == x_expected
    @test solution.trajectories[:, :, 2] == x_expected
    @test isempty(solution.events)
    @test solution.success

    # Save only mean
    model = create_dummy("ode")
    solution = simulate_population(model, 2; save_trajectories=false)
    @test solution.time == t_expected
    @test solution.mean == x_expected
    @test isempty(solution.trajectories)
    @test isempty(solution.events)
    @test solution.success

end

@testset "plot_solution" begin

    @test plot_solution isa Function

end
