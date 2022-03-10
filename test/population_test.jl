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

@testset "select_time" begin

    t = [0.0, 1.0, 2.0]
    m = [1.0 0.0; 2.0 2.0; 3.0 4.0]
    U = Array{Float64, 3}(undef, 3, 2, 2)
    U[:, :, 1] = [0.5 0.0; 1.5 1.5; 2.5 3.5]
    U[:, :, 2] = [1.5 0.0; 2.5 2.5; 3.5 4.5]
    events = [0.0 1.0; 1.2 1.8]
    
    # Remove offset
    t1 = [5.0, 6.0, 7.0]
    events1 = [5.0 6.0; 6.2 6.8]
    solution = OscillatorPopulation.PopulationSolution(t1, m, U, events1, true)
    @test solution.time == t1
    @test solution.events == events1
    solution2 = select_time(solution, remove_offset=true)
    @test all(solution2.time .≈ t)
    @test all(solution2.mean .≈ m)
    @test all(solution2.trajectories .≈ U)
    @test all(solution2.events .≈ events)
    @test all(solution2.success)
    
    # Set offset
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, offset=3.6)
    @test all(solution2.time .≈ t .+ 3.6)
    @test all(solution2.mean .≈ m)
    @test all(solution2.trajectories .≈ U)
    @test all(solution2.events .≈ events .+ 3.6)
    @test all(solution2.success)
    
    # Min time (in light)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, min_time=1.5, remove_offset=false)
    @test all(solution2.time .≈ [2.0])
    @test all(solution2.mean .≈ [3.0 4.0])
    @test all(solution2.trajectories .≈ U[3:3, :, :])
    @test all(solution2.events .≈ [1.5 1.8])
    @test all(solution2.success)
    
    # Min time (in darkness)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, min_time=1.1, remove_offset=false)
    @test all(solution2.time .≈ [2.0])
    @test all(solution2.mean .≈ [3.0 4.0])
    @test all(solution2.trajectories .≈ U[3:3, :, :])
    @test all(solution2.events .≈ [1.2 1.8])
    @test all(solution2.success)
    
    # Min time (light onset)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, min_time=1.2, remove_offset=false)
    @test all(solution2.time .≈ [2.0])
    @test all(solution2.mean .≈ [3.0 4.0])
    @test all(solution2.trajectories .≈ U[3:3, :, :])
    @test all(solution2.events .≈ [1.2 1.8])
    @test all(solution2.success)
    
    # Min time (light offset)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, min_time=1.0, remove_offset=false)
    @test all(solution2.time .≈ [1.0, 2.0])
    @test all(solution2.mean .≈ [2.0 2.0; 3.0 4.0])
    @test all(solution2.trajectories .≈ U[2:3, :, :])
    @test all(solution2.events .≈ [1.2 1.8])
    @test all(solution2.success)
    
    # Max time (in light)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, max_time=1.5, remove_offset=false)
    @test all(solution2.time .≈ [0.0, 1.0])
    @test all(solution2.mean .≈ [1.0 0.0; 2.0 2.0])
    @test all(solution2.trajectories .≈ U[1:2, :, :])
    @test all(solution2.events .≈ [0.0 1.0; 1.2 1.5])
    @test all(solution2.success)
    
    # Max time (in darkness)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, max_time=1.1, remove_offset=false)
    @test all(solution2.time .≈ [0.0, 1.0])
    @test all(solution2.mean .≈ [1.0 0.0; 2.0 2.0])
    @test all(solution2.trajectories .≈ U[1:2, :, :])
    @test all(solution2.events .≈ [0.0 1.0])
    @test all(solution2.success)
    
    # Max time (light onset)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, max_time=1.2, remove_offset=false)
    @test all(solution2.time .≈ [0.0, 1.0])
    @test all(solution2.mean .≈ [1.0 0.0; 2.0 2.0])
    @test all(solution2.trajectories .≈ U[1:2, :, :])
    @test all(solution2.events .≈ [0.0 1.0])
    @test all(solution2.success)
    
    # Max time (light offset)
    solution = OscillatorPopulation.PopulationSolution(t, m, U, events, true)
    solution2 = select_time(solution, max_time=1.8, remove_offset=false)
    @test all(solution2.time .≈ [0.0, 1.0])
    @test all(solution2.mean .≈ [1.0 0.0; 2.0 2.0])
    @test all(solution2.trajectories .≈ U[1:2, :, :])
    @test all(solution2.events .≈ [0.0 1.0; 1.2 1.8])
    @test all(solution2.success)

end
