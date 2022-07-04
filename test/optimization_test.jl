@testset "estimate_initial_conditions" begin

    model = create_dummy("ode")
    parameters = (
        ["a", "b"],
        [0 0; 1 3; -1 2]
    )

    # Free-running oscillator
    initial_conditions = estimate_initial_conditions(model, 3, (:DD, 2);
        input_parameter_name="c", parameters=parameters)
    @test initial_conditions ≈ [
        1.0 0.0;
        3.0 6.0;
       -1.0 4.0
    ]

    # Entrained oscillator
    initial_conditions = estimate_initial_conditions(model, 3, (:LD, 2, 2, 2);
        input_parameter_name="a", parameters=parameters)
    @test initial_conditions ≈ [
        1.0 0.0;
        5.0 24.0;
       -3.0 16.0
    ]

end

@testset "create_data_objective" begin

    model = create_dummy("ode")
    
    t = [5, 6]
    x = [1 2; 5 4]
    cost_function = create_data_objective(model, t, x; input_parameter_name="c")
    err = cost_function([1, 2, 0])
    @test err ≈ 157  # sum(([6 7; 10 12] .- [1 5; 2 4]).^2)

    t = [1, 2]
    x = [1 2; 5 4]
    cost_function = create_data_objective(model, t, x; parameter_names=["b"], input_parameter_name="c")
    err = cost_function([1])
    @test err ≈ 10  # sum(([2 3; 1 2] .- [1 5; 2 4]).^2)

end

@testset "optimize" begin

    cost_function = x -> (x[1] - 5)^2 + (x[2] - 4)^2
    
    search_range = [(4, 6), (3, 5)]
    initial_population = [5.1 4.1; 4.9 3.9; 4.8 4.2; 5.2 3.8]
    best_candidate, final_population = optimize(cost_function;
        search_range=search_range,
        max_steps=10_000,
        initial_population=initial_population,
        trace_mode=:silent
    )
    @test sum((best_candidate .- [5, 4]).^2) < 0.1
    @test size(final_population) == (4, 2)

    search_range = [(4, 6), (3, 5)]
    best_candidate, final_population = optimize(cost_function;
        search_range=search_range,
        max_steps=10_000,
        trace_mode=:silent
    )
    @test sum((best_candidate .- [5, 4]).^2) < 0.1
    @test size(final_population) == (50, 2)

end
