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

@testset "damped_sine" begin

    t = [0.2, 1.1, 2.8]
    p = [3, 2, 0.5, 1.5]
    x = damped_sine(t, p)
    @test x ≈ [-1.539213246532194, 0.12482573436989124, -0.009413673996984902]

end

@testset "polynomial" begin

    t = [0.8, 1.4, 2.9]
    p = Float64[]
    x = polynomial(t, p)
    @test x ≈ [0.0, 0.0, 0.0]

    t = [0.8, 1.4, 2.9]
    p = [3]
    x = polynomial(t, p)
    @test x ≈ [3.0, 3.0, 3.0]

    t = [0.8, 1.4, 2.9]
    p = [2, 6, 1]
    x = polynomial(t, p)
    @test x ≈ [7.44, 12.36, 27.81]

end

@testset "fit_curve" begin
    
    t = 0:0.01:10
    p = [3, 0.1, 0.5, 1.5]
    x = damped_sine(t, p)
    p0 = [2.9, 0.15, 0.6, 1.4]
    p_estimate = fit_curve(damped_sine, t, x, p0)
    @test p ≈ p_estimate

end

@testset "binary_search" begin
    
    fun = x -> x ^ 2
    target_value = 144
    search_range = [0, 100]
    x = binary_search(fun, target_value, search_range; tolerance=0)
    @test x ≈ 12

    # Target value is outside of the search range
    target_value = 61504
    @test_throws OscillatorPopulationError binary_search(fun, target_value, search_range)

end
