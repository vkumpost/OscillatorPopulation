@testset "create_data_objective" begin

    model = create_dummy("ode")
    
    t = [5, 6]
    x = [1 2; 5 4]
    cost_function = create_data_objective(model, t, x; input_parameter_name="c")
    err = cost_function([1, 2, 0])
    @test err ≈ 39.25  # sum(([6 7; 10 12] .- [1 5; 2 4]).^2) / 4

    t = [1, 2]
    x = [1 2; 5 4]
    cost_function = create_data_objective(model, t, x; parameter_names=["b"], input_parameter_name="c")
    err = cost_function([1])
    @test err ≈ 2.5  # sum(([2 3; 1 2] .- [1 5; 2 4]).^2) / 4

end

@testset "create_entrainable_oscillator_objective" begin
    
    @test hasmethod(create_entrainable_oscillator_objective,
        Tuple{Any, Any},
        (:trajectories, :parameter_names, :entrainment_extrema,
        :input_parameter_name, :show_plots)
    )

end

@testset "create_desynchronization_objective" begin

    @test hasmethod(create_desynchronization_objective,
        Tuple{Any, Any},
        (:trajectories, :noise_parameter_name, :input_parameter_name, :show_plots)
    )

    @test hasmethod(create_desynchronization_objective,
        Tuple{Any, Any, Any},
        (:trajectories, :noise_parameter_name, :input_parameter_name, :show_plots)
    )

end
