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
