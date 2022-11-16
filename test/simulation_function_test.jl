@testset "create_simulation_function" begin
    
    model = create_dummy("ode")
    events = create_events_cycle(10, 0.5)
    set_input!(model, events, "c")
    set_solver!(model, saveat=1.0)
    fun = create_simulation_function(; transient=0.5, variable_x=2, variable_y=1)
    names = fun()
    @test names == ["minimum", "maximum", "amplitude", "rms", "winding_number",
        "autocorrelation", "phase_coherence", "mean_phase",
        "phase_coherence_cxcorr", "mean_phase_cxcorr",
        "phase_coherence_population", "collective_phase",
        "phase_coherence_population_cxcorr", "collective_phase_cxcorr"
    ]
    
    properties = fun(model)
    @test properties[1] == 4
    @test properties[2] == 8
    @test properties[3] == 4

    fun = create_simulation_function(["rms", "mean_phase_cxcorr"])
    names = fun()
    @test names == ["rms", "mean_phase_cxcorr"]

end
