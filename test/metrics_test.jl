@testset "estimate_phase" begin
    
    t = 1:20
    x_original = fill(0.0, length(t))
    events = [0 2; 4 6; 8 10; 12 14; 16 18]

    # Entrained signal
    x = copy(x_original)
    x[7] = 1
    x[11] = 1
    x[15] = 1

    phase, phase_error = estimate_phase(t, x, events)
    @test phase ≈ 0.75
    @test phase_error ≈ 0

    phase, phase_error = estimate_phase(t, x, events; normalize=false)
    @test phase ≈ 3
    @test phase_error ≈ 0

    # Skipped period
    x = copy(x_original)
    x[7] = 1
    x[11] = 0
    x[15] = 1

    phase, phase_error = estimate_phase(t, x, events)
    @test isnan(phase)
    @test isnan(phase_error)

    # Unentrained signal
    x = copy(x_original)
    x[5] = 1
    x[10] = 1
    x[15] = 1

    phase, phase_error = estimate_phase(t, x, events)
    @test phase ≈ 0.5
    @test phase_error ≈ 0.25

    # Multiple peaks
    x = copy(x_original)
    x[5] = 1
    x[7] = 0.5
    x[10] = 1
    x[13] = 0.5
    x[15] = 1

    phase, phase_error = estimate_phase(t, x, events)
    @test phase ≈ 0.5
    @test phase_error ≈ 0.25

end

@testset "create_simulation_function" begin
    
    model = create_dummy("ode")

    fun = create_simulation_function(; transient=0.5, variable=2)
    names = fun()
    @test names == ["minimum", "maximum", "phase", "phase_error"]
    
    properties = fun(model)
    @test properties[1] == 4
    @test properties[2] == 8
    @test isnan(properties[3])
    @test isnan(properties[4])

end
