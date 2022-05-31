@testset "estimate_phase_array" begin
    
    t = 1:20
    x_original = fill(0.0, length(t))
    events = [0 2; 4 6; 8 10; 12 14; 16 18]

    # Entrained signal
    x = copy(x_original)
    x[7] = 1
    x[11] = 1
    x[15] = 1

    phase_arr = estimate_phase_array(t, x, events)
    @test phase_arr ≈ [0.75, 0.75, 0.75]

    phase_arr = estimate_phase_array(t, x, events; normalize=false)
    @test phase_arr ≈ [3, 3, 3]

    # Skipped period
    x = copy(x_original)
    x[7] = 1
    x[11] = 0
    x[15] = 1

    phase_arr = estimate_phase_array(t, x, events)
    @test phase_arr[1] ≈ 0.75
    @test isnan(phase_arr[2])
    @test phase_arr[3] ≈ 0.75

    # Unentrained signal
    x = copy(x_original)
    x[5] = 1
    x[10] = 1
    x[15] = 1

    phase_arr = estimate_phase_array(t, x, events)
    @test phase_arr ≈ [0.25, 0.5, 0.75]

    # Multiple peaks
    x = copy(x_original)
    x[5] = 1
    x[7] = 0.5
    x[10] = 1
    x[13] = 0.5
    x[15] = 1

    phase_arr = estimate_phase_array(t, x, events)
    @test phase_arr ≈ [0.25, 0.5, 0.75]

end

@testset "estimate_order_parameter" begin
    
    # Throw error if the phase is outside the allowed interval [0, 1]
    phase_array = [0.1, -0.1, 0.5]
    @test_throws OscillatorPopulationError estimate_order_parameter(phase_array)

    phase_array = [0.1, 0.5, 1.5]
    @test_throws OscillatorPopulationError estimate_order_parameter(phase_array)

    # Coherent population
    phase_array = [0.2, 0.2, 0.2]
    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
    @test phase_coherence ≈ 1
    @test collective_phase ≈ 0.2

    # Incoherent population
    phase_array = [0.0, 0.25, 0.5, 0.75]
    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
    @test phase_coherence < 1e-10

    # Slightly desynchronized population
    phase_array = [0.2, 0.3, 0.4]
    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
    @test phase_coherence ≈ 0.872677996249965
    @test collective_phase ≈ 0.3

end


@testset "estimate_period" begin

    t = 0:0.01:20
    x = sin.(2π .* (1/0.1) .* t)
    T, _ = estimate_period(t, x)
    @test T ≈ 0.1

    x = sin.(2π .* (1/5) .* t)
    T, _ = estimate_period(t, x)
    @test T ≈ 5

    T, _ = estimate_period(t, x; n_lags=1000)
    @test T ≈ 5

    T, Te = estimate_period(t, x; n_lags=400)
    @test isnan(T)
    @test isnan(Te)

end

@testset "estimate_winding_number" begin

    # Half loop
    x = [0, 1, 0]
    y = [-1, 0, 1]
    winding_number = estimate_winding_number(x, y)
    @test winding_number ≈ 0.5

    # Half loop in the opposite direction
    x = [0, -1, 0]
    y = [-1, 0, 1]
    winding_number = estimate_winding_number(x, y)
    @test winding_number ≈ 0.5

    # One loop
    x = [ 0, 1, 0, -1,  0]
    y = [-1, 0, 1,  0, -1]
    winding_number = estimate_winding_number(x, y)
    @test winding_number ≈ 1

    # One loop around, one outside
    x = [0, 1, 0, -1,  0,  1,  1,  0,  0]
    y = [-1, 0, 1, 0, -1, -1, -2, -2, -1]
    winding_number = estimate_winding_number(x, y)
    @test winding_number ≈ 1

end

@testset "estimate_winding_number_period" begin
    
    # Half loop
    x = [0, 1, 0]
    y = [-1, 0, 1]
    period = estimate_winding_number_period(x, y, 0.5)
    @test period ≈ 1
    period = estimate_winding_number_period(x, y, 2)
    @test period ≈ 4
    period = estimate_winding_number_period(x, y, 0.25)
    @test period ≈ 0.5

end

@testset "create_simulation_function" begin
    
    model = create_dummy("ode")
    events = create_events_cycle(10, 0.5)
    set_input!(model, events, "c")
    fun = create_simulation_function(; transient=0.5, variable=2, variable_2=1)
    names = fun()
    @test names == ["minimum", "maximum", "winding_number", "phase_coherence",
        "collective_phase"]
    
    properties = fun(model)
    @test properties[1] == 4
    @test properties[2] == 8
    @test 0 < properties[3]
    @test isnan(properties[4])
    @test isnan(properties[5])

end
