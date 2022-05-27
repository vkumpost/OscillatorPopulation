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
