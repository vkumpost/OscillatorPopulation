@testset "rsquared" begin
    
    x = [1, 2, 3, 4]
    R = rsquared(x, x)
    @test R ≈ 1

    x = [1, 2, 3, 4]
    y = [1, 2, 2, 3]
    R = rsquared(x, y)
    @test 0 < R < 1

    x = [0, 1, 0, -1, 0]
    y = [0, 0, 0,  0, 0]
    R = rsquared(x, y)
    @test R ≈ 0

    x = [0, 1, 0, -1, 0]
    y = [0, 0, -2, 20, 50]
    R = rsquared(x, y)
    @test R < 0

end

@testset "cmean" begin

    x = [-2, -1, 0, 1, 2, 3]
    m = cmean(x)
    @test abs(m) < 1e-15

    x = [0.25, 0.75, 1.0, 0.0]
    m = cmean(x)
    @test abs(m) < 1e-15

    x = [-0.5, 0.5, 1.5]
    m = cmean(x)
    @test abs(m - 0.5) < 1e-15

end

@testset "cstd" begin

    x = [-2.1, -1.1, -0.1, 0.9, 1.9]
    sd = cstd(x)
    @test sd ≈ 0

    sd1 = cstd([-0.4, 0.5, 2.6])
    sd2 = cstd([0.3, 0.5, 0.8])
    @test sd1 < sd2

end

@testset "window_xcorr" begin

    # Cross-correlation
    x = [0, 1, 2, -1]
    y = [0, 1, 1, 0, 0, 5, 8, 12, -2, 3, 0, 0, 1, 2]
    r = window_xcorr(x, y)
    @test r ≈ [0.41138392038293375, 0, -1.6455356815317348, -0.925613820861601,
        -0.7199218606701341, 4.628069104308004, -1.1313057810530678,
        -0.5142299004786671, 0.5142299004786671, -0.6170758805744005,
        -0.30853794028720033] ./ 4

    x = [0, 1, 2, -1]
    y = [1, 2, 1, 5]
    r = window_xcorr(x, y)
    @test r ≈ [-2.2505813202525626] ./ 4

    x = [0, 1, 2, -1]
    y = [1, 2, 1]
    @test_throws OscillatorPopulationError window_xcorr(x, y)

end

@testset "cxcorr" begin

    x = [0, 0, 1, 1, 0]
    y = [0, 1, 1, 0, 0]
    r = cxcorr(x, y)
    @test r ≈ [0.6666666666666663, -2.6666666666666665, -2.6666666666666665,
        0.6666666666666663, 3.999999999999999] ./ 5

    r = cxcorr(y, x)
    @test r ≈ [0.6666666666666663, 3.999999999999999, 0.6666666666666663,
        -2.6666666666666665, -2.6666666666666665] ./ 5

end

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

    # Pass a matrix
    X = fill(0.0, length(t), 4)
    X[7, 1] = 1
    X[11, 1] = 1
    X[15, 1] = 1
    X[7, 2] = 1
    X[11, 2] = 0
    X[15, 2] = 1
    X[5, 3] = 1
    X[7, 3] = 0.5
    X[10, 3] = 1
    X[13, 3] = 0.5
    X[15, 3] = 1
    X[5, 4] = 1
    X[7, 4] = 0.5
    X[10, 4] = 1
    X[13, 4] = 0.5
    X[15, 4] = 1
    phase_arr = estimate_phase_array(t, X, events)
    phase_arr[isnan.(phase_arr)] .= -0.1  # replace NaNs for easier test evaluation
    @test phase_arr ≈ [
        0.75 0.75 0.25 0.25;
        0.75 -0.1 0.50 0.50;
        0.75 0.75 0.75 0.75
    ]

end

@testset "estimate_phase_array_peaks" begin

    @test estimate_phase_array_peaks isa Function

end

@testset "estimate_phase_array_cxcorr" begin

    @test estimate_phase_array_cxcorr isa Function

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

    # Matrix input
    phase_array = [
        0.2 0.2 0.2;
        0.2 0.3 0.4
    ]
    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
    @test phase_coherence ≈ (0.872677996249965 + 1)/2
    @test collective_phase ≈ 0.25

    # NaNs
    phase_array = [0.2, NaN, 0.3, 0.4]
    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
    @test phase_coherence ≈ 0.872677996249965
    @test collective_phase ≈ 0.3

    phase_coherence, collective_phase = estimate_order_parameter(phase_array; skip_nan=false)
    @test isnan(phase_coherence)
    @test isnan(collective_phase)

    phase_array = [NaN, NaN]
    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
    @test isnan(phase_coherence)
    @test isnan(collective_phase)

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
    winding_number = estimate_winding_number(x, y; remove_mean=false)
    @test winding_number ≈ 0.5

    # Half loop in the opposite direction
    x = [0, -1, 0]
    y = [-1, 0, 1]
    winding_number = estimate_winding_number(x, y; remove_mean=false)
    @test winding_number ≈ 0.5

    # One loop
    x = [ 0, 1, 0, -1,  0]
    y = [-1, 0, 1,  0, -1]
    winding_number = estimate_winding_number(x, y; remove_mean=false)
    @test winding_number ≈ 1

    # One loop around, one outside
    x = [0, 1, 0, -1,  0,  1,  1,  0,  0]
    y = [-1, 0, 1, 0, -1, -1, -2, -2, -1]
    winding_number = estimate_winding_number(x, y; remove_mean=false)
    @test winding_number ≈ 1

    # remove_mean option
    x = [ 0, 1, 0, -1,  0] .+ 10
    y = [-0.5, 0, 1,  0, -0.5]
    winding_number = estimate_winding_number(x, y; remove_mean=true)
    @test winding_number ≈ 1
    winding_number = estimate_winding_number(x, y; remove_mean=false)
    @test winding_number < 1e-10

end

@testset "estimate_period_winding_number" begin
    
    # Half loop
    x = [0, 1, 0]
    y = [-1, 0, 1]
    period = estimate_period_winding_number(x, y, 0.5; remove_mean=false)
    @test period ≈ 1
    period = estimate_period_winding_number(x, y, 2; remove_mean=false)
    @test period ≈ 4
    period = estimate_period_winding_number(x, y, 0.25; remove_mean=false)
    @test period ≈ 0.5

    # remove_mean option
    x = [ 0, 1, 0, -1,  0] .+ 10
    y = [-0.5, 0, 1,  0, -0.5]
    period = estimate_period_winding_number(x, y, 2; remove_mean=true)
    @test period ≈ 2
    period = estimate_period_winding_number(x, y, 2; remove_mean=false)
    @test period > 2

end

@testset "create_simulation_function" begin
    
    model = create_dummy("ode")
    events = create_events_cycle(10, 0.5)
    set_input!(model, events, "c")
    set_solver!(model, saveat=1.0)
    fun = create_simulation_function(; transient=0.5, variable_x=2, variable_y=1)
    names = fun()
    @test names == ["minimum", "maximum", "amplitude", "rms",
        "winding_number", "autocorrelation",
        "phase_coherence", "mean_phase",
        "phase_coherence_cxcorr", "mean_phase_cxcorr",
        "phase_coherence_population", "collective_phase"]
    
    properties = fun(model)
    @test properties[1] == 4
    @test properties[2] == 8
    @test properties[3] == 4

    fun = create_simulation_function(["rms", "mean_phase_cxcorr"])
    names = fun()
    @test names == ["rms", "mean_phase_cxcorr"]

end
