@testset "create_events" begin

    # Events description is an empty array
    events_description = []
    events = create_events(events_description)
    @test events isa Matrix
    @test isempty(events)

    # Check that light period of :LL followed by LD will merge
    events_description = [(:LL, 2), (:LD, 3, 1, 2), (:DD, 4), (:LL, 1)]
    events = create_events(events_description)
    @test events == [0 3; 5 6; 8 9; 15 16]

    # Tailing DD is ignored
    events_description = [(:DD, 2), (:LL, 3), (:DD, 1)]
    events = create_events(events_description)
    @test events == [2 5]

    # Nonexisting desriptor throws an error
    events_description = [(:DD, 2), (:err, 3)]
    @test_throws OscillatorPopulationError create_events(events_description)

    # Calling the function without input parameters
    events = create_events()
    @test events isa Matrix
    @test isempty(events)

    # Specifying a single-type event using variable args...
    events = create_events(:LD, 5, 2, 1)
    @test events == [0 2; 3 5; 6 8; 9 11; 12 14]

end

@testset "create_events_cycle" begin
    
    # Total duration excatly divisible by period
    events = create_events_cycle(9, 3)
    @test events == [0 1.5; 3 4.5; 6 7.5]

    # Total duration not excatly divisible by period (shorter period)
    events = create_events_cycle(9, 2)
    @test events == [0 1; 2 3; 4 5; 6 7; 8 9]

    # Total duration not excatly divisible by period (longer period)
    events = create_events_cycle(9, 4)
    @test events == [0 2; 4 6; 8 10]

    # Change duty cycle
    events = create_events_cycle(9, 4, 0.75)
    @test events == [0 3; 4 7; 8 11]

    # Duty cycle is 0
    events = create_events_cycle(9, 4, 0)
    @test events isa Matrix
    @test isempty(events)

    # Duty cycle is 1
    events = create_events_cycle(9, 4, 1)
    @test events == [0 9]

    # Duty cycle is outside of the required interval [0, 1]
    @test_throws OscillatorPopulationError create_events_cycle(9, 4, 10)
    @test_throws OscillatorPopulationError create_events_cycle(9, 4, -0.1)

end

@testset "events_to_function" begin
    
    # Convert events to a function
    events = [0 1.5; 3 4.5]
    fun = events_to_function(events)

    # Events starts will evaluate to 1
    @test fun(0) == 1
    @test fun(3) == 1

    # Events ends will evaluate to 0
    @test fun(1.5) == 0
    @test fun(4.5) == 0

    # Times inside of events will evaluate to 1
    @test fun(1) == 1
    @test fun(3.01) == 1

    # Times outside of events will evaluate to 0
    @test fun(-1) == 0
    @test fun(2) == 0
    @test fun(4.6) == 0

    # Empty events will always evaluate to 0
    events = Matrix{Float64}(undef, 0, 2)
    fun = events_to_function(events)
    @test fun.([-1, 0, 1.5, 3.01, 4.6]) == [0, 0, 0, 0, 0]

end

@testset "plot_events" begin
    
    # Check that the function is exported
    @test plot_events isa Function

end

@testset "create_callback" begin
    
    # =========================================================================
    # Continuous callback

    # Create a saving callback to check that the parameter is being modified
    saved_values = SavedValues(Float64, Tuple{Float64,Float64})
    saving_callback = SavingCallback(
        (u, t, integrator) -> (integrator.p[1], integrator.p[2]),
        saved_values,
        saveat=0:0.1:10
    )  
    f = (u, p, t) -> 1.01*u
    u0 = 1/2
    tspan = (0.0, 10.0)
    p = [3, 2]
    prob = ODEProblem(f, u0, tspan, p)

    events = [2 4; 6.5 8.3]
    events_callback = create_callback(events, 1; callback_type="continuous")
    @test events_callback isa CallbackSet
    callback = CallbackSet(saving_callback, events_callback)
    sol = solve(prob, DP5(), saveat=0.1, callback=callback)

    # First parameter is modified
    p1 = [x[1] for x in saved_values.saveval]

    indices = (0 .< saved_values.t .< 2) .|| (4 .< saved_values.t .< 6.5) .||
        (8.3 .< saved_values.t .< 10.0)
    @test all(p1[indices] .== 0)

    indices = (2 .< saved_values.t .< 4) .|| (6.5 .< saved_values.t .< 8.3)
    @test all(p1[indices] .== 3)

    # Second parameter is not modified
    p2 = [x[2] for x in saved_values.saveval]
    @test all(p2 .== 2)

    # =========================================================================
    # Continuous callback - events at the tspan limits

    # Create a saving callback to check that the parameter is being modified
    saved_values = SavedValues(Float64, Tuple{Float64,Float64})
    saving_callback = SavingCallback(
        (u, t, integrator) -> (integrator.p[1], integrator.p[2]),
        saved_values,
        saveat=0:0.1:10
    )       
    f = (u, p, t) -> 1.01*u
    u0 = 1/2
    tspan = (0.0, 10.0)
    p = [3, 2]
    prob = ODEProblem(f, u0, tspan, p)

    events = [0 4; 6.5 10.0]
    events_callback = create_callback(events, 1; callback_type="continuous")
    @test events_callback isa CallbackSet
    callback = CallbackSet(saving_callback, events_callback)
    sol = solve(prob, DP5(), saveat=0.1, callback=callback)

    # First parameter is modified
    p1 = [x[1] for x in saved_values.saveval]

    indices = (4 .< saved_values.t .< 6.5)
    @test all(p1[indices] .== 0)

    indices = (0 .< saved_values.t .< 4) .|| (6.5 .< saved_values.t .< 10.0)
    @test all(p1[indices] .== 3)

    # Second parameter is not modified
    p2 = [x[2] for x in saved_values.saveval]
    @test all(p2 .== 2)

    # =========================================================================
    # Continuous callback - empty events

    # Create a saving callback to check that the parameter is being modified
    saved_values = SavedValues(Float64, Tuple{Float64,Float64})
    saving_callback = SavingCallback(
        (u, t, integrator) -> (integrator.p[1], integrator.p[2]),
        saved_values,
        saveat=0:0.1:10
    )       
    f = (u, p, t) -> 1.01*u
    u0 = 1/2
    tspan = (0.0, 10.0)
    p = [3, 2]
    prob = ODEProblem(f, u0, tspan, p)

    events = Matrix{Float64}(undef, 0, 2)
    events_callback = create_callback(events, 1; callback_type="continuous")
    @test events_callback isa CallbackSet
    callback = CallbackSet(saving_callback, events_callback)
    sol = solve(prob, DP5(), saveat=0.1, callback=callback)

    # First parameter is modified
    p1 = [x[1] for x in saved_values.saveval]

    indices = (0 .< saved_values.t .< 10)
    @test all(p1[indices] .== 0)

    # Second parameter is not modified
    p2 = [x[2] for x in saved_values.saveval]
    @test all(p2 .== 2)

    # =========================================================================
    # Discrete callback

    # Create a saving callback to check that the parameter is being modified
    saved_values = SavedValues(Float64, Tuple{Float64,Float64})
    saving_callback = SavingCallback(
        (u, t, integrator) -> (integrator.p[1], integrator.p[2]),
        saved_values,
        saveat=0:0.1:10
    )       
    f = (u, p, t) -> 1.01*u
    u0 = 1/2
    tspan = (0.0, 10.0)
    p = [3, 2]
    prob = ODEProblem(f, u0, tspan, p)

    events = [2 4; 6.5 8.3]
    events_callback = create_callback(events, 2; callback_type="discrete")
    @test events_callback isa DiscreteCallback
    callback = CallbackSet(saving_callback, events_callback)
    sol = solve(prob, Euler(), dt=0.1, callback=callback)

    # First parameter is not modified
    p1 = [x[1] for x in saved_values.saveval]
    @test all(p1 .== 3)

    # First parameter is modified
    p2 = [x[2] for x in saved_values.saveval]

    indices = (0 .< saved_values.t .< 2) .|| (4 .< saved_values.t .< 6.5) .||
        (8.3 .< saved_values.t .< 10.0)
    @test all(p2[indices] .== 0)
    
    indices = (2 .< saved_values.t .< 4) .|| (6.5 .< saved_values.t .< 8.3)
    @test all(p2[indices] .== 2)

    # =========================================================================
    # Unknown callback type
    
    events = [1 2]
    i = 1
    @test_throws OscillatorPopulationError create_callback(events, i;
        callback_type="discrete_02")

end

@testset "detect_events" begin

    # 0s at the beginning and the end
    t = [0, 1, 2, 3, 4, 5, 6]
    x = [0, 0, 1, 1, 0, 1, 0]
    events = detect_events(t, x)
    @test events ≈ [1.5 3.5; 4.5 5.5]

    # 1s at the beginning and 0s at the end
    t = [0, 1, 2, 3, 4, 5, 6]
    x = [1, 1, 0, 1, 0, 1, 0]
    events = detect_events(t, x)
    @test events ≈ [0 1.5; 2.5 3.5; 4.5 5.5]

    # 0s at the beginning and 1s at the end
    t = [0, 1, 2, 3, 4, 5, 6]
    x = [0, 0, 0, 1, 0, 1, 1]
    events = detect_events(t, x)
    @test events ≈ [2.5 3.5; 4.5 6]

    # 1s at the beginning and the end
    t = [0, 1, 2, 3, 4, 5, 6]
    x = [1, 0, 1, 1, 0, 1, 1]
    events = detect_events(t, x)
    @test events ≈ [0 0.5; 1.5 3.5; 4.5 6]

    # One big event
    t = [0, 1, 2, 3, 4, 5, 6]
    x = [1, 1, 1, 1, 1, 1, 1]
    events = detect_events(t, x)
    @test events ≈ [0 6]

    # No event
    t = [0, 1, 2, 3, 4, 5, 6]
    x = [0, 0, 0, 0, 0, 0, 0]
    events = detect_events(t, x)
    @test events isa Matrix
    @test isempty(events)

end
