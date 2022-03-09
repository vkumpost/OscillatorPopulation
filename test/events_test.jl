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

    # Change photoperiod
    events = create_events_cycle(9, 4, 0.75)
    @test events == [0 3; 4 7; 8 11]

    # Photoperiod is 0
    events = create_events_cycle(9, 4, 0)
    @test events isa Matrix
    @test isempty(events)

    # Photoperiod is 1
    events = create_events_cycle(9, 4, 1)
    @test events == [0 9]

    # Photoperiod is outside of the required interval [0, 1]
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
