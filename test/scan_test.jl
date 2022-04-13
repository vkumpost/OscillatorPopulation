@testset "_find_all_combinations" begin

    _find_all_combinations = OscillatorPopulation._find_all_combinations

    vectors = [[1, 0, 2]]
    combinations = _find_all_combinations(vectors)
    @test combinations == hcat([1, 0, 2])

    vectors = [[1, 0, 2], [7, 6]]
    combinations = _find_all_combinations(vectors)
    @test combinations == [1 7; 1 6; 0 7; 0 6; 2 7; 2 6]

    vectors = [[1, 0, 2], ['a'], [7, 6]]
    combinations = _find_all_combinations(vectors)
    @test combinations == [1 'a' 7; 1 'a' 6; 0 'a' 7; 0 'a' 6; 2 'a' 7; 2 'a' 6]

end

@testset "_binary_boundary_search" begin

    _binary_boundary_search = OscillatorPopulation._binary_boundary_search

    # Estimate the target value
    X_target = -13
    fun = x -> x > X_target
    err = 0.1
    x_range = [-17, 84]
    X = _binary_boundary_search(fun, x_range, err)
    @test abs(X_target - X) < err
    
    # Target is outside of the range to the left
    fun = x -> x > -100
    X = _binary_boundary_search(fun, x_range, err)
    @test X == x_range[1]

    # Target is outside of the range to the right
    fun = x -> x > 100
    X = _binary_boundary_search(fun, x_range, err)
    @test X == x_range[2]

end

@testset "scan" begin
    
    model = create_dummy("ode")
    parameters = [
        "a" => [8, 1, 7],
        "c" => [5, 7]
    ]
    simulation_function = function (model=nothing)
        if isnothing(model)
            return ["a+c", "a-c", "b"]
        else
            a = model.problem.p[1]
            b = model.problem.p[2]
            c = model.problem.p[3]
            return [a+c, a-c, b]
        end
    end
    df = scan(model, parameters, simulation_function)
    
    names = OscillatorPopulation.DataFrames.names(df)
    @test names == ["a", "c", "a+c", "a-c", "b"]
    @test df[:, 1] == [ 8,  8,  1,  1,  7,  7]
    @test df[:, 2] == [ 5,  7,  5,  7,  5,  7]
    @test df[:, 3] == [13, 15,  6,  8, 12, 14]
    @test df[:, 4] == [ 3,  1, -4, -6,  2,  0]
    @test df[:, 5] == [ 2,  2,  2,  2,  2,  2]

end

@testset "scan_arnold" begin
    
    model = create_dummy("ode")
    model.problem = remake(model.problem, tspan=(0.0, 20.0))
    simulation_function = function (model=nothing)
        if isnothing(model)
            return ["2nd_event_end", "input_parameter", "input_value"]
        else
            second_events_end = model.input[1][2, 2]
            input_parameter = Float64(model.input[2][1])  # Int('c') == 99
            input_value = model.problem.p[3]
            return [second_events_end, input_parameter, input_value]
        end
    end
    input_amplitudes = [0.5, 1.0, 1.5]
    input_periods = [5.0]
    input_photoperiods = [0.2, 0.5]
    input_parameter = "c"
    df = scan_arnold(model, simulation_function;
        input_amplitudes=input_amplitudes,
        input_periods=input_periods,
        input_photoperiods=input_photoperiods,
        input_parameter=input_parameter)

    names = OscillatorPopulation.DataFrames.names(df)
    @test names == ["input_amplitude", "input_period", "input_photoperiod",
        "2nd_event_end", "input_parameter", "input_value"]
    @test df[:, 1] == [0.5, 0.5, 1.0, 1.0, 1.5, 1.5]
    @test df[:, 2] == [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    @test df[:, 3] == [0.2, 0.5, 0.2, 0.5, 0.2, 0.5]
    @test df[:, 4] == [6.0, 7.5, 6.0, 7.5, 6.0, 7.5]
    @test df[:, 5] == [ 99,  99,  99,  99,  99,  99]  # Int('c') == 99
    @test df[:, 6] == [0.5, 0.5, 1.0, 1.0, 1.5, 1.5]

end

@testset "plot_arnold" begin

    @test plot_arnold isa Function

end

@testset "estimate_prc" begin

    @test estimate_prc isa Function

end
