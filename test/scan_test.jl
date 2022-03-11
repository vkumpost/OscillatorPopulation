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
