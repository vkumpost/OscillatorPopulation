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
