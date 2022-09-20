@testset "smooth" begin
    
    # Default span is 5
    y = [-2, -5, 1, -1, 3, 0, 6]
    yy = smooth(y)
    @test all(yy .≈ [-2, -2, -0.8, -0.4, 1.8, 3, 6])

    # Span 1 retunrs the same vector
    yy = smooth(y, span=1)
    @test all(yy .≈ [-2, -5, 1, -1, 3, 0, 6])

    # Setting span longer that the vector length
    yy = smooth(y, span=100)
    @test all(yy .≈ [-2, -2, -0.8, 0.2857142857142857, 1.8, 3, 6])

end

@testset "generate_random_values" begin

    seed = 283
    μ = 9
    σ = 4
    n = 1_000_000
    values = generate_random_values(μ, σ, n; seed=seed)
    @test abs(mean(values) - μ) < 0.01
    @test abs(std(values) - σ) < 0.01
    @test length(values) == n

    # Same seed same results
    values_new = generate_random_values(μ, σ, n; seed=seed)
    @test values == values_new

    # Different seed different results
    values_new = generate_random_values(μ, σ, n; seed=seed+1)
    @test values != values_new

    # Tructunated population
    values = generate_random_values(μ, σ, n; seed=seed, lower=8, upper=12)
    @test all(8 .< values .< 12)

end

@testset "benchmark" begin

    @test hasmethod(benchmark, Tuple{Function}, (:n_repeats,))

end
