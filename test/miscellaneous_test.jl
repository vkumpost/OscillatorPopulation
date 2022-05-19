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
