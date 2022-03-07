using DifferentialEquations
using OscillatorPopulation
using Test


@testset "exceptions" begin
    include("exceptions_test.jl")
end

@testset "model_library" begin
    include("model_library_test.jl")
end
