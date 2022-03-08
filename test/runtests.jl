using DifferentialEquations
using OscillatorPopulation
using Test

@testset "events" begin
    include("events_test.jl")
end

@testset "exceptions" begin
    include("exceptions_test.jl")
end

@testset "model_library" begin
    include("model_library_test.jl")
end

@testset "model" begin
    include("model_test.jl")
end
