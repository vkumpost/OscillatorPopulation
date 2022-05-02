using DataFrames
using DifferentialEquations
using OscillatorPopulation
using StatsBase
using Test

@testset "data" begin
    include("data_test.jl")
end

@testset "dummy_model" begin
    include("dummy_model.jl")
end

@testset "events" begin
    include("events_test.jl")
end

@testset "exceptions" begin
    include("exceptions_test.jl")
end

@testset "metrics" begin
    include("metrics_test.jl")
end

@testset "model_library" begin
    include("model_library_test.jl")
end

@testset "model" begin
    include("model_test.jl")
end

@testset "population" begin
    include("population_test.jl")
end

@testset "scan" begin
    include("scan_test.jl")
end
