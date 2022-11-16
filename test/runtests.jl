using DataFrames
using DifferentialEquations
using OscillatorPopulation
using StatsBase
using Test

@testset "biolum_data" begin
    include("biolum_data_test.jl")
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

@testset "input_output" begin
    include("input_output_test.jl")
end

@testset "metrics" begin
    include("metrics_test.jl")
end

@testset "miscellaneous" begin
    include("miscellaneous_test.jl")
end

@testset "model_library" begin
    include("model_library_test.jl")
end

@testset "model" begin
    include("model_test.jl")
end

@testset "optimization" begin
    include("optimization_test.jl")
end

@testset "population" begin
    include("population_solution_test.jl")
end

@testset "scan" begin
    include("scan_test.jl")
end

@testset "simulation_function" begin
    include("simulation_function_test.jl")
end
