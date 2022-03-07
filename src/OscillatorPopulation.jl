module OscillatorPopulation

using DifferentialEquations


export OscillatorPopulationError
include("exceptions.jl")

export load_model
include("model_library.jl")

export simulate_model
include("model.jl")


end  # module
