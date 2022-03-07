module OscillatorPopulation

using DifferentialEquations


export OscillatorPopulationError
include("exceptions.jl")

export load_model
include("model_library.jl")

include("model.jl")


end # module
