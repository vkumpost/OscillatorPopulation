module OscillatorPopulation


using DifferentialEquations
using PyPlot


export create_events, create_events_cycle, plot_events
include("events.jl")

export OscillatorPopulationError
include("exceptions.jl")

export load_model
include("model_library.jl")

export simulate_model
include("model.jl")


end  # module
