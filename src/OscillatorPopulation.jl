module OscillatorPopulation


using DifferentialEquations
using PyPlot


export create_events, create_events_cycle, events_to_function, plot_events,
    create_callback
include("events.jl")

export OscillatorPopulationError
include("exceptions.jl")

export load_model
include("model_library.jl")

export set_initial_conditions!, set_timespan!, set_output!, set_solver!,
    get_parameter_index, get_parameter, set_parameter!, set_input!,
    simulate_model
include("model.jl")


end  # module
