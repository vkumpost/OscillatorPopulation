module OscillatorPopulation

using DifferentialEquations
using PyPlot

import ProgressMeter

export create_events, create_events_cycle, events_to_function, plot_events,
    create_callback
include("events.jl")

export OscillatorPopulationError
include("exceptions.jl")

export load_model
include("model_library.jl")

export set_initial_conditions!, set_timespan!, set_output!, set_solver!,
    get_parameter_index, get_parameter, set_parameter!, set_input!, print_info,
    simulate_model
include("model.jl")

export simulate_population, plot_solution, select_time
include("population.jl")

include("scan.jl")

end  # module
