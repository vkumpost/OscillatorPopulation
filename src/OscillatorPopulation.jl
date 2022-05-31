"""
Package for a simulation of populations of uncpupled oscillators under periodic
    pacing.

**Data Processing**
- `load_biolum`: Load bioluminescence data.
- `detect_events`: Detect events in the bioluminescence data.
- `biolum_zscore_traces`: Apply Z-score to bioluminescence traces.
- `biolum_mean`: Calculate mean of the bioluminescence traces.

**Events**
- `create_events`: Create events represented by a matrix.
- `create_events_cycle`: Ceate a cycle with specific lengh and period.
- `events_to_function`: Convert the events matrix into a function.
- `plot_events`: Plot events.
- `create_callback`: Create a callback to apply events to a model.
- `detect_events`: Detect events in a binary vector.

**Exceptions**
- `OscillatorPopulationError`: A general exeption if something goes wrong.

**Metrics**
- `estimate_phase_array`: Estimate entrainment phase at each cycle.
- `estimate_phase`: Estimate average entrainment phase.
- `estimate_period`: Estimate period.
- `estimate_winding_number`: Estimate winding number.
- `estimate_winding_number_period`: Estimate period based on the winding number.
- `create_simulation_function`: Generate a simulation function.

**Miscellaneous**
- `smooth`: Smooth a vector using a moving average filter.

**Model Library**
- `load_model`: Load a model from the library.

**Model Functions**
- `set_initial_conditions!`: Set initial conditions.
- `set_timespan!`: Set the time span for the simulation.
- `set_output!`: Set output function.
- `set_solver!`: Set numerical solver and its parameters.
- `get_parameter_index`: Get index of a parameter.
- `get_parameter`: Get value of a parameter.
- `set_parameter!`: Set value of a parameter.
- `set_input!`: Set input events.
- `print_info`: Print model information (initial conditions, parameter values, ...)
- `simulate_model`: Perform a numerical simulation.

**Population Simulation**
- `simulate_population`: Simulate an uncoupled population.
- `plot_solution`: Plot results of a simulation.
- `select_time`: Select specific time of a simulation result.

**Parameter Scans**
- `scan`: Scan parameter values.
- `scan_arnold`: Scan input-signal parameters (amplitude, period, duty cycle).
- `plot_arnold`: Plot the output of `scan_arnold`.
- `estimate_prc`: Estimate a phase response curve.
"""
module OscillatorPopulation

using DataFrames
using DifferentialEquations
using PyPlot
using Random
using Statistics
using StatsBase

import ProgressMeter
import XLSX

include("FindPeaks/FindPeaks.jl")
using .FindPeaks

export load_biolum, detect_events, biolum_zscore_traces, biolum_mean
include("data.jl")

export create_events, create_events_cycle, events_to_function, plot_events,
    create_callback, detect_events
include("events.jl")

export OscillatorPopulationError
include("exceptions.jl")

export estimate_phase_array, estimate_phase, estimate_period,
    estimate_winding_number, estimate_winding_number_period,
    create_simulation_function
include("metrics.jl")

export smooth
include("miscellaneous.jl")

export load_model
include("model_library.jl")

export set_initial_conditions!, set_timespan!, set_output!, set_solver!,
    get_parameter_index, get_parameter, set_parameter!, set_input!, print_info,
    simulate_model
include("model.jl")

export simulate_population, plot_solution, select_time
include("population.jl")

export scan, scan_arnold, plot_arnold, estimate_prc
include("scan.jl")

end  # module
