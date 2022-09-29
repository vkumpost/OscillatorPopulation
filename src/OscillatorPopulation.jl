"""
Package for a simulation of populations of uncpupled oscillators under periodic
    pacing.

**Data Processing**
- `load_biolum`: Load bioluminescence data.
- `detect_events`: Detect events in the bioluminescence data.
- `biolum_zscore_traces`: Apply Z-score to bioluminescence traces.
- `biolum_mean`: Calculate mean of the bioluminescence traces.
- `save_data`: Save a dataframe as a csv file.
- `load_data`: Load a dataframe from a csv file.
- `save_figure`: Save figure as a file.

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
- `rsquared`: The coefficient of determination.
- `cmean`: Circular mean.
- `cstd`: Circular standard deviation.
- `window_xcorr`: Moving-window cross-correlation.
- `cxcorr`: Circular cross-correlation.
- `estimate_phase_array`: Estimate entrainment phase at each cycle.
- `estimate_phase_array_peaks`: Estimate entrainment phase using peak detection.
- `estimate_phase_array_cxcorr`: Estimate entrainment phase using circular
    cross-correlation.
- `estimate_order_parameter`: Estimate Kuramoto's order parameter.
- `estimate_period`: Estimate period.
- `estimate_winding_number`: Estimate winding number.
- `estimate_period_winding_number`: Estimate period based on the winding number.
- `create_simulation_function`: Generate a simulation function.

**Miscellaneous**
- `smooth`: Smooth a vector using a moving average filter.
- `generate_random_values`: Generate a vector of random numbers.
- `benchmark`: Estimate evaluation time for a function.
- `find_closest`: Find the index of the closest value in an array.

**Model Library**
- `kfr`: Kim-Forger function.
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

**Optimization**
- `estimate_initial_conditions`: Preestimate initial conditions for a model.
- `create_data_objective`: Create a cost function to fit a model to data.
- `create_entrainable_oscillator_objective`: A cost function to find parameters for an entrainmeble oscillator.
- `create_desynchronization_objective`: A cost function to find noise intensity based on the population-level damping rate.
- `optimize`: Minimize a cost function.
- `damped_sine`: A damped sine function.
- `polynomial`: Polynomial function.
- `fit_curve`: Fit curve to data.

**Population Simulation**
- `simulate_population`: Simulate an uncoupled population.
- `plot_solution`: Plot results of a simulation.
- `select_time`: Select specific time of a simulation result.
- `select_subset`: Select a subset of trajectories.

**Parameter Scans**
- `scan`: Scan parameter values.
- `scan_arnold`: Scan input-signal parameters (amplitude, period, duty cycle).
- `plot_arnold`: Plot the output of `scan_arnold`.
- `estimate_prc`: Estimate a phase response curve.
- `estimate_T_prc`: Estimate T-cycle phase response curve.
"""
module OscillatorPopulation

using DataFrames
using DelimitedFiles
using DifferentialEquations
using Distributions
using FFTW
using PyPlot
using Random
using Statistics
using StatsBase

import BlackBoxOptim
import LsqFit
import ProgressMeter
import XLSX

include("FindPeaks/FindPeaks.jl")
using .FindPeaks

export load_biolum, detect_events, biolum_zscore_traces, biolum_mean, save_data,
    load_data, save_figure
include("data.jl")

export create_events, create_events_cycle, events_to_function, plot_events,
    create_callback, detect_events
include("events.jl")

export OscillatorPopulationError
include("exceptions.jl")

export rsquared, cmean, cstd, window_xcorr, cxcorr, estimate_phase_array,
    estimate_phase_array_peaks, estimate_phase_array_cxcorr,
    estimate_order_parameter, estimate_period, estimate_winding_number,
    estimate_period_winding_number, create_simulation_function
include("metrics.jl")

export smooth, generate_random_values, benchmark, find_closest
include("miscellaneous.jl")

export kfr, load_model
include("model_library.jl")

export set_initial_conditions!, set_timespan!, set_output!, set_solver!,
    get_parameter_index, get_parameter, set_parameter!, set_input!, print_info,
    simulate_model
include("model.jl")

export estimate_initial_conditions, create_data_objective,
    create_entrainable_oscillator_objective, create_desynchronization_objective,
    optimize, damped_sine, polynomial,
    fit_curve
include("optimization.jl")

export simulate_population, plot_solution, select_time, select_subset
include("population.jl")

export scan, scan_arnold, plot_arnold, estimate_prc, estimate_T_prc
include("scan.jl")

end  # module
