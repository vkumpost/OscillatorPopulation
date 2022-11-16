# OscillatorPopulation

Package for simulation of populations of deterministic and stochastic oscillators under external forcing. The package implements various metrics for the quantification of entrainment (winding number, correlation, peak detection), provides methods for the numerical exploration of the oscillator models (parameter scan, Arnold tongue, phase response curve), and optimization routines for fitting the models to the experimental data. This package aims to encompass all methods developed and used in the context of my Ph.D. thesis.

## Installation

The package can be installed by running
```julia
import Pkg
Pkg.add(url="https://github.com/vkumpost/OscillatorPopulation")
```

To make sure everything is ready to go we can run package tests
```julia
Pkg.test("OscillatorPopulation")
```

If the compilation fails due to `PyCall` package not being able to find a Python
installation, see the [installation guide for the PyPlot.jl package](https://github.com/JuliaPy/PyPlot.jl#installation).

## List of implemented methods

The complete list of exported functions can be obtain from Julia REPL by calling 
```julia
using OscillatorPopulation
?OscillatorPopulation
```

**Input and Output**
- `save_data`: Save dataframe as a csv file.
- `load_data`: Load dataframe from a csv file.
- `save_figure`: Save figure on the disk.

**Luminescence Data**
- `load_biolum`: Load bioluminescence data.
- `detect_events`: Detect events in the bioluminescence data.
- `biolum_zscore_traces`: Apply Z-score to bioluminescence traces.
- `biolum_mean`: Calculate mean of the bioluminescence traces.

**Events**
- `create_events`: Create events represented by a matrix.
- `create_events_cycle`: Ceate a cycle with specific lengh and period.
- `events_to_function`: Convert events matrix into a function.
- `plot_events`: Plot events.
- `create_callback`: Create a callback to apply events to a model.
- `detect_events`: Detect events in a binary vector.

**Exceptions**
- `OscillatorPopulationError`: A general exeption if something goes wrong.

**Metrics**
- `rsquared`: Coefficient of determination.
- `cmean`: Circular mean.
- `cstd`: Circular standard deviation.
- `window_xcorr`: Moving-window cross-correlation.
- `cxcorr`: Circular cross-correlation.
- `estimate_phase_array`: Estimate entrainment phase at each cycle.
- `estimate_phase_array_peaks`: Estimate entrainment phase using peak detection.
- `estimate_phase_array_cxcorr`: Estimate entrainment phase using circular cross-correlation.
- `estimate_order_parameter`: Estimate Kuramoto's order parameter.
- `estimate_period`: Estimate period.
- `estimate_winding_number`: Estimate winding number.
- `estimate_period_winding_number`: Estimate period based on the winding number.
- `create_simulation_function`: Generate a simulation function.

**Miscellaneous**
- `smooth`: Smooth a vector using moving average filter.
- `generate_random_values`: Generate a vector of random numbers.
- `benchmark`: Estimate evaluation time for a function.
- `find_closest`: Find the index of the closest value in an array.

**Model Library**
- `kfr`: Kim-Forger function.
- `load_model`: Load a model from the library.

**Model Functions**
- `set_initial_conditions!`: Set initial conditions.
- `set_timespan!`: Set time span for the simulation.
- `set_output!`: Set output function.
- `set_solver!`: Set numerical solver and its parameters.
- `get_parameter_index`: Get index of a parameter.
- `get_parameter`: Get value of a parameter.
- `set_parameter!`: Set value of a parameter.
- `set_input!`: Set input events.
- `print_info`: Print model information (initial conditions, parameter values, ...)
- `simulate_model`: Perform numerical simulation.

**Optimization**
- `estimate_initial_conditions`: Preestimate initial conditions for a model.
- `create_data_objective`: Create a cost function to fit a model to data.
- `create_entrainable_oscillator_objective`: Cost function to find parameters of an entrainmeble oscillator.
- `create_desynchronization_objective`: Cost function to find noise intensity based on the population-level damping rate.
- `optimize`: Minimize a cost function.
- `damped_sine`: Damped sine function.
- `polynomial`: Polynomial function.
- `fit_curve`: Fit curve to data.
- `binary_search`: Binary search.

**Population Solution**
- `simulate_population`: Simulate uncoupled population.
- `plot_solution`: Plot simulation output.
- `select_time`: Select specific time range from the simulation output.
- `select_subset`: Select a subset of trajectories from the simulation output.

**Parameter Scans**
- `scan`: Scan parameter values.
- `scan_arnold`: Scan input-signal parameters (amplitude, period, duty cycle).
- `plot_arnold`: Plot the output of `scan_arnold`.
- `select_arnold_row`: Select a specific row from the output of `scan_arnold`.
- `estimate_prc`: Estimate phase response curve.
- `estimate_T_prc`: Estimate T-cycle phase response curve.
