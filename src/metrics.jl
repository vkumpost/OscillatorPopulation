"""
`rsquared`

Calculate the coefficient of determination.

**Arguments**
- `x`: Vector of real data.
- `y`: Vector of predicted values.

**Returns**
- `R`: The coefficient of determination as a number on the interval `[-∞, 1]`.

**Reference**
- https://doi.org/10.1016/0022-1694(70)90255-6
"""
function rsquared(x, y)

    μ = mean(x)
    SSres = sum((x .- y).^2)
    SStot = sum((x .- μ).^2)
    R = 1 - SSres/SStot
    return R

end


"""
`cmean`

Compute circular mean of an array.

**Arguments**
- `x`: Array of angles normalized to the interval [0, 1].

**Returns**
- `m`: Circular mean of `x`.
"""
function cmean(x)
    x_rad = x * 2π
    S = mean(sin.(x_rad))
    C = mean(cos.(x_rad))
    m_rad = atan(S, C) 
    m = m_rad / 2π
    return m
end


"""
`xcorr`

Compute cross-correlation. This is a particular implementation of the
cross-correlation function that assumes that `x` is shorter than `y`. `x` is
moved along `y` and the cross-correlation function is calculated as weighted
moving average, where `x` are the weights and `y` is the averaged signal.

**Arguments**
- `x`: A signal array.
- `y`: A signal array.

**Returns**
- `r`: Cross-correlation of `x` and `y`.
"""
function xcorr(x, y)

    nx = length(x)
    ny = length(y)

    if nx > ny
        msg = "x must be shorter than y!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    nr = ny - nx + 1

    r = fill(NaN, nr)
    lags = 0:(nr-1)
    for lag in lags
        y_window = y[(lag+1):(lag+nx)]
        r[lag + 1] = sum(x .* y_window)
    end

    return r

end


"""
`xcorr`

Compute autocorrelation.

**Arguments**
- `x`: A signal array.

**Keyword Arguments**
- `window_length`: Length of the window taken from the beginning of `x` and used
    as a moving vector to calculate the cross-correlation (see `xcorr(x, y)`).
    The length is given as fraction of the length of `x`. Default value is 0.5.

**Returns**
- `r`: Autocorrelation of `x`.
"""
function xcorr(x; window_length=0.5)

    nx = length(x)
    n_window = round(Int64, window_length * nx)
    x_window = x[1:n_window]
    r = xcorr(x_window, x)
    return r

end


"""
`cxcorr`

Compute circular cross-correlation.

**Arguments**
- `x`: Array representing one period of a signal.
- `y`: Array representing one period of a signal.

**Returns**
- `r`: Circular cross-correlation of `x` and `y`. Normalized so
    the maximum of `r` is 1.
"""
function cxcorr(x, y)

    nx = length(x)
    xx = vcat(x, fill(0, nx))
    yy = vcat(y, y)
    r = crosscov(xx, yy, 0:(nx-1); demean=false)
    r ./= maximum(r)

    return r

end

"""
`estimate_phase_array`

Estimate entrainment phase. The input parameters are expected to represent a
model output entrained to a periodic square input. The first and last event are
discarted.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `method`: Choose method to estimate the phase. Possible options are 
    `"peak_height"` (default), `"peak_prominence"`, and `"cxcorr"`.
- `smooth_span`: Length of the moving average filter (in samples) used to average
    `x` before the phase is estimated. Default value is 1 (no smoothing).
- `show_plots`: [NOT IMPLEMENTED]

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data. If 
    `x` is a matrix, `phase_array` is also a matrix, where each column contains
    estimated phases for the individual data vectors.
"""
function estimate_phase_array(t, x, events; method="peak_height",
    smooth_span=1, show_plots=false)

    if smooth_span > 1
        x = smooth(x, span=smooth_span)
    end

    if method == "peak_height"
        phase_array = estimate_phase_array_peaks(t, x, events; use_prominence=false)
    elseif method == "peak_prominence"
        phase_array = estimate_phase_array_peaks(t, x, events; use_prominence=true)
    elseif method == "cxcorr"
        phase_array = estimate_phase_array_cxcorr(t, x, events)
    end

    return phase_array

end


"""
`estimate_phase_array_peaks`

Estimate entrainment phase based on the peak prominence.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `use_prominence`: If `true`, the algorithm will use peak prominences
    instead of peak heights to determine the heights peaks. Default value is
    `false`.

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data.
"""
function estimate_phase_array_peaks(t, x, events; use_prominence=false)

    pr = findpeaks(x, t)
    if use_prominence
        pks = peakprominences(pr)
    else
        pks = peakheights(pr)
    end
    locs = peaklocations(pr)

    n_events = size(events, 1)
    phase_array = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        idx = window_start .<= locs .< window_end
        window_locs = locs[idx]
        window_pks = pks[idx]

        if isempty(window_pks)
            # If there are no peaks in the window, set phase to NaN
            push!(phase_array, NaN)
        else
            # If there are peaks in the window, estimate phase
            idx_peak = argmax(window_pks)
            loc = window_locs[idx_peak]
            phase = loc - window_start
            phase /= (window_end - window_start)
            push!(phase_array, phase)
        end
    end

    return phase_array

end


"""
`estimate_phase_array_cxcorr`

Estimate entrainment phase based on the circular cross-correlation.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data.
"""
function estimate_phase_array_cxcorr(t, x, events)

    fun = events_to_function(events)
    y = fun.(t)

    n_events = size(events, 1)
    phase_array = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        window_idices = window_start .<= t .< window_end
        if sum(window_idices) == 0
            push!(phase_array, NaN)
        else
            r = cxcorr(x[window_idices], y[window_idices])
            _, max_index = findmax(r)
            phase = (max_index - 1) / sum(window_idices)
            push!(phase_array, phase)
        end
    end

    return phase_array

end


function estimate_phase_array(t, X::Matrix, events; kwargs...)

    # Call estimate_phase_array on each column of a matrix
    n_trajectories = size(X, 2)
    n_events = size(events, 1) - 2
    phase_matrix = fill(NaN, n_events, n_trajectories)
    for i_trajectory = 1:n_trajectories
        x = X[:, i_trajectory]
        phase_array = estimate_phase_array(t, x, events; kwargs...)
        phase_matrix[:, i_trajectory] = phase_array
    end
    return phase_matrix

end


"""
`estimate_order_parameter`

Estimate period of a signal based on its autocorrelation function.

**Argument**
- `phase_array`: Array of phases normalized to interval [0, 1].

**Keyword Arguments**
- `skip_nan`: Skip NaN values. If `false`, `NaN` values are propageted to the
    output. Default value is `true`.

**Returns**
- `phase_coherence`: A number between 0 (no coherence) to 1 (full coherence).
    If `phase_array` is a matrix, phase coherence is calculated as a mean of
    the phase coherences calculated from individual rows.
- `collective_phase`: A number between 0 and 1 indicating the average phase.
    If `phase_array` is a matrix, collective phase is circular mean of the
    indivudal phases.
"""
function estimate_order_parameter(phase_array; skip_nan=true)

    if minimum(phase_array) < 0 || maximum(phase_array) > 1
        msg = "Phase must be on interval [0, 1]!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    nan_indices = isnan.(phase_array)
    if !skip_nan && any(nan_indices)
        return NaN, NaN
    else
        phase_array = phase_array[.!nan_indices]
        if isempty(phase_array)
            return NaN, NaN
        end
    end

    n = length(phase_array)
    complex_order_parameter = 0im
    for phase in phase_array
        complex_order_parameter += exp(2π * phase * im)
    end
    complex_order_parameter /= n
    phase_coherence = abs(complex_order_parameter)
    collective_phase = angle(complex_order_parameter)  # ∈ [-π, π]
    if collective_phase < 0  # map from [-π, 0] on [π, 2π]
        collective_phase += 2π
    end
    collective_phase /= 2π  # map from [0, 2π] on [0, 1]
    return phase_coherence, collective_phase

end


function estimate_order_parameter(phase_matrix::Matrix)

    # Call estimate_order_parameter on each row of a matrix
    n = size(phase_matrix, 1)
    phase_coherence_arr = fill(NaN, n)
    collective_phase_arr = fill(NaN, n)
    for i = 1:n
        phase_array = phase_matrix[i, :]
        phase_coherence, collective_phase = estimate_order_parameter(phase_array)
        phase_coherence_arr[i] = phase_coherence
        collective_phase_arr[i] = collective_phase
    end
    phase_coherence_mean = mean(phase_coherence_arr)
    collective_phase_mean = cmean(collective_phase_arr)
    return phase_coherence_mean, collective_phase_mean

end


"""
`estimate_period(t, x; kwargs...)`

Estimate period of a signal based on its autocorrelation function.

**Argument**
- `t`: Time vector.
- `x`: Data vector.

**Keyword Arguments**
- `n_lags`: Number of lags used for autocorrelation. Default is `length(x)`.
- `show_plots`: Visualize the result. Default is false.

**Returns**
- `period`: Estimated period (the position of the first peak).
- `peak`: Height of the first peak.
"""
function estimate_period(t, x; n_lags=length(x), show_plots=false)

    # If x is constat, return NaNs
    if all(x[1] .== x)
        return (NaN, NaN)
    end

    # Normalize input
    x = zscore(x)

    # Calculate autocorrelation
    n = min(n_lags, length(x))  # number of lags to calculate
    R = autocor(x, 0:(n-1))  # autocorrelation function

    # Find peaks of the autocorrelation function
    pr = findpeaks(R, (0:(n-1)) .* mean(diff(t)))  # , sortstr="descend", sortref="prominence")

    # If there are no peaks, return NaNs
    if length(pr) == 0
        return (NaN, NaN)
    end

    # Get the highest peaks (peak height and location = period)
    peak = peakheights(pr)[1]
    period = peaklocations(pr)[1]

    # Visualize the results
    if show_plots
        fig, axs = subplots(2)
        axs[1].plot(t, x, color="black")
        axs[1].set_title("Time Series")
        axs[2].plot((0:(n-1)) .* mean(diff(t)), R, color="black")
        axs[2].plot([period, period], [0, peak], color="red")
        axs[2].set_title("peak = $(round(peak, digits=2)), " *
            "period = $(round(period, digits=2))")
        fig.tight_layout()
    end

    return (period, peak)

end


"""
`estimate_winding_number`

Estimate the winding number around the origin.

**Argument**
- `x`: Data vector for the x-coordinate.
- `y`: Data vector for the y-coordinate.

**Keyword Argument**
- `remove_mean`: Subtract mean from `x` and `y` before estimating the winding
    number. Default value is `true`.

**Returns**
- `winding_number`: Estimated winding number.
"""
function estimate_winding_number(x, y; remove_mean=true)

    if remove_mean
        x = x .- mean(x)
        y = y .- mean(y)
    end
    n = length(x)
    total_angle = 0
    for i = 1:n-1
        total_angle += atan(x[i]*y[i+1] - y[i]*x[i+1], x[i]*x[i+1] + y[i]*y[i+1])
    end
    winding_number = abs(total_angle / 2π)
    return winding_number

end


"""
`estimate_period_winding_number`

Estimate period based on the winding number.

**Argument**
- `x`: Data vector for the x-coordinate.
- `y`: Data vector for the y-coordinate.
- `time_duration`: Time duration of `x` and `y`.

**Keyword Argument**
- `remove_mean`: Subtract mean from `x` and `y` before estimating the winding
    number. Default value is `true`.

**Returns**
- `period`: Estimated period.
"""
function estimate_period_winding_number(x, y, time_duration; remove_mean=true)

    winding_number = estimate_winding_number(x, y; remove_mean=remove_mean)
    period = time_duration / winding_number
    return period

end


function _estimate_entrainment_properties()
    property_names = ["minimum", "maximum", "amplitude", "rms",
            "winding_number", "autocorrelation",
            "phase_coherence", "mean_phase",
            "phase_coherence_cxcorr", "mean_phase_cxcorr",
            "phase_coherence_population", "collective_phase"]
    return property_names
end


function _estimate_entrainment_properties(solution; smooth_span=1, variable_x=1, variable_y=2, property_names=nothing)

    if isnothing(property_names)
        property_names = _estimate_entrainment_properties()
    end
    n_properties = length(property_names)

    # Extract time, state, and events
    t = solution.time
    x = solution.mean[:, variable_x]
    y = solution.mean[:, variable_y]
    U = solution.trajectories[:, variable_x, :]
    events = solution.events

    if smooth_span > 1
        x = smooth(x, span=smooth_span)
        y = smooth(y, span=smooth_span)
    end

    # Initialize vector for the estimated properties
    property_values = fill(NaN, n_properties)

    # Variables to save properties, so they do not be estimated repeatedly
    phase_coherence = nothing
    mean_phase = nothing
    phase_coherence_cxcorr = nothing
    mean_phase_cxcorr = nothing
    phase_coherence_population = nothing
    collective_phase = nothing

    # Iterate property names
    for (i_property, property_name) in enumerate(property_names)

        # Estimate property
        if property_name == "minimum"
            property_values[i_property] = minimum(x)

        elseif property_name == "maximum"
            property_values[i_property] = maximum(x)

        elseif property_name == "amplitude"
            property_values[i_property] = maximum(x) .- minimum(x)

        elseif property_name == "rms"
            property_values[i_property] = sqrt(mean((x .- mean(x)) .^ 2))

        elseif property_name == "winding_number"
            time_duration = maximum(t) - minimum(t)
            winding_number_period = estimate_period_winding_number(x, y, time_duration)
            input_period = events[3, 1] - events[2, 1]
            property_values[i_property] = winding_number_period / input_period

        elseif property_name == "autocorrelation"
            autocorrelation_period, _ = estimate_period(t, x)
            input_period = events[3, 1] - events[2, 1]
            property_values[i_property] = autocorrelation_period / input_period

        elseif property_name == "phase_coherence"
            if isnothing(phase_coherence)
                phase_array = estimate_phase_array(t, x, events; method="peak_height")
                phase_coherence, mean_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = phase_coherence
            
        elseif property_name == "mean_phase"
            if isnothing(mean_phase)
                phase_array = estimate_phase_array(t, x, events; method="peak_height")
                phase_coherence, mean_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = mean_phase

        elseif property_name == "phase_coherence_cxcorr"
            if isnothing(phase_coherence_cxcorr)
                phase_array = estimate_phase_array(t, x, events; method="cxcorr")
                phase_coherence_cxcorr, mean_phase_cxcorr = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = phase_coherence_cxcorr

        elseif property_name == "mean_phase_cxcorr"
            if isnothing(mean_phase_cxcorr)
                phase_array = estimate_phase_array(t, x, events; method="cxcorr")
                phase_coherence_cxcorr, mean_phase_cxcorr = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = mean_phase_cxcorr

        elseif property_name == "phase_coherence_population"
            if isnothing(phase_coherence_population)
                phase_array = estimate_phase_array(t, U, events; method="peak_height", smooth_span=smooth_span)
                phase_coherence_population, collective_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = phase_coherence_population

        elseif property_name == "collective_phase"    
            if isnothing(collective_phase)
                phase_array = estimate_phase_array(t, U, events; method="peak_height", smooth_span=smooth_span)
                phase_coherence_population, collective_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = collective_phase

        else
            msg = "Property `$(property_name)` is not valid!"
            err = OscillatorPopulationError(msg)
            throw(err)
            
        end
    end

    return property_values

end


"""
`create_simulation_function`

Generate a function that simulates a model population and apply metrics to the
    numerical solution. 

**Arguments**
- `property_names`: Names of metrics as an array of strings.

**Keyword Arguments**
- `transient`: The proportion of the solution to be discrated to avoid transient
    effects (default 0.9).
- `trajectories`: Number of trajectories in the population (default 1).
- `cycle_samples`: The number of samples within one cycle (default 10).
- `smooth_span`: Smooth the time series using running average of length
    `smooth_span`. Default value is 1 (no smoothing).
- `variable_x`: Trajectory that is used to calculate the metrics (default 1).
- `variable_y`: The second variable for the phase plane (default 2).
- `single_cells`: If `true`, parameters for single cells are also estimated. 
    Default value is `false`.
- `subpopulations`: An array indicating the size of subpoopulations, for which
    the parameters are also estimated. For example `[1, 10, 100]` will calculate
    the parameter values also for populations of size 1, 10, and 100. Default
    value is `Int64[]`.
- `average_subpopulations`: If `true` (default), the subpopulations are averaged
    over all subpopulations that can fit the main population. For example
    subpopulation of size `100` can be extracted 10-times from the main
    population of 1000 without overlaps.
- `show_plots`: Show plots for visual verification (default `false`).
- `kwargs...`: Keyword arguments passed to `simulate_population`.

**Returns**
- `simulation_function`: Function that takes `Model` as a parameter and returns
    an array of metrics calculated from the numerical solution. If `Model` is
    not passed, the function returns the metric names as an array of strings.
"""
function create_simulation_function(property_names=nothing; transient=0.9,
    trajectories=1, cycle_samples=10, smooth_span=1, variable_x=1, variable_y=2,
    single_cells=false, subpopulations=Int64[], average_subpopulations=true,
    show_plots=false, kwargs...)

    simulation_function = function (model=nothing; show_plots=show_plots)

        # Return property names, if the model was not passed
        if isnothing(model)
            if isnothing(property_names)
                property_names = _estimate_entrainment_properties()
            end
            property_names_output = deepcopy(property_names)
            if single_cells
                for i in 1:trajectories
                    extended_names = [x * "_$(i)" for x in property_names]
                    append!(property_names_output, extended_names)    
                end
            end
            for subpopulation in subpopulations
                extended_names = [x * "_n$(subpopulation)" for x in property_names]
                append!(property_names_output, extended_names)    
            end
            return property_names_output
        end

        # Simulate the population
        max_time = model.problem.tspan[end]
        time_step = diff(model.input[1][2:3, 1])[1] / cycle_samples
        min_time = transient * max_time
        _, idx = findmin(abs.(min_time .- model.input[1][:, 1]))
        min_time = model.input[1][idx, 1]
        model2 = deepcopy(model)
        set_solver!(model2, saveat=min_time:time_step:max_time)
        solution = simulate_population(model2, trajectories; kwargs...)
        solution = select_time(solution, min_time=min_time)

        # solution = simulate_population(model, trajectories; kwargs...)
        # min_time = transient * maximum(solution.time)
        # solution = select_time(solution, min_time=min_time)

        if show_plots
            _, ax = subplots()
            t = solution.time
            x = solution.mean[:, variable_x]
            events = solution.events
            ax.plot(t, x; color="black")
            plot_events(events)
            ax.set_title("Mean")
        end

        property_values = _estimate_entrainment_properties(
            solution;
            smooth_span=smooth_span,
            variable_x=variable_x,
            variable_y=variable_y,
            property_names=property_names
        )

        if single_cells
            for i in 1:trajectories
                solution_subset = select_subset(solution, i)
                property_values_subset = _estimate_entrainment_properties(
                    solution_subset;
                    smooth_span=smooth_span,
                    variable_x=variable_x,
                    variable_y=variable_y,
                    property_names=property_names
                )
                append!(property_values, property_values_subset)
            end
        end

        for subpopulation in subpopulations

            property_values_subset_mean = nothing
            counter = 0
            if average_subpopulations
                sub_start = 1
                sub_end = subpopulation
            else
                sub_start = trajectories - subpopulation + 1
                sub_end = trajectories
            end
            while sub_end <= trajectories
                solution_subset = select_subset(solution, sub_start:sub_end)
                property_values_subset = _estimate_entrainment_properties(
                    solution_subset;
                    smooth_span=smooth_span,
                    variable_x=variable_x,
                    variable_y=variable_y,
                    property_names=property_names
                )
                if isnothing(property_values_subset_mean)
                    property_values_subset_mean = property_values_subset
                else
                    property_values_subset_mean .+= property_values_subset
                end
                counter += 1
                sub_start = sub_end + 1
                sub_end = sub_start + subpopulation - 1
            end
            append!(property_values, property_values_subset_mean ./ counter)
        end
            
        return property_values

    end

    return simulation_function

end
