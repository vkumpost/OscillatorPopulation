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
`cxcorr`

Compute circular cross-correlation.

**Arguments**
- `x`: Array representing one period of a signal.
- `y`: Array representing one period of a signal.

**Keyword Arguments**
- `use_fft`: If `true`, the cross-correlation is calculated using the fast
    Fourier transform algorithm. Default value is `false`.

**Returns**
- `r`: Circular cross-correlation of `x` and `y`.
"""
function cxcorr(x, y; use_fft=false)

    if use_fft
        
        n = length(x)

        # p = plan_rfft(zeros(n); flags=FFTW.PATIENT)
        # x_fft = p*x
        # y_fft = p*y
        x_fft = rfft(x)
        y_fft = rfft(y)
        
        S = x_fft .* conj(y_fft)

        # nS = length(S)
        # ip = plan_irfft(zeros(ComplexF64, nS), n; flags=FFTW.PATIENT)
        # r = ip*S
        r = irfft(S, n)

    else
        n = length(x)
        r = fill(NaN, n)
        lags = 0:(n-1)
        for lag in lags
            x_circ = x[[(lag+1):end...; 1:lag...]]
            r[lag + 1] = sum(x_circ .* y)
        end
    end

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
    `"peak_prominence"` (default), `"cxcorr"`.
- `smooth_span`: Length of the moving average filter (in samples) used to average
    `x` before the phase is estimated.
- `show_plots`: [NOT IMPLEMENTED]

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data. If 
    `x` is a matrix, `phase_array` is also a matrix, where each column contains
    estimated phases for the individual data vectors.
"""
function estimate_phase_array(t, x, events; method="peak_prominence",
    smooth_span=1, show_plots=false)

    if smooth_span > 1
        x = smooth(x, span=smooth_span)
    end

    if method == "peak_prominence"
        phase_array = estimate_phase_array_peak_prominence(t, x, events)
    elseif method == "cxcorr"
        phase_array = estimate_phase_array_cxcorr(t, x, events)
    end

    return phase_array

end


"""
`estimate_phase_array_peak_prominence`

Estimate entrainment phase based on the peak prominence.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data.
"""
function estimate_phase_array_peak_prominence(t, x, events)

    pr = findpeaks(x, t)
    pks = peakprominences(pr)
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
- `period`: Estimated period (the position of the highest).
- `peak`: Height of the highest peak.
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
    pr = findpeaks(R, (0:(n-1)) .* mean(diff(t)), sortstr="descend")

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
`estimate_winding_number_period`

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
function estimate_winding_number_period(x, y, time_duration; remove_mean=true)

    winding_number = estimate_winding_number(x, y; remove_mean=remove_mean)
    period = time_duration / winding_number
    return period

end


function _estimate_entrainment_properties()
    property_names = ["minimum", "maximum", "amplitude", "rms",
            "winding_number", "phase_coherence", "mean_phase",
            "phase_coherence_population", "collective_phase"]
    return property_names
end


function _estimate_entrainment_properties(solution; variable_x=1, variable_y=2, property_names=nothing)

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

    # Initialize vector for the estimated properties
    property_values = fill(NaN, n_properties)

    # Variables to save properties, so they do not be estimated repeatedly
    phase_coherence = nothing
    mean_phase = nothing
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
            input_period = events[3, 1] - events[2, 1]
            winding_number_period = estimate_winding_number_period(x, y, time_duration)
            property_values[i_property] = input_period / winding_number_period

        elseif property_name == "phase_coherence"
            if isnothing(phase_coherence)
                phase_array = estimate_phase_array(t, x, events; method="cxcorr")
                phase_coherence, mean_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = phase_coherence
            
        elseif property_name == "mean_phase"
            if isnothing(mean_phase)
                phase_array = estimate_phase_array(t, x, events; method="cxcorr")
                phase_coherence, mean_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = mean_phase

        elseif property_name == "phase_coherence_population"
            if isnothing(phase_coherence_population)
                phase_array = estimate_phase_array(t, U, events; method="cxcorr")
                phase_coherence_population, collective_phase = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = phase_coherence_population

        elseif property_name == "collective_phase"    
            if isnothing(collective_phase)
                phase_array = estimate_phase_array(t, U, events; method="cxcorr")
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
- `variable_x`: Trajectory that is used to calculate the metrics (default 1).
- `variable_y`: The second variable for the phase plane (default 2).
- `single_cells`: If `true`, parameters for single cells are also estimated. 
    Default value is `false`.
- `show_plots`: Show plots for visual verification (default `false`).
- `kwargs...`: Keyword arguments passed to `simulate_population`.

**Returns**
- `simulation_function`: Function that takes `Model` as a parameter and returns
    an array of metrics calculated from the numerical solution. If `Model` is
    not passed, the function returns the metric names as an array of strings.
"""
function create_simulation_function(property_names=nothing; transient=0.9,
    trajectories=1, variable_x=1, variable_y=2, single_cells=false,
    show_plots=false, kwargs...)

    simulation_function = function (model=nothing; show_plots=show_plots)

        # Return property names, if the model was not passed
        if isnothing(model)
            property_names = _estimate_entrainment_properties()
            property_names_output = deepcopy(property_names)
            if single_cells
                for i in 1:trajectories
                    extended_names = [x * "_$(i)" for x in property_names]
                    append!(property_names_output, extended_names)    
                end
            end
            return property_names_output
        end

        # Simulate the population
        solution = simulate_population(model, trajectories; kwargs...)
        min_time = transient * maximum(solution.time)
        solution = select_time(solution, min_time=min_time)

        if show_plots
            _, ax = subplots()
            t = solution.time
            x = solution.mean[:, variable_x]
            events = solution.events
            ax.plot(t, x; color="black")
            plot_events(events)
            ax.set_title("Mean")
        end

        property_values = _estimate_entrainment_properties(solution;
            variable_x=variable_x, variable_y=variable_y,
            property_names=property_names)
        if single_cells
            for i in 1:trajectories
                solution_subset = select_subset(solution, i)
                property_values_subset = _estimate_entrainment_properties(
                    solution_subset;
                    variable_x=variable_x, variable_y=variable_y,
                    property_names=property_names)
                append!(property_values, property_values_subset)
            end
        end
        return property_values

    end

    return simulation_function

end
