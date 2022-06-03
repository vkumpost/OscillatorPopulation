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
`estimate_phase_array`

Estimate entrainment phase. The input parameters are expected to represent a
model output entrained to a periodic square input. The first and last event are
discarted.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `normalize`: If `true`, the estimated phase is normalized by the event period.
- `smooth_span`: Length of the moving average filter (in samples) used to average
    `x` before the phase is estimated.
- `show_plots`: If `true`, a figure is generated for a visual control.

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data. If 
    `x` is a matrix, `phase_array` is also a matrix, where each column contains
    estimated phases for the individual data vectors.
"""
function estimate_phase_array(t, x, events; normalize=true, smooth_span=1, show_plots=false)

    if smooth_span > 1
        x = smooth(x, span=smooth_span)
    end
    pr = findpeaks(x, t)
    pks = peakprominences(pr)
    locs = peaklocations(pr)

    if show_plots
        fig, ax = subplots()
        ax.plot(t, x, color="black")
        ax.plot(t[pr.indices], x[pr.indices], "o", color="blue")
        plot_events(events, ax=ax)
    end

    n_events = size(events, 1)
    phase_array = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        idx = window_start .< locs .< window_end
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
            if normalize
                phase /= (window_end - window_start)
            end
            push!(phase_array, phase)
        end

        if show_plots
            ax.plot(loc, x[pr.indices[idx][idx_peak]], "o", color="red")
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

**Returns**
- `phase_coherence`: A number between 0 (no coherence) to 1 (full coherence).
    If `phase_array` is a matrix, phase coherence is calculated as a mean of
    the phase coherences calculated from individual rows.
- `collective_phase`: A number between 0 and 1 indicating the average phase.
    If `phase_array` is a matrix, collective phase is circular mean of the
    indivudal phases.
"""
function estimate_order_parameter(phase_array)

    if minimum(phase_array) < 0 || maximum(phase_array) > 1
        msg = "Phase must be on interval [0, 1]!"
        err = OscillatorPopulationError(msg)
        throw(err)
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
- `variable`: Trajectory that is used to calculate the metrics (default 1).
- `variable_2`: The second variable for the phase plane (default 2).
- `show_plots`: Show plots for visual verification (default `false`).
- `kwargs...`: Keyword arguments passed to `simulate_population`.

**Returns**
- `simulation_function`: Function that takes `Model` as a parameter and returns
    an array of metrics calculated from the numerical solution. If `Model` is
    not passed, the function returns the metric names as an array of strings.
"""
function create_simulation_function(property_names=nothing; transient=0.9,
    trajectories=1, variable=1, variable_2=2, show_plots=false, kwargs...)

    if isnothing(property_names)
        property_names = ["minimum", "maximum", "amplitude", "winding_number",
            "phase_coherence", "collective_phase"]
    end
    n_properties = length(property_names)

    simulation_function = function (model=nothing; show_plots=show_plots)

        # Return property names, if the model was not passed
        if isnothing(model)
            return property_names
        end

        # Simulate the population
        solution = simulate_population(model, trajectories; kwargs...)
        min_time = transient * maximum(solution.time)
        solution = select_time(solution, min_time=min_time)

        # Extract time, state, and events
        t = solution.time
        x = solution.mean[:, variable]
        y = solution.mean[:, variable_2]
        events = solution.events
        time_duration = maximum(t) - minimum(t)
        input_period = events[3, 1] - events[2, 1]

        if show_plots
            _, ax = subplots()
            ax.plot(t, x; color="black")
            plot_events(events)
            ax.set_title("Variable $variable")
        end

        # Initialize vector for the estimated properties
        property_values = fill(NaN, n_properties)

        # Variables to save properties, so they do not be estimated repeatedly
        phase_coherence = nothing
        collective_phase = nothing

        # Iterate property names
        for (i_property, property_name) in enumerate(property_names)

            # Estimate property
            if property_name == "minimum"
                property_values[i_property] = minimum(x)

            elseif property_name == "maximum"
                property_values[i_property] = maximum(x)

            elseif property_name == "amplitude"
                property_values[i_property] = maximum(x) - minimum(x)

            elseif property_name == "winding_number"
                winding_number_period = estimate_winding_number_period(x, y, time_duration)
                property_values[i_property] = input_period / winding_number_period

            elseif property_name == "phase_coherence"
                if isnothing(phase_coherence)
                    phase_array = estimate_phase_array(t, x, events)
                    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
                end
                property_values[i_property] = phase_coherence
                
            elseif property_name == "collective_phase"
                if isnothing(collective_phase)
                    phase_array = estimate_phase_array(t, x, events)
                    phase_coherence, collective_phase = estimate_order_parameter(phase_array)
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

    return simulation_function

end
