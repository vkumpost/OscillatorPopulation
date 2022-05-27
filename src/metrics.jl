"""
`estimate_phase_array`

Estimate entrainment phase. The input parameters are expected to represent a
model output entrained to a periodic square input. The first and last event are
discarted.

**Arguments**
- `t`: Time vector.
- `x`: Data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `normalize`: If `true`, the estimated phase is normalized by the event period.
- `smooth_span`: Length of the moving average filter (in samples) used to average
    `x` before the phase is estimated.
- `show_plots`: If `true`, a figure is generated for a visual control.

**Returns**
- `phase_arr`: Array of phases estimated at each cycle of the input data.
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
    phase_arr = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        idx = window_start .< locs .< window_end
        window_locs = locs[idx]
        window_pks = pks[idx]

        if isempty(window_pks)
            # If there are no peaks in the window, set phase to NaN
            push!(phase_arr, NaN)
        else
            # If there are peaks in the window, estimate phase
            idx_peak = argmax(window_pks)
            loc = window_locs[idx_peak]
            phase = loc - window_start
            if normalize
                phase /= (window_end - window_start)
            end
            push!(phase_arr, phase)
        end

        if show_plots
            ax.plot(loc, x[pr.indices[idx][idx_peak]], "o", color="red")
        end
    end

    return phase_arr

end


"""
`estimate_phase`

Estimate entrainment phase. The input parameters are expected to represent a
model output entrained to a periodic square input. The first and last event are
discarted.

**Arguments**
- `t`: Time vector.
- `x`: Data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `normalize`: If `true`, the estimated phase is normalized by the event period.
- `smooth_span`: Length of the moving average filter (in samples) used to average
    `x` before the phase is estimated.
- `show_plots`: If `true`, a figure is generated for a visual control.

**Returns**
- `phase`: Mean phase calculated over the events cycles.
- `phase_error`: Standard deviation for the mean phase.
"""
function estimate_phase(args...; kwargs...)

    phase_arr = estimate_phase_array(args...; kwargs...)

    phase = mean(phase_arr)
    phase_error = std(phase_arr)

    return phase, phase_error

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

Estimate the winding number number around the origin.

**Argument**
- `x`: Data vector for the x-coordinate.
- `y`: Data vector for the y-coordinate.

**Returns**
- `winding_number`: Estimated winding number.
"""
function estimate_winding_number(x, y)

    n = length(x)
    total_angle = 0
    for i = 1:n-1
        total_angle += atan(x[i]*y[i+1] - y[i]*x[i+1], x[i]*x[i+1] + y[i]*y[i+1])
    end
    winding_number = abs(total_angle / 2π)
    return winding_number

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
- `show_plots`: Show plots for visual verification (default `false`).
- `kwargs...`: Keyword arguments passed to `simulate_population`.

**Returns**
- `simulation_function`: Function that takes `Model` as a parameter and returns
    an array of metrics calculated from the numerical solution. If `Model` is
    not passed, the function returns the metric names as an array of strings.
"""
function create_simulation_function(property_names=nothing; transient=0.9,
    trajectories=1, variable=1, show_plots=false, kwargs...)

    if isnothing(property_names)
        property_names = ["minimum", "maximum", "phase", "phase_error"]
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
        events = solution.events

        if show_plots
            _, ax = subplots()
            ax.plot(t, x; color="black")
            plot_events(events)
            ax.set_title("Variable $variable")
        end

        # Initialize vector for the estimated properties
        property_values = fill(NaN, n_properties)

        # Variables to save properties, so they do not be estimated repeatedly
        phase = nothing
        phase_error = nothing

        # Iterate property names
        for (i_property, property_name) in enumerate(property_names)

            # Estimate property
            if property_name == "minimum"
                property_values[i_property] = minimum(x)

            elseif property_name == "maximum"
                property_values[i_property] = maximum(x)

            elseif property_name == "phase"
                if isnothing(phase)
                    phase, phase_error = estimate_phase(t, x, events; show_plots=show_plots)
                end
                property_values[i_property] = phase

            elseif property_name == "phase_error"
                if isnothing(phase_error)
                    phase, phase_error = estimate_phase(t, x, events; show_plots=show_plots)
                end
                property_values[i_property] = phase_error

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
