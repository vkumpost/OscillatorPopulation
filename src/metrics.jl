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
- `show_plots`: If `true`, a figure is generated for a visual control.

**Returns**
- `phase`: Mean phase calculated over the events cycles.
- `phase_error`: Standard deviation for the mean phase.
"""
function estimate_phase(t, x, events; normalize=true, show_plots=false)

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
            push!(phase_arr, NaN)
            break
        end
        idx_peak = argmax(window_pks)
        loc = window_locs[idx_peak]
        phase = loc - window_start
        if normalize
            phase /= (window_end - window_start)
        end
        push!(phase_arr, phase)

        if show_plots
            ax.plot(loc, x[pr.indices[idx][idx_peak]], "o", color="red")
        end
    end

    phase = mean(phase_arr)
    phase_error = std(phase_arr)

    return phase, phase_error

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
