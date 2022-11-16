function _estimate_entrainment_properties()
    property_names = ["minimum", "maximum", "amplitude", "rms",
            "winding_number", "autocorrelation",
            "phase_coherence", "mean_phase",
            "phase_coherence_cxcorr", "mean_phase_cxcorr",
            "phase_coherence_population", "collective_phase",
            "phase_coherence_population_cxcorr", "collective_phase_cxcorr"]
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
    phase_coherence_population_cxcorr = nothing
    collective_phase_cxcorr = nothing

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

        elseif property_name == "phase_coherence_population_cxcorr"
            if isnothing(phase_coherence_population_cxcorr)
                phase_array = estimate_phase_array(t, U, events; method="cxcorr", smooth_span=smooth_span)
                phase_coherence_population_cxcorr, collective_phase_cxcorr = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = phase_coherence_population_cxcorr

        elseif property_name == "collective_phase_cxcorr"    
            if isnothing(collective_phase_cxcorr)
                phase_array = estimate_phase_array(t, U, events; method="cxcorr", smooth_span=smooth_span)
                phase_coherence_population_cxcorr, collective_phase_cxcorr = estimate_order_parameter(phase_array)
            end
            property_values[i_property] = collective_phase_cxcorr

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
        if model.problem isa JumpProblem
            max_time = model.problem.prob.tspan[2]
        else
            max_time = model.problem.tspan[2]
        end
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
