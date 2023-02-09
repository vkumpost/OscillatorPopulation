"""
`_find_all_combinations`

Find all possible combinations among the vectors of values. 

**Arguments**
- `vectors`: Vector of vectors representing different value sets.

**Returns**
- `combinations`: Matrix containing all combinations of values from the
    individual vectors. Each row represents one combination.

**Example**
```
julia> vectors = [[1, 2], [3, 4, 5], [6]];
julia> combinations = _find_all_combinations(vectors);
julia> print(combinations);
[1 3 6; 1 4 6; 1 5 6; 2 3 6; 2 4 6; 2 5 6]
```
"""
function _find_all_combinations(vectors)

    # Number of input vectors
    n = length(vectors)  

    # The last column of the combinations matrix is the last input vector
    combinations = hcat(vectors[n])

    # Iterate the remaining vectors
    for i = (n-1):-1:1

        # Get the current vector and its length
        x = vectors[i]
        nx = length(x)

        # Get the current size of the combinations matrix
        nc = size(combinations, 1)

        # Repeat the current combinations for each element in x 
        combinations = repeat(combinations, nx)

        # Add an element from x to each repetition of combinations
        x = vcat([fill(xi, nc) for xi in x]...)
        combinations = hcat(x, combinations)

    end

    return combinations

end


"""
`_binary_boundary_search`

Estimate a boundary value `X` at which `fun(X)` changes its value from `0` to `1`.
That is, `fun(x) = 0` for all `x < X` and `fun(x) = 1` for all `x > X`.

**Arguments**
- `fun`: Function that returns `0` if `x < X` and `1` if `x > X`.
- `x_range`: Search range as `[x_min, x_max]`.
- `err`: Maximal deviation of the result from the real value.

**Returns**
- `X`: Boundary value, at which `fun(x)` changes from `0` to `1`.
"""
function _binary_boundary_search(fun, x_range, err)

    x1 = x_range[1]
    x2 = x_range[2]

    y1 = fun(x1)
    if y1 == 1
        return x1
    end

    y2 = fun(x2)
    if y2 == 0
        return x2
    end

    while abs(x1 - x2) > 2err

        t0 = (x1 + x2) / 2
        x0 = fun(t0)

        if x0 == y1
            y1 = x0
            x1 = t0
        else
            y2 = x0
            x2 = t0
        end

    end

    X = (x1 + x2) / 2
    return X

end


"""
`scan`

Scan parameters of a model.

**Arguments**
- `model`: Model.
- `parameters`: A vector of pairs. Pairs indicate parameter names and values.
- `simulation_function`: A function that takes in `Model`, and returns
    some descriptive parameters. If the function is called without any arguments,
    it returns the names of the parameters as a vector of strings. For example:
    ```
    julia> simulation_function(model)
    [1.1, 2.2, 3.3]
    julia> simulation_function()
    ["Period", "Amplitude", "Phase"]
    ```

**Keyword Arguments**
- `catch_errors`: If `true`, errors from the `simulation_function` will be catched
    and parameters set to `NaN` for the corresponding entry.
- `show_progress`: If `true`, show a progress bar in the terminal.

**Returns**
- `df`: Dataframe with results. First columns correspond to the input parameters
    followed by columns returned by `simulation_function` for the corresponding
    parameters.
"""
function scan(model::Model, parameters::Vector, simulation_function::Function;
    catch_errors=true, show_progress::Bool=false)

    # Get values from the parameter dictionary
    parameter_names = [x[1] for x in parameters]
    parameter_values = [x[2] for x in parameters]

    # Find all possible parameter combinations
    parameter_value_combinations = _find_all_combinations(parameter_values)
    n = size(parameter_value_combinations, 1)
    
    # Initialize summary matrix
    summary_names = simulation_function()
    scan_summary = fill(NaN, n, length(summary_names))

    # Initialize progress bar
    if show_progress
        progressmeter = ProgressMeter.Progress(n; barlen=20)
    end

    # Iterate all parameter value combinations
    lk = ReentrantLock()
    Threads.@threads for i = 1:n

        # Protect the original model from overwriting
        model2 = deepcopy(model)

        # Set model parameters
        set_parameter!(model2, parameter_names, parameter_value_combinations[i, :])

        # Calculate parameters for the simulation
        try
            scan_summary[i, :] = simulation_function(model2)
        catch err
            @warn "An error occured for $(parameter_names) = $(parameter_values)"
            if !catch_errors
                throw(err)
            else
                scan_summary[i, :] .= NaN
            end
        end

        # Update progress bar
        if show_progress
            lock(lk) do
                ProgressMeter.next!(progressmeter)
            end
        end

    end

    # Build the output dataframe
    matrix = hcat(parameter_value_combinations, scan_summary)
    names = vcat(parameter_names, summary_names)
    df = DataFrame(matrix, names)

    return df

end


"""
`scan_arnold`

Estimate arnold tongue and/or onion.

**Arguments**
- `model`: Model.
- `simulation_function`: A function that takes in `Model`, and returns
    some descriptive parameters. If the function is called without any arguments,
    it returns the names of the parameters as a vector of strings. For example:
    ```
    julia> simulation_function(model)
    [1.1, 2.2, 3.3]
    julia> simulation_function()
    ["Period", "Amplitude", "Phase"]
    ```

**Keyword Arguments**
- `input_amplitudes`: Input amplitudes to scan.
- `input_periods`: Input periods to scan.
- `input_duty_cycles`: Input duty cycles to scan.
- `input_parameter`: Name of the parameter that controls the input amplitude.
- `catch_errors`: If `true`, errors from the `simulation_function` will be catched
    and parameters set to `NaN` for the corresponding entry.
- `show_progress`: If `true`, show progress in the terminal.

**Returns**
- `df`: Dataframe with results. First columns correspond to the input amplitude,
    period, and duty cycle. The following columns correspond to the values
    returned by `simulation_function` for the corresponding parameters.
"""
function scan_arnold(model, simulation_function; input_amplitudes=[1.0],
    input_periods=[1.0], input_duty_cycles=[0.5], input_parameter="I",
    catch_errors=true, show_progress=false)

    # Find all possible combinations of input amplitudes, periods and duty cycles
    vectors = [input_amplitudes, input_periods, input_duty_cycles]
    input_value_combinations = _find_all_combinations(vectors)
    n = size(input_value_combinations, 1)
    
    # Initialize summary matrix
    scan_header = simulation_function()
    scan_results = fill(NaN, n, length(scan_header))

    # Initialize progress bar
    if show_progress
        progressmeter = ProgressMeter.Progress(n; barlen=20)
    end

    # Iterate all parameter value combinations
    lk = ReentrantLock()
    Threads.@threads for i = 1:n

        # Protect the original model from overwriting
        model2 = deepcopy(model)

        # Extract input amplitude, period, and duty cycle
        input_amplitude = input_value_combinations[i, 1]
        input_period = input_value_combinations[i, 2]
        input_duty_cycle = input_value_combinations[i, 3]

        # Generate input
        if model.problem isa JumpProblem
            end_time = model.problem.prob.tspan[2]
        else
            end_time = model.problem.tspan[2]
        end
        events = create_events_cycle(end_time, input_period, input_duty_cycle)
        set_input!(model2, events, input_parameter)
        set_parameter!(model2, input_parameter, input_amplitude)

        # Calculate parameters of the solution
        try
            scan_results[i, :] = simulation_function(model2)
        catch err
            @warn "An error occured for [I, T, D] = [$input_amplitude, $input_period, $input_duty_cycle]"
            if !catch_errors
                throw(err)
            else
                scan_results[i, :] .= NaN
            end
        end

        # Update progress bar
        if show_progress
            lock(lk) do
                ProgressMeter.next!(progressmeter)
            end
        end

    end

    # Build the output dataframe
    matrix = hcat(input_value_combinations, scan_results)
    input_header = ["input_amplitude", "input_period", "input_duty_cycle"]
    names = vcat(input_header, scan_header)
    df = DataFrame(matrix, names)

    return df

end


"""
`plot_arnold`

Plot Arnold tongue or onion.

**Arguments**
- `df`: `DataFrame` with `input_period`, `input_amplitude`, and
    `input_duty_cycle` columns returned by `scan_arnold` function.
- `type`: `"tongue"` for Arnold tongue, `"onion"` for Arnold onion, or
    `"duty_cycle_tongue"` for duty cycle Arnold tongue.

**Keyword Arguments**
- `property_name`: Name of the column to plot.
- `error_name`: Column with error values.
- `error_range`: Minimal and maximal value, specified as an array, allowed for
    `error_name` column.
- `fixed_value`: Fixed value for the third dimension of the "Arnold" space.
- `color_limits`: Color limits specified as `[cmin, cmax]`.
- `show_colorbar`: If `true`, show a colorbar.
- `colorbar_label`: Label for the colorbar.
- `colormap`: Select colormap. For example `"twilight"` for cyclic data. Default
    value is `"viridis"`.
- `figure_title`: String to be displayed as the figure title.
- `ax`: PyPlot axes.
"""
function plot_arnold(df::DataFrame, type="tongue"; property_name=nothing, error_name=nothing,
    error_range=nothing, fixed_value=nothing, color_limits=nothing,
    show_colorbar=true, colorbar_label=nothing, colormap="viridis", 
    figure_title=nothing, ax=gca())

    if type == "tongue"
        x_axis_name = "input_period"
        y_axis_name = "input_amplitude"
        fixed_name = "input_duty_cycle"
    elseif type == "onion"
        x_axis_name = "input_period"
        y_axis_name = "input_duty_cycle"
        fixed_name = "input_amplitude"
    elseif type == "duty_cycle_tongue"
        x_axis_name = "input_duty_cycle"
        y_axis_name = "input_amplitude"
        fixed_name = "input_period"
    end

    # Select data for the given input photoperiod
    if isnothing(fixed_value)
        fixed_value = df[1, fixed_name]
    end
    idx = df[:, fixed_name] .== fixed_value
    if sum(idx) == 0
        throw("No rows for input_duty_cycle=$fixed_value")
    end
    df = df[idx, :]

    # Prepare axes - x-axis: input period, y-axis: input amplitude
    n_x_values = length(unique(df[:, x_axis_name]))
    n_y_values = length(unique(df[:, y_axis_name]))
    if type in ["tongue", "duty_cycle_tongue"]
        x_values_matrix = Matrix(transpose(reshape(df[:, x_axis_name], n_x_values, n_y_values)))
        y_values_matrix = Matrix(transpose(reshape(df[:, y_axis_name], n_x_values, n_y_values)))
    else
        x_values_matrix = Matrix(reshape(df[:, x_axis_name], n_y_values, n_x_values))
        y_values_matrix = Matrix(reshape(df[:, y_axis_name], n_y_values, n_x_values))
    end

    # Estimate which points are entrained
    if !isnothing(error_name)
        error_matrix = reshape(df[:, error_name], n_y_values, n_x_values)
        idx = error_range[1] .<= error_matrix .<= error_range[2]
    else
        idx = fill(true, n_y_values, n_x_values)
    end

    # Plot Arnold tongue
    if isnothing(property_name)
        property_name = "Entrainment"
        Z = fill(1, n_y_values, n_x_values)
        Z[idx] .= 0
        h = ax.pcolor(x_values_matrix, y_values_matrix, Z, rasterized=true, cmap=colormap)
    else
        if type in ["tongue", "duty_cycle_tongue"]
            Z = Matrix(transpose(reshape(df[:, property_name], n_x_values, n_y_values)))
        else
            Z = Matrix(reshape(df[:, property_name], n_y_values, n_x_values))
        end
        Z[.!idx] .= NaN
        h = ax.pcolor(x_values_matrix, y_values_matrix, Z, rasterized=true, cmap=colormap)
    end

    # Set color limits
    if !isnothing(color_limits)
        h.set_clim(color_limits[1], color_limits[2])
    end
    
    # Plot colorbar
    if show_colorbar
        cbar = colorbar(h, ax=ax)
        if isnothing(colorbar_label)
            cbar_label = replace(property_name, "_" => " ")
            cbar_label = uppercase(cbar_label[1]) * cbar_label[2:end]
            cbar.set_label(cbar_label, labelpad=0)
        else
            cbar.set_label(colorbar_label, labelpad=0)
        end
    end

    # Label axes
    if !isnothing(figure_title)
        ax.set_title(figure_title, loc="left", pad=0)
    end

    if type == "tongue"
        if isnothing(figure_title)
            ax.set_title("Arnold tongue", loc="left", pad=0)
        end
        ax.set_xlabel("Input period", labelpad=0)
        ax.set_ylabel("Input amplitude", labelpad=0)
    elseif type == "onion"
        if isnothing(figure_title)
            ax.set_title("Arnold onion", loc="left", pad=0)
        end
        ax.set_xlabel("Input period", labelpad=0)
        ax.set_ylabel("Input duty cycle", labelpad=0)
    elseif type == "duty_cycle_tongue"
        if isnothing(figure_title)
            ax.set_title("Arnold tongue (duty cycle)", loc="left", pad=0)
        end
        ax.set_xlabel("Input duty cycle", labelpad=0)
        ax.set_ylabel("Input amplitude", labelpad=0)
    end

    return h

end


"""
`select_arnold_row`

Find a specific row in the Arnold tongue DataFrame.

**Arguments**
- `df`: `DataFrame` representing Arnold tongue.
- `amplitude`: Input amplitude of the Arnold tongue.
- `period`: Input period of the Arnold tongue.
- `duty_cycle`: Input duty cycle of the Arnold tongue.

**Returns**
- `selected_row`: Row from `df` that is the closest to the requested arguments.
"""
function select_arnold_row(df::DataFrame, amplitude, period, duty_cycle)

    df_copy = deepcopy(df)
    df_copy[!, 1:3] = Float64.(df_copy[!, 1:3])
    df_copy[!, "input_amplitude"] = abs.(df_copy[!, "input_amplitude"] .- amplitude)
    df_copy[!, "input_period"] = abs.(df_copy[!, "input_period"] .- period)
    df_copy[!, "input_duty_cycle"] = abs.(df_copy[!, "input_duty_cycle"] .- duty_cycle)
    indices = sortperm(df_copy, ["input_amplitude", "input_period", "input_duty_cycle"])
    selected_row = DataFrame(df[indices, :][1, :])
    return selected_row

end


"""
`estimate_prc`

Estimate the phase response curve.

**Arguments**
- `model`: Model.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population.
- `n_pulses`: Number of pulses to apply in the cycle.
- `frp`: Estimated free-running period.
- `pacing_offset`: Number of pacing periods used to entrain the model.
- `pacing_length`: Number of pacing periods after `pacing_offset` used to
    estimate the entrained phase of the model.
- `offset`: Number of periods to skip after the pacing before calculating the
    reference cycle.
- `response_length`: Number of periods after the reference cycle used to estimate
    the phase shift.
- `input_parameter`: Name of the parameter on which is the input signal acting.
- `pulse_length`: Length of the input pulse.
- `minprominence`: Minimal prominence threshold to detect a peak.
- `show_progress`: If `true`, show a progress bar in the terminal.
- `show_plots`: If `true`, show plots visualizing the PRC estimation.

**Returns**
- `PRC`: Dataframe representing the estimated PRC. The first column correspond to
    the times of the cycle at which a pulse was applied, normalized to interval
    [0, 1]. The second column represents phase shift normalized to interval
    [-0.5, 0.5]. Negative phase shift means that the peak occured before the free
    running reference.
"""
function estimate_prc(model; trajectories=1, n_pulses=10, frp=1, pacing_offset=10,
    pacing_length=10, offset=3, response_length=3, input_parameter="I",
    pulse_length=nothing, minprominence=0.1, show_progress=false, show_plots=false)

    # Copy model so the original is not modified
    model_prc = deepcopy(model)

    # Set pulse length, if not passed
    if isnothing(pulse_length)
        pulse_length = 0.5*frp
    end

    # Set maximal integration time, 2 is for the reference cycle and reserve
    max_time = (pacing_offset + pacing_length + offset + response_length + 2) * frp
    set_timespan!(model_prc, max_time)

    # Set input events
    events = create_events(:LD, pacing_offset + pacing_length, 0.5frp, 0.5frp)
    set_input!(model_prc, events, input_parameter)

    # Simulate reference
    solution = simulate_population(model_prc, trajectories; save_trajectories=false)
    if show_plots
        t = solution.time
        x = solution.mean[:, 1]
        sol_events = solution.events
        fig, axs = subplots(2, 1)
        axs[1].plot(t, x; color="black")
        plot_events(sol_events, ax=axs[1])
        axs[1].set_title("Full simulation")
    end
    
    # Remove offset
    solution = select_time(solution, min_time=pacing_offset * frp)
    if show_plots
        t = solution.time
        x = solution.mean[:, 1]
        sol_events = solution.events
        axs[2].plot(t, x; color="black")
        plot_events(sol_events, ax=axs[2])
    end

    # Find peaks
    t = solution.time
    x = solution.mean[:, 1]
    pr = findpeaks(x, t, sortstr="descend", sortref="prominence", npeaks=1)
    ref_prom = peakprominences(pr)[1]
    pr = findpeaks(x, t, minprominence=minprominence * ref_prom)
    pks = peakheights(pr)
    locs = peaklocations(pr)
    if show_plots
        axs[2].plot(locs, pks, "o", color="red")
        axs[2].set_title("Removed transient period")
    end

    # Estimate entrainment phase
    sol_events = solution.events
    phase_arr = Float64[]
    for i in 1:pacing_length
        time_start = sol_events[i, 1]
        time_end = time_start + frp
        loc_arr = locs[time_start .<= locs .< time_end]
        if length(loc_arr) > 1
            throw("Too many peaks in the cycle, please debug!")
        elseif length(loc_arr) == 1
            push!(phase_arr, loc_arr[1] - time_start)
        end
    end
    entrainment_phase = mean(phase_arr)

    ## Find the reference cycle ===============================================
    events_end = pacing_length * frp
    idx = locs .> events_end
    locs = locs[idx]
    pks = pks[idx]
    reference_start = locs[offset+1] - entrainment_phase
    reference_end = locs[offset+2] - entrainment_phase
    reference_start_original = reference_start + (pacing_offset * frp)
    reference_end_original = reference_end + (pacing_offset * frp)
    reference_frp = reference_end - reference_start
    trajectory_start = reference_end + 2*pulse_length
    trajectory_start_original = trajectory_start + (pacing_offset * frp)
    idx = t .> trajectory_start
    trajectory_time = t[idx]
    trajectory_reference = x[idx]
    n_trajectory_reference = length(trajectory_reference)
    if show_plots
        axs[2].plot([events_end, events_end], [minimum(x), maximum(x)], "--", color="blue")
        axs[2].plot(locs, pks, "x", color="blue")
        axs[2].plot([reference_start, reference_start], [minimum(x), maximum(x)], color="blue")
        axs[2].plot([reference_end, reference_end], [minimum(x), maximum(x)], color="blue")
        axs[2].plot(trajectory_time, trajectory_reference, color="blue")
        axs[1].plot([trajectory_start_original, trajectory_start_original], [minimum(x), maximum(x)], "--", color="blue")
        axs[1].plot([reference_start_original, reference_start_original], [minimum(x), maximum(x)], color="blue")
        axs[1].plot([reference_end_original, reference_end_original], [minimum(x), maximum(x)], color="blue")
        fig.tight_layout()
    end
    
    ## Simulate pulse responses ===============================================
    # Create a vector with all pulse times along the cycle
    pulse_times = Vector(range(reference_start_original, reference_end_original, length=n_pulses))

    # Prepare figures and axes for the visual control
    if show_plots
        fig, ax_arr = subplots(length(pulse_times), 1)
        fig_corr, ax_arr_corr = subplots(length(pulse_times), 1)
    end

    # Initialize progress bar
    if show_progress
        progressmeter = ProgressMeter.Progress(length(pulse_times); barlen=20)
    end

    # Estimated the phase shifts for the different pulse times along the cycle
    phase_shifts = fill(NaN, length(pulse_times))
    for (i, pulse_time) in enumerate(pulse_times)

        ## Simulation with a pulse ===
        # Add a perturbation pulse to the original pacing events
        events_pulse = vcat(events, [pulse_time pulse_time+pulse_length])
        set_input!(model_prc, events_pulse, input_parameter)
        
        # Simulate the model with the perturbation pulse
        solution = simulate_population(model_prc, trajectories; save_trajectories=false)
        t = solution.time
        x = solution.mean[:, 1]

        # Extract the last part to be compared with the reference
        nx = length(x)
        trajectory = x[nx-n_trajectory_reference+1:nx]

        # Plot the simulation for the visual verification
        if show_plots
            sol_events = solution.events
            ax_arr[i].plot(t, x; color="black")
            ax_arr[i].plot(trajectory_time .+ (pacing_offset * frp), trajectory; color="blue")
            plot_events(events_pulse, ax=ax_arr[i])
        end

        ## Phase shift calculation ===
        # Calculate cross-correlation between the perturbated and reference trajectory
        lags = (-length(trajectory)+1) : (length(trajectory)-1)
        R = crosscor(trajectory, trajectory_reference, lags)

        # Map lags to times and find peaks in the cross-correlation function
        R_time = trajectory_time .- trajectory_time[1]
        R_time = [-R_time[end:-1:2]..., 0, R_time[2:end]...]
        pr = findpeaks(R, R_time)

        # Find peak in the cross-correlation function that is the closest to 0
        idx = argmin(abs.(peaklocations(pr)))
        pk = peakheights(pr)[idx]
        loc = peaklocations(pr)[idx]

        # Normalize the phase shift by the reference FRP
        phase_shifts[i] = loc / reference_frp

        # Plot the cross-correlation function and the chosen peak
        if show_plots
            ax_arr_corr[i].plot(R_time, R, color="black")
            ax_arr_corr[i].plot(loc, pk, "o", color="red")
        end

        # Update progress bar
        if show_progress
            ProgressMeter.next!(progressmeter)
        end

    end

    # Normalize pulse times on the unit cycle
    cycle_times = (pulse_times .- pulse_times[1]) ./ reference_frp

    # Save the estimated PRC as a DataFrame
    PRC = DataFrame(pulse_time=cycle_times, phase_shift=phase_shifts)

    # Plot the estimated PRC
    if show_plots
        fig, ax = subplots()
        ax.plot(PRC[!, "pulse_time"], PRC[!, "phase_shift"], "o")
    end

    return PRC

end


"""
`estimate_T_prc`

Estimate T-cycle phase response curve.

**Arguments**
- `model`: Model.
- `T_cycles`: Array of periods for which the phase response will be estimated.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population.
- `input_parameter`: Name of the parameter on which is the input signal acting.
- `show_progress`: If `true`, show a progress bar in the terminal.

**Returns**
- `df`: Dataframe representing the estimated T-cycle PRC. The first column
    correspond to T_cycles, second to the estimated phases, and third to phase
    coherence for the individual estimations.
"""
function estimate_T_prc(model, T_cycles; trajectories=1, input_parameter="I",
    show_progress=false)

    # Prepare arrays to save the results
    n_T_cycles = length(T_cycles)
    phase_coherence_array = fill(NaN, n_T_cycles)
    collective_phase_array = fill(NaN, n_T_cycles)

    # Initialize progress bar
    if show_progress
        progressmeter = ProgressMeter.Progress(n_T_cycles; barlen=20)
    end

    # Iterate T cycles
    for i_T_cycle = 1:n_T_cycles

        T_cycle = T_cycles[i_T_cycle]
        
        # Create forcing events
        if model.problem isa JumpProblem
            end_time = model.problem.prob.tspan[2]
        else
            end_time = model.problem.tspan[2]
        end
        events = create_events_cycle(end_time, T_cycle)
        set_input!(model, events, input_parameter)

        # Perform simulation
        solution = simulate_population(model, trajectories)
        solution = select_time(solution, min_time=end_time/2)
        t = solution.time
        x = solution.mean[:, 1]
        events = solution.events

        # Estimate the phase of entrainment
        phase_array = estimate_phase_array(t, x, events)
        phase_coherence, collective_phase = estimate_order_parameter(phase_array)
        phase_coherence_array[i_T_cycle] = phase_coherence
        collective_phase_array[i_T_cycle] = collective_phase

        # Update progress bar
        if show_progress
            ProgressMeter.next!(progressmeter)
        end

    end

    # Output the results as a dataframe
    df = DataFrame(T_cycle=T_cycles, phase=collective_phase_array, phase_coherence=phase_coherence_array)
    return df

end


"""
`estimate_jet_lag`

Estimate phase response on a jet-lag-like phase reversal.

**Arguments**
- `model`: Model.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population.
- `input_parameter`: Name of the parameter on which is the input signal acting.
- `offset_days`: Number of days before the start of the experiment to let the oscillator entrain.
- `pre_days`: Number of days to record before the jet lag.
- `post_days`: Number of days to record after the jet lag.
- `phase_estimator`: Select phase estimator for `estimate_phase_array`.
- `show_plots`: If `true`, show plots visualizing the PRC estimation.

**Returns**
- `phases`: Vector of recorded phases centered by the default phase of entrainment.
- `days`: Time vector for the phases. Negative days indicate days before the jet lag and positive days indicate days after the jet lag.
"""
function estimate_jet_lag(model::Model; trajectories=1, input_parameter="I", offset_days=20, pre_days=5, post_days=20, phase_estimator="peak_height", show_plots=false)
    
    # Protect the original model from overwritting
    model = deepcopy(model)

    # Create jet lag events
    events = create_events([
        (:LD, offset_days + pre_days, 0.5, 0.5),  # day-night cycle before jet lag
        (:LD, 1, 1, 0.5),  # jet lag
        (:LD, post_days + 1, 0.5, 0.5),  # day-night cycle after jet lag
    ])

    # Set the timespan and input signal of the model
    time_end = events[end, 2] + 0.5
    set_timespan!(model, time_end)
    set_input!(model, events, input_parameter)

    # Simulate the model
    solution = simulate_population(model, trajectories; show_progress=true)

    # Extract the phase for the days before the jet lag
    events_pre = solution.events[offset_days:(offset_days + pre_days + 1), :]
    idx = events_pre[1, 1] .<= solution.time .<= events_pre[end, 2] + 0.5
    t_pre = solution.time[idx]
    x_pre = solution.mean[idx, 1]
    phase_array_pre = estimate_phase_array(t_pre, x_pre, events_pre, method=phase_estimator)
    entrainment_phase = cmean(phase_array_pre)
    phase_array_pre .-= entrainment_phase

    # Extract the phase for the days after the jet lag
    events_post = solution.events[(offset_days + pre_days + 1):end, :]
    idx = events_post[1, 1] .<= solution.time .<= events_post[end, 2] + 0.5
    t_post = solution.time[idx]
    x_post = solution.mean[idx, 1]
    phase_array_post = estimate_phase_array(t_post, x_post, events_post, method=phase_estimator)
    phase_array_post .-= entrainment_phase

    # Create the output arrays of phases and days (time vector)
    days = Vector(-pre_days:post_days)
    phases = vcat(phase_array_pre, [NaN], phase_array_post)

    # Plot the results (simulated time series and estimated phase)
    if show_plots
        fig, ax = subplots()
        plot_solution(solution, ax=ax)
        ax.plot(t_pre, x_pre)
        ax.plot(t_post, x_post)

        fig, ax = subplots()
        ax.plot([days[1] - 0.5, days[end] + 0.5], [0, 0], "k--")
        ax.plot([days[pre_days + 1], days[pre_days + 1]], [minimum(phases[.!isnan.(phases)]) - 0.1, maximum(.!isnan.(phases)) + 0.1], "k--")
        ax.plot(days, phases, "ko")
    end

    return phases, days

end
