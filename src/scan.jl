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
`scan`

Scan parameters of a model.

**Arguments**
- `model`: Model.
- `parameters`: A vector of pairs. Pairs indicate parameter names and values.
- `simulation_function`: A function that takes in `PopulationSolution`, and returns
    some descriptive parameters. If the function is called without any arguments,
    it returns the names of the parameters as a vector of strings. For example:
    ```
    julia> simulation_function(solution)
    [1.1, 2.2, 3.3]
    julia> simulation_function()
    ["Period", "Amplitude", "Phase"]
    ```

**Keyword Arguments**
- `show_progress`: If `true`, show a progress bar in the terminal.

**Returns**
- `df`: Dataframe with results. First columns correspond to the input parameters
    followed by columns returned by `simulation_function` for the corresponding
    parameters.
"""
function scan(model::Model, parameters::Vector, simulation_function::Function;
    show_progress::Bool=false)

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

        # Evaluate the sumary function
        simulation_summary = simulation_function(model2)

        # Save the model summary to the overall matrix
        scan_summary[i, :] = simulation_summary

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
- `simulation_function`: A function that takes in `PopulationSolution`, and returns
    some descriptive parameters. If the function is called without any arguments,
    it returns the names of the parameters as a vector of strings. For example:
    ```
    julia> simulation_function(solution)
    [1.1, 2.2, 3.3]
    julia> simulation_function()
    ["Period", "Amplitude", "Phase"]
    ```

**Keyword Arguments**
- `input_amplitudes`: Input amplitudes to scan.
- `input_periods`: Input periods to scan.
- `input_photoperiods`: Input photoperiods to scan.
- `input_parameter_name`: Name of the parameter that controls the input amplitude.
- `show_progress`: If true, show progress in the terminal.

**Returns**
- `df`: Dataframe with results. First columns correspond to the input amplitude,
    period, and photoperiod. The following columns correspond to the values
    returned by `simulation_function` for the corresponding
    parameters.
"""
function scan_arnold(model, simulation_function; input_amplitudes=[1.0],
    input_periods=[1.0], input_photoperiods=[0.5], input_parameter="I",
    show_progress=false)

    # Find all possible combinations of input amplitudes, periods and photoperiods
    vectors = [input_amplitudes, input_periods, input_photoperiods]
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
    for i = 1:n  # Threads.@threads 

        # Protect the original model from overwriting
        model2 = deepcopy(model)

        # Extract input amplitude, period, and photoperiod
        input_amplitude = input_value_combinations[i, 1]
        input_period = input_value_combinations[i, 2]
        input_photoperiod = input_value_combinations[i, 3]

        # Generate input
        end_time = model.problem.tspan[2]
        events = create_events_cycle(end_time, input_period, input_photoperiod)
        set_input!(model2, events, input_parameter)
        set_parameter!(model2, input_parameter, input_amplitude)

        # Calculate parameters of the solution
        try
            scan_results[i, :] = simulation_function(model2)
        catch err
            @warn "An error occured for [$input_amplitude, $input_period, $input_photoperiod]"
            throw(err)
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
    input_header = ["input_amplitude", "input_period", "input_photoperiod"]
    names = vcat(input_header, scan_header)
    df = DataFrame(matrix, names)

    return df

end


"""
`estimate_prc`

Estimate arnold tongue and/or onion.

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
    pulse_length=nothing, show_plots=false)

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
    pr = findpeaks(x, t, minprominence=0.1ref_prom)
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

    # Find the reference cycle
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
    
    # Simulate pulse responses
    pulse_times = Vector(range(reference_start_original, reference_end_original, length=n_pulses))
    if show_plots
        fig, ax_arr = subplots(length(pulse_times), 1)
        fig_corr, ax_arr_corr = subplots(length(pulse_times), 1)
    end
    phase_shifts = fill(NaN, length(pulse_times))
    for (i, pulse_time) in enumerate(pulse_times)

        # Simulation with a pulse
        events_pulse = vcat(events, [pulse_time pulse_time+pulse_length])
        set_input!(model_prc, events_pulse, input_parameter)
        
        solution = simulate_population(model_prc, trajectories; save_trajectories=false)
        t = solution.time
        x = solution.mean[:, 1]
        idx = t .> trajectory_start_original
        trajectory = x[idx]
        if show_plots
            sol_events = solution.events
            ax_arr[i].plot(t, x; color="black")
            ax_arr[i].plot(trajectory_time .+ (pacing_offset * frp), trajectory; color="blue")
            plot_events(events_pulse, ax=ax_arr[i])
        end

        # Phase shift calculation´
        lags = (-length(trajectory)+1) : (length(trajectory)-1)
        R = crosscor(trajectory_reference, trajectory, lags)
        R_time = trajectory_time .- trajectory_time[1]
        R_time = [-R_time[end:-1:2]..., 0, R_time[2:end]...]
        pr = findpeaks(R, R_time, sortstr="descend", sortref="prominence")
        pk = peakheights(pr)[1]
        loc = peaklocations(pr)[1]
        phase_shifts[i] = loc / reference_frp
        if show_plots
            ax_arr_corr[i].plot(R_time, R, color="black")
            ax_arr_corr[i].plot(loc, pk, "o", color="red")
        end

    end

    cycle_times = (pulse_times .- pulse_times[1]) ./ reference_frp
    PRC = DataFrame(pulse_time=cycle_times, phase_shift=phase_shifts)
    if show_plots
        fig, ax = subplots()
        ax.plot(PRC[!, "pulse_time"], PRC[!, "phase_shift"])
    end

    return PRC

end