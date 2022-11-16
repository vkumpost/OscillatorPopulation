"""
`create_data_objective`

Create a cost function for fitting a model to data.

**Arguments**
- `model`: `Model`.
- `t`: Time vector.
- `x`: Data vector.

**Optional Arguments**
- `events`: Events matrix.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population. Default value is 1.
- `parameter_names`: Names of the model parameters to be optimized. By default,
    all model parameters are optimized.
- `initial_conditions`: Matrix indicating initial conditions for each trajectory
    of the population. Columns are variables and rows are time points. By
    default, the initial conditions from `model` are applied.
- `pre_events_description`: Events described by a tuple (see `create_events`) to
    be applied before the fitting to reduce the effect of the initial conditions.
    Be default, no additional events are used.
- `input_parameter_name`: Default `I`.
- `show_plots`: If `true`, plot the model simulation vs the data. Default value
    is `false`.
- `save_plot_filename`: A string indicating the path to the image file, where the
    figure should be saved. By default is the figure not saved.

**Returns**
- `cost_function`: A cost function that takes an array of parameter values as an
    input argument and returns the squared error of the simulation from the data.
"""
function create_data_objective(model, t, x, events=nothing; trajectories=1,
    parameter_names=nothing, initial_conditions=nothing,
    pre_events_description=nothing, show_plots=false, input_parameter_name="I",
    save_plot_filename=nothing)

    # Create copies of the input parameters
    model = deepcopy(model)
    t = deepcopy(t)
    x = deepcopy(x)

    if isnothing(events)
        events = create_events()
    else
        events = deepcopy(events)
    end

    # If no parameter names are passed, use all parameters of the model
    if isnothing(parameter_names)
        parameter_names = model.parameter_names
    end

    # Extend the simulation by additional events preceding the data-fitting time
    if !isnothing(pre_events_description)
        if pre_events_description[1] in (:DD, :LL)
            pre_end_time = pre_events_description[2]
        else  # :LD
            pre_end_time = pre_events_description[2] * (pre_events_description[3] + pre_events_description[4])
        end
        pre_events = create_events(pre_events_description...)

        fit_events = vcat(pre_events, events .+ pre_end_time)
        fit_t = t .+ pre_end_time
        fit_x = x
    else
        pre_end_time = 0
        fit_events = events
        fit_t = t
        fit_x = x
    end

    # Configurate the model for data fitting
    set_timespan!(model, maximum(fit_t))
    set_input!(model, fit_events, input_parameter_name)
    set_solver!(model, saveat=fit_t)

    # Cost function
    function cost_function(parameter_values; show_plots=show_plots, save_plot_filename=save_plot_filename)
        
        # Set parameters of the model
        model2 = deepcopy(model)
        set_parameter!(model2, parameter_names, parameter_values)

        # Simulate the model
        solution = nothing
        try
            solution = simulate_population(model2, trajectories; initial_conditions=initial_conditions)
            solution = select_time(solution, min_time=pre_end_time)
        catch err
            @warn "Error during simulation:\n$(err)"
            return Inf
        end
        if !solution.success
            @warn "Unsuccessful simulation!"
            return Inf
        end

        # Show plot of the simulated trajectory vs the actual data
        if show_plots || !isnothing(save_plot_filename)
            fig, ax = subplots()
            ax.plot(t, x, color="black", label="Data")
            ax.plot(solution.time, solution.mean, color="blue", label="Simulation")
            plot_events(solution.events, ax=ax)

            if !isnothing(save_plot_filename)
                fig.savefig(save_plot_filename)
            end

            if !show_plots
                close(fig)
            end

        end

        # Return the mean squared error
        return mean( (x .- solution.mean) .^2 )

    end

    return cost_function

end


"""
`create_entrainable_oscillator_objective`

Create a cost function to find parameters for an entrainable limit cycle
oscillator.

**Arguments**
- `model`: `Model`.
- `reference_period`: Target period for oscillations.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population. Default value is 1.
- `parameter_names`: Names of the model parameters to be optimized. By default,
    all model parameters are optimized.
- `entrainment_extrema`: Array with two elements indicating the minimum and
    maximum of the entrained oscillations. By default, entrainment extrema are
    not considered as an  objective.
- `input_parameter_name`: Name of the input parameter. Default value is `I`.
- `show_plots`: If `true`, plot the model simulation vs the data. Default value
    is `false`.

**Returns**
- `cost_function`: A cost function that takes an array of parameter values as an
    input argument and returns the squared error of the simulation from the data.
"""
function create_entrainable_oscillator_objective(model, reference_period;
    trajectories=1, parameter_names=nothing, entrainment_extrema=nothing,
    input_parameter_name="I", show_plots=false)

    model = deepcopy(model)

    if isnothing(parameter_names)
        parameter_names = model.parameter_names
    end

    max_time = reference_period * 20
    min_time = reference_period * 10
    set_timespan!(model, max_time)

    function cost_function(parameter_values; show_plots=show_plots)

        model2 = deepcopy(model)

        # Free running period
        set_parameter!(model2, parameter_names, parameter_values)
        set_parameter!(model2, input_parameter_name, 0.0)

        solution = nothing
        try
            solution = simulate_population(model2, trajectories)
        catch e
            return Inf
        end
        if !solution.success
            return Inf
        end
        solution = select_time(solution, min_time=min_time)
        t_dd = solution.time
        x_dd = solution.mean[:, 1]
        amp_dd = (maximum(x_dd) - minimum(x_dd))
        if amp_dd < 0.1
            return Inf
        end

        period, _ = estimate_period(t_dd, x_dd)
        E1 = abs(period - reference_period)
        if E1 < reference_period * 0.1
            E1 = 0.0
        else
            E1 /= reference_period
        end

        # Entrainment
        set_parameter!(model2, parameter_names, parameter_values)
        events_original = create_events_cycle(max_time, reference_period)
        
        set_input!(model2, events_original, input_parameter_name)
        solution = nothing
        try
            solution = simulate_population(model2, trajectories)
        catch e
            return Inf
        end
        if !solution.success
            return Inf
        end
        solution = select_time(solution, min_time=min_time)
        t_ld = solution.time
        x_ld = solution.mean[:, 1]
        amd_ld = maximum(x_ld) - minimum(x_ld)
        events_ld = solution.events
        phase_array = estimate_phase_array_peaks(t_ld, x_ld, events_ld)
        if length(phase_array) < 5
            return Inf
        else
            phase_coherence1, collective_phase1 = estimate_order_parameter(phase_array)
        end

        set_input!(model2, events_original .+ period/2, input_parameter_name)
        solution = nothing
        try
            solution = simulate_population(model2, trajectories)
        catch e
            return Inf
        end
        if !solution.success
            return Inf
        end
        solution = select_time(solution, min_time=min_time)
        t_ld2 = solution.time
        x_ld2 = solution.mean[:, 1]
        events_ld2 = solution.events
        phase_array = estimate_phase_array_peaks(t_ld2, x_ld2, events_ld2)
        if length(phase_array) < 5
            return Inf
        else
            phase_coherence2, collective_phase2 = estimate_order_parameter(phase_array)
        end

        E2 = abs(collective_phase1 - collective_phase2) / reference_period

        E_total = E1 + E2

        if !isnothing(entrainment_extrema)
            E3 = abs(entrainment_extrema[1] - minimum(solution.mean[:, 1]))
            if E3 < entrainment_extrema[1] * 0.1
                E3 = 0.0
            else
                E3 /=  entrainment_extrema[1]
            end
            E4 = abs(entrainment_extrema[2] - maximum(solution.mean[:, 1]))
            if E4 < entrainment_extrema[2] * 0.1
                E4 = 0.0
            else
                E4 /= entrainment_extrema[2]
            end
            E_total = E_total + E3 + E4
        end

        if show_plots
            fig, ax_array = subplots(3)
            ax_array[1].plot(t_dd, x_dd, color="black")
            ax_array[2].plot(t_ld, x_ld, color="black")
            plot_events(events_ld, ax=ax_array[2])
            ax_array[3].plot(t_ld2, x_ld2, color="black")
            plot_events(events_ld2, ax=ax_array[3])
            fig.tight_layout()
        end

        return E_total

    end

    return cost_function

end


"""
`create_desynchronization_objective`

Create a function that estimates the damping rate.

**Arguments**
- `model`: `Model`.
- `forcing_period`: Period of the forcing signal to entrain the population.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population. Default value is 1.
- `noise_parameter_name`: Name of the noise parameter. Default value is `σ`.
- `input_parameter_name`: Name of the input parameter. Default values is `I`.
- `show_plots`: If `true`, plot the model simulation vs the data. Default value
    is `false`.

**Returns**
- `cost_function`: A cost function that takes an array of parameter values as an
    input argument and returns the estimated damping rate.
"""
function create_desynchronization_objective(model, forcing_period;
    trajectories=1, noise_parameter_name="σ", input_parameter_name="I",
    show_plots=false
    )

    # Create a copy of the input model
    model = deepcopy(model)

    events = create_events(:LD, 10, forcing_period/2, forcing_period/2)
    min_time = 10*forcing_period
    max_time = min_time + 5*forcing_period

    # Configurate the model for data fitting
    set_timespan!(model, max_time)
    set_input!(model, events, input_parameter_name)
    set_solver!(model, saveat=range(min_time, max_time, 100))

    # Cost function
    function cost_function(noise_intensity; show_plots=show_plots)
        
        # Set parameters of the model
        model2 = deepcopy(model)
        set_parameter!(model2, noise_parameter_name, noise_intensity[1])

        # Simulate the model
        solution = nothing
        try
            solution = simulate_population(model2, trajectories)
            solution = select_time(solution, min_time=min_time)
        catch err
            @warn "Error during simulation:\n$(err)"
            return Inf
        end
        if !solution.success
            @warn "Unsuccessful simulation!"
            return Inf
        end

        t = solution.time
        x = zscore(solution.mean[:, 1])
        T0 = forcing_period
        A0 = maximum(x)
        p0 = [A0, 0.0, T0, 0.0]
        p = fit_curve(damped_sine, t, x, p0)
        estimated_damping_rate = p[2]

        xx = damped_sine(t, p)
        R = rsquared(x, xx)

        # Show plot of the simulated trajectory vs the actual data
        if show_plots

            _, ax = subplots()
            ax.plot(t, x, color="black", label="Simulation")
            ax.plot(t, xx, color="blue", label="Curve fit")
            plot_events(solution.events, ax=ax)

        end

        # Return the squared error
        if R > 0.8
            return estimated_damping_rate
        else
            return Inf
        end

    end

    return cost_function

end


"""
`create_desynchronization_objective`

Create a cost function for fitting noise intensity to the population-level 
damping rate.

**Arguments**
- `model`: `Model`.
- `damping_rate`: Damping rate of the damped sine fitted to the population-level
    data.
- `forcing_period`: Period of the forcing signal to entrain the population.

**Keyword Arguments**
- `trajectories`: Number of trajectories in the population. Default value is 1.
- `noise_parameter_name`: Name of the noise parameter. Default value is `σ`.
- `input_parameter_name`: Name of the input parameter. Default values is `I`.
- `show_plots`: If `true`, plot the model simulation vs the data. Default value
    is `false`.

**Returns**
- `cost_function`: A cost function that takes an array of parameter values as an
    input argument and returns the squared error of the simulated damping rate
    from the reference damping rate.
"""
function create_desynchronization_objective(model, damping_rate, forcing_period;
    trajectories=1, noise_parameter_name="σ", input_parameter_name="I",
    show_plots=false
    )

    cost_function_tmp = create_desynchronization_objective(model, forcing_period;
        trajectories=trajectories, noise_parameter_name=noise_parameter_name,
        input_parameter_name=input_parameter_name, show_plots=show_plots
    )
    
    # Cost function
    function cost_function(parameter_values; show_plots=show_plots)
        
        estimated_damping_rate = cost_function_tmp(parameter_values)

            return (estimated_damping_rate - damping_rate) ^2

    end

    return cost_function

end
