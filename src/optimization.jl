"""
`estimate_initial_conditions`

Estimate initial conditions for a population of entrained or free-running
oscillators.

**Arguments**
- `model`: Model.
- `trajectories`: Number of trajectories to simulate.
- `events_description`: A tuple passed to `create_events`.

**Keyword Arguments**
- `input_parameter_name`: Name of the forcing parameter. Default values is `I`.
- `kwargs...`: Any additional arguments are passed to `simulate_population`.

**Returns**
- `initial_conditions`: Matrix of initial conditions. Columns are variables and
    rows are time points. If forcing is applied, the initial conditions are
    aligned to the beginning of the last event.
"""
function estimate_initial_conditions(model, trajectories, events_description;
    input_parameter_name="I", kwargs...)

    model = deepcopy(model)

    if events_description[1] in (:DD, :LL)
        max_time = events_description[2]
    else  # :LD
        max_time = events_description[2] * (events_description[3] + events_description[4])
    end
    events = create_events(events_description...)

    set_timespan!(model, max_time)
    set_input!(model, events, input_parameter_name)
    set_solver!(model, saveat=max_time)

    # Set output to extract the individual state variables
    variable_indices = 1:length(model.variable_names)
    set_output!(model, variable_indices)
        
    # Run simulation and extract the initial conditions
    solution = simulate_population(model, trajectories; kwargs...)
    initial_conditions = solution.trajectories[end, :, :]
    initial_conditions = Matrix(transpose(initial_conditions))

    return initial_conditions

end


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

        # Return the squared error
        return sum( (x .- solution.mean) .^2 )

    end

    return cost_function

end


"""
`optimize`

Optimize a cost function.

**Arguments**
- `cost_function`: Function to be minimized. It takes a vector of parameter
    values as input and outputs the cost function value.

**Keyword Arguments**
- `search_range`: Vector of tuples that specify the search range for the
    individual parameters. This keyword argument must be passed!
- `max_steps`: The maximal number of steps that the optimizer can take.
- `initial_population`: Initial population passed to the optimizer. The rows
    are individuals and columns are parameter values.
- `trace_mode`: `:compact` (defualt), `:silent` or `:verbose`.

**Returns**
- `best_candidate`: A vector representing the best candidate solution.
- `final_population`: The final population of candidates.
"""
function optimize(cost_function; search_range=nothing, max_steps=nothing,
    initial_population=nothing, trace_mode=nothing)

    optimizer_kwargs = Dict{Symbol, Any}(
        :TraceMode => :compact,
    )

    if !isnothing(search_range)
        optimizer_kwargs[:SearchRange] = search_range
        optimizer_kwargs[:NumDimensions] = length(search_range)
    end

    if !isnothing(max_steps)
        optimizer_kwargs[:MaxSteps] = max_steps
    end

    if !isnothing(initial_population)
        optimizer_kwargs[:PopulationSize] = size(initial_population, 1)
        optimizer_kwargs[:Population] = Matrix(transpose(initial_population))
    end

    if !isnothing(trace_mode)
        optimizer_kwargs[:TraceMode] = trace_mode
    end

    res = BlackBoxOptim.bboptimize(cost_function; optimizer_kwargs...)

    best_candidate = BlackBoxOptim.best_candidate(res)
    final_population = Matrix(transpose(BlackBoxOptim.population(res).individuals))

    return best_candidate, final_population

end
