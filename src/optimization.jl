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
- `input_parameter_name`: Default `I`.
- `show_plots`: If `true`, plot the model simulation vs the data. Default value
    is `false`.

**Returns**
- `cost_function`: A cost function that takes an array of parameter values as an
    input argument and returns the squared error of the simulation from the data.
"""
function create_entrainable_oscillator_objective(model, reference_period;
    trajectories=1, parameter_names=nothing, input_parameter_name="I",
    show_plots=false)

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
        solution = select_time(solution, min_time=min_time)
        t_dd = solution.time
        x_dd = solution.mean[:, 1]
        y_dd = solution.mean[:, 2]
        amp_dd = (maximum(x_dd) - minimum(x_dd))
        if amp_dd < 0.1
            return Inf
        end

        period = estimate_period_winding_number(x_dd, y_dd, min_time)
        E1 = abs(period - reference_period) / reference_period

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

        E2 = abs(collective_phase1 - collective_phase2) / period

        if show_plots
            fig, ax_array = subplots(3)
            ax_array[1].plot(t_dd, x_dd, color="black")
            ax_array[2].plot(t_ld, x_ld, color="black")
            plot_events(events_ld, ax=ax_array[2])
            ax_array[3].plot(t_ld2, x_ld2, color="black")
            plot_events(events_ld2, ax=ax_array[3])
            fig.tight_layout()
        end

        return E1 + E2

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
    function cost_function(parameter_values; show_plots=show_plots)
        
        # Set parameters of the model
        model2 = deepcopy(model)
        set_parameter!(model2, noise_parameter_name, parameter_values[1])

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

        # Show plot of the simulated trajectory vs the actual data
        if show_plots

            _, ax = subplots()
            ax.plot(t, x, color="black", label="Simulation")
            ax.plot(t, damped_sine(t, p), color="blue", label="Curve fit")
            plot_events(solution.events, ax=ax)

        end

        # Return the squared error
        return (estimated_damping_rate - damping_rate) ^2

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
- `population_size`: Population size.
- `initial_population`: Initial population passed to the optimizer. The rows
    are individuals and columns are parameter values. This option overwrites 
    `population_size`, if both are passed.
- `trace_mode`: `:compact` (defualt), `:silent` or `:verbose`.

**Returns**
- `best_candidate`: A vector representing the best candidate solution.
- `final_population`: The final population of candidates.
"""
function optimize(cost_function; search_range=nothing, max_steps=nothing,
    population_size=nothing, initial_population=nothing, trace_mode=nothing)

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

    if !isnothing(population_size)
        optimizer_kwargs[:PopulationSize] = population_size
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


"""
`damped_sine`

Evaluated a damped sine function.

**Arguments**
- `t`: Time point or vector.
- `p`: Vector of parameter values of form `p = [A, d, T, θ]` where `A` is
    amplitude, `d` damping ratio, `T` period and `θ` phase.
        
**Returns**
- `x`: Values of the damped sine for the time points `t`.
"""
function damped_sine(t, p)
    A, d, T, θ = p
    x = A .* exp.(-d.*t) .* sin.( (2*pi*t) ./ T .+ θ)
    return x
end


"""
`polynomial`

Evaluated a polynomial function.

**Arguments**
- `t`: Time point or vector.
- `p`: Vector of parameter values of form `p = [a₀, a₁, a₂, ...]` where the
    parameters correspond to polynomial `a₀ + a₁t + a₂t² + ...`.
        
**Returns**
- `x`: Values of the damped sine for the time points `t`.
"""
function polynomial(t, p)

    nt = length(t)
    np = length(p)

    x = fill(0.0, nt)

    for i = 1:nt
        for j = 1:np
            x[i] += p[j] * t[i]^(j-1)
        end
    end

    return x

end


"""
`fit_curve`

Fit curve to data.

**Arguments**
- `fun`: Function that takes two arguments, time vector and parameter values and
    returns the evaluated function values for each time point.
- `t`: Time vector.
- `x`: Fitting data vector.
- `p0`: Initial guess for the parameter values.
        
**Returns**
- `p`: Fitted parameter values.
"""
function fit_curve(fun::Function, t, x, p0)

    lsqfit = LsqFit.curve_fit(fun, t, x, p0)
    return lsqfit.param

end
