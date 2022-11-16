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
- `progress_filename`: Name of a csv file to save the optimization progress
    (fitness as a function of the number of steps).

**Returns**
- `best_candidate`: A vector representing the best candidate solution.
- `final_population`: The final population of candidates.
"""
function optimize(cost_function; search_range=nothing, max_steps=nothing,
    population_size=nothing, initial_population=nothing, trace_mode=nothing,
    progress_filename=nothing)

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

    if !isnothing(progress_filename)

        function create_callback_function()

            best_fitness = Inf
            counter = 0
            df = DataFrame(:Step=>Float64[], :Fitness=>Float64[])

            function callback_function(oc)

                counter += 1
                current_fitness = BlackBoxOptim.best_fitness(oc)

                if current_fitness < best_fitness
                    
                    best_fitness = current_fitness
                    push!(df, [counter, best_fitness])
                    save_data(df, progress_filename; force=true)
                    
                end

            end
    
            return callback_function
    
        end

        optimizer_kwargs[:CallbackFunction] = create_callback_function()
        optimizer_kwargs[:CallbackInterval] = 0.0

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


"""
`binary_search`

Use binary search to find `x` for which `fun` evaluates to `target_value`.

**Arguments**
- `fun`: A monotonically increasing function which takes one number as a parameter.
- `target_value`: Targer value to which should `fun` evaluate.
- `search_range`: A two-element array specifying the search range.

**Keyword Arguments**
- `tolerance`: Maximal tolerance from the target value.
- `max_steps`: Maximal number of steps.
- `verbose`: If `true`, print the estimated value and error at each step.

**Returns**
- `x`: A number for which `fun(x) ≈ target_value`.
"""
function binary_search(fun, target_value, search_range; tolerance=√eps(),
    max_steps=1_000_000, verbose=false)
    
    # Get search boundaries
    x_left = search_range[1]
    x_right = search_range[2]

    # Check that the target value lies in the search interval
    if !(fun(x_left) < target_value < fun(x_right))
        msg = "Target value is not in the search interval!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    # Initialized loop variables
    estimated_value = fun(x_left)
    x_new = x_left
    step_counter = 0

    while (abs(estimated_value - target_value) > tolerance) && (step_counter < max_steps)
        
        # Guess new value in the middle of the search interval
        x_new = (x_left + x_right) / 2
        estimated_value = fun(x_new)

        # Adjust the search interval
        if estimated_value < target_value
            x_left = x_new
        else
            x_right = x_new
        end

        # Increase step counter
        step_counter += 1

        # Print step number, estimated value and error
        if verbose
            println("Step $(step_counter), Estimate $(round(x_new, digits=10)), Error $(round(abs(estimated_value - target_value), digits=10))")
        end

    end

    return x_new

end
