"""
`estimate_initial_conditions`

Estimate initial conditions for a population of entrained or free-running
oscillators.

**Arguments**
- `model`: Model.

**Optional Arguments**
- `trajectories`: Number of trajectories to simulate. Default value is `1000`.

**Keyword Arguments**
- `max_time`: Time for which is the model simulated before the initial
    conditions are calculated. Should be long enough to minimize the effect
    of the original initial conditions. If not passed, the `max_time` is taken
    from the currect timespan of the model. 
- `forcing_period`: Period of the forcing signal. If not passsed, no forcing is
    applied.
- `input_parameter_name`: Name of the forcing parameter. Default values is `I`.
- `kwargs...`: Any additional arguments are passed to `simulate_population`.

**Returns**
- `initial_conditions`: Matrix of initial conditions. Columns are variables and
    rows are time points. If forcing is applied, the initial conditions are
    aligned to the beginning of the last event.
"""
function estimate_initial_conditions(model, trajectories=1000;
    max_time=nothing, forcing_period=nothing, input_parameter_name="I",
    kwargs...)

    # Do not modify the original model
    model = deepcopy(model)

    # Set the maximal simulation time
    if !isnothing(max_time)
        set_timespan!(model, max_time)
    end

    # Set the forcing period
    if !isnothing(forcing_period)
        t_end = model.problem.tspan[end]
        events = create_events_cycle(t_end, forcing_period)

        # Align the maximal simulation time to the beginning of the last event
        set_timespan!(model, events[end, 1])
        set_input!(model, events, input_parameter_name)
    end

    # Set output to extract the individual state variables
    variable_indices = 1:length(model.variable_names)
    set_output!(model, variable_indices)

    # Run simulation and extract the initial conditions
    solution = simulate_population(model, trajectories; kwargs...)
    initial_conditions = solution.trajectories[end, :, :]
    initial_conditions = Matrix(transpose(initial_conditions))

    return initial_conditions

end
