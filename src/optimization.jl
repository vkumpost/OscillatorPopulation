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
