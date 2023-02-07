"""
`Model`

A struct representing an oscillatory model.

**Fields**
- `variable_names`: Variable names (e.g. `["X", "Y", "Z"]`).
- `parameter_names`: Parameter names (e.g. `["v", "d1", "d2"]`).
- `problem`: `ODEProblem` or `SDEProblem` or `JumpProblem`.
- `solver_algorithm`: Solver to be used to solve `problem` (e.g. `SOSRI()`).
- `solver_parameters`: Parameters passed to the solver (e.g. `(saveat=0.1,)`)
- `input`: Input to the model represented as a tuple with two elements. The
    first element is a matrix representing the input events and the second
    element is a string that specifies, which parameter is manipulated by the
    events.
- `output`: Function that transforms the solution into a matrix: rows are time
    steps and columns are variables.
"""
mutable struct Model
    variable_names::Vector{String}
    parameter_names::Vector{String}
    problem::Union{ODEProblem, SDEProblem, JumpProblem}
    solver_algorithm::Any
    solver_parameters::NamedTuple
    input::Tuple{Matrix{Float64}, String}
    output::Function
end


"""
`_get_problem_property`

Extract a property of a problem.

**Arguments**
- `problem`: `ODEProblem`, `SDEProblem`, or 'JumpProblem`.
- `property_name`: Property name (e.g. `u0`, `p`, `tspan`).

**Returns**
- `property`: Property of a problem.
"""
function _get_problem_property(problem::Union{ODEProblem, SDEProblem, JumpProblem},
    property_name)

    # Property name must be a symbol
    property_name = Symbol(property_name)

    # For JumpProblem access problem.prob, otherwise problem directly
    if problem isa JumpProblem
        property = getproperty(problem.prob, property_name)
    else
        property = getproperty(problem, property_name)
    end

    # Return the found property
    return property

end


"""
`_get_problem_property`

Extract a property of a problem that is stored part of a model.

**Arguments**
- `model`: `Model`.
- `property_name`: Property name (e.g. `u0`, `p`, `tspan`).

**Returns**
- `property`: Property of a problem.
"""
function _get_problem_property(model::Model, property_name)

    return _get_problem_property(model.problem, property_name)

end


"""
`_set_problem_property!`

Set a propery of a model's problem.

**Arguments**
- `model`: `Model`.

**Keyword Arguments**
- `kwargs...`: Name-value pairs of the properties to be set (e.g.
    `u0=[1.0, 0.0]`, `tspan=(0.0, 10.0)`).

**Returns**
- `model`: In-place edited model.
"""
function _set_problem_property!(model::Model; kwargs...)

    # Remake problem with the new properties
    new_problem = remake(model.problem; kwargs...)

    # Replace the problem in the model
    model.problem = new_problem
    
    return model

end


"""
`set_initial_conditions!`

Set initial conditions for a model.

**Arguments**
- `model`: `Model`.
- `u0`: Initial conditions.

**Returns**
- `model`: In-place edited model.
"""
function set_initial_conditions!(model::Model, u0)

    model = _set_problem_property!(model; u0=u0)
    return model

end


"""
`set_timespan!`

Set timespan for a model.

**Arguments**
- `model`: `Model`.
- `tspan`: Timespan.

**Returns**
- `model`: In-place edited model.
"""
function set_timespan!(model::Model, tspan::Tuple)

    model = _set_problem_property!(model; tspan=tspan)
    return model

end


"""
`set_timespan!`

Set timespan from 0 to some a certain time.

**Arguments**
- `model`: `Model`.
- `max_time`: Maximal time of integration.

**Returns**
- `model`: In-place edited model.
"""
function set_timespan!(model::Model, max_time::Number)

    model = set_timespan!(model, (0.0, max_time))
    return model

end


"""
`set_output!`

Set output for a model.

**Arguments**
- `model`: `Model`.
- `output`: Output function that transforms `DifferentialEqualtions` solution
    into a matrix: rows are time steps and columns are variables.

**Returns**
- `model`: In-place edited model.
"""
function set_output!(model::Model, output::Function)

    model.output = output
    return model

end


"""
`set_output!`

Select output variable for a model.

**Arguments**
- `model`: `Model`.
- `index`: Index of the variable that should serve as output.

**Returns**
- `model`: In-place edited model.
"""
function set_output!(model::Model, index)

    if length(index) == 1
        model.output = sol -> Matrix(Matrix(sol[index, :]')')
    else
        model.output = sol -> Matrix(sol[index, :]')
    end
    return model

end


"""
`set_solver!`

Set solver algorithm and its properties.

**Arguments**
- `model`: Model.

**Optional Arguments**
- `algorithm`: Solver algorithm (e.g. `DP5()` or `SOSRI()`).

**Keyword Arguments**
- `merge_kwargs`: If `false`, the original arguments are deleted.
- `kwargs...`: Solver parameters.

**Returns**
- `model`: In-place edited model.
"""
function set_solver!(model::Model, algorithm=nothing; merge_kwargs=true, kwargs...)
 
    # Protect callback for the input function
    if :callback in keys(kwargs)
        msg = "Callback is reserved for the input function!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    if isnothing(algorithm)
        # If the algorithm was not specified, use to original one
        solver_algorithm = model.solver_algorithm
    else
        # If the algorithm was specified, use the new one
        solver_algorithm = algorithm
    end

    if length(kwargs) > 0
        # If kwargs were passed, add them to the solver_parameters

        if merge_kwargs
            # Merge kwargs with the existing parameters
            solver_parameters = merge(model.solver_parameters, NamedTuple(kwargs))
        else
            # Replace the old parameters with the new ones
            if :callback in keys(model.solver_parameters)
                # Preserve callback
                callback = model.solver_parameters[:callback]
                solver_parameters = merge(NamedTuple(kwargs), (callback=callback,))
            else
                solver_parameters = NamedTuple(kwargs)
            end
        end

    else
        # If no kwargs were passed, just copy the original ones
        solver_parameters = model.solver_parameters

    end

    # Set solver algorithm and parameters
    model.solver_algorithm = solver_algorithm
    model.solver_parameters = solver_parameters

    return model

end


"""
`get_parameter_index`

Find index of a parameter.

**Arguments**
- `model`: Model.
- `parameter_name`: Name of the parameter.

**Returns**
- `index`: Index indicating the position of the parameter.
"""
function get_parameter_index(model::Model, parameter_name::String)

    parameter_names = model.parameter_names
    index = findfirst(parameter_names .== parameter_name)
    if isnothing(index)
        msg = "Parameter `$(parameter_name)` is not among the model parameters!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end
    return index

end


"""
`get_parameter`

Get a parameter value.

**Arguments**
- `model`: Model.
- `parameter_name`: Name of the parameter.

**Returns**
- `parameter_value`: Value of the parameter.
"""
function get_parameter(model::Model, parameter_name::String)
    
    index = get_parameter_index(model, parameter_name)
    parameter_values = _get_problem_property(model, "p")
    return parameter_values[index]

end


"""
`set_parameter!`

Set a parameter value.

**Arguments**
- `model`: Model.
- `name`: Name of the parameter. Can also be in the form of `d1=d2=d3`, which
    sets all parameters d1, d2, and d3 to the same value.
- `value`: New value for the parameter.

**Returns**
- `model`: In-place edited model.
"""
function set_parameter!(model::Model, name::String, value::Number)

    # Get array of parameter values
    parameter_values = _get_problem_property(model, "p")

    if occursin("=", name)
        # If "=" occurs in `name`, split the string into individual parameters
        name_array = String.(split(name, "="))
        for name_split in name_array
            # Set each parameter to `value`
            index = get_parameter_index(model, name_split)
            parameter_values[index] = value
        end
    else
        # Set the desired parameter to `value`
        index = get_parameter_index(model, name)
        parameter_values[index] = value
    end

    # Rebuild the model with the new parameter value
    _set_problem_property!(model, p=parameter_values)
    return model

end


"""
`set_parameter!`

Set parameter values.

**Arguments**
- `model`: Model.
- `name_array`: Parameter names. Can also be in the form of "d1=d2=d3", which
    sets all parameters d1, d2, and d3 to the same value.
- `value_array`: New values for the parameters.

**Returns**
- `model`: In-place edited model.
"""
function set_parameter!(model::Model, name_array::Array, value_array::Array)

    # The input arrays must have the same length
    if length(name_array) != length(value_array)
        msg = "The array of the parameter names and values must have the same length!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    # Set the parameters
    for (name, value) in zip(name_array, value_array)
        set_parameter!(model, name, value)
    end

    return model

end


"""
`set_input!`

Set input to the model.

**Arguments**
- `model`: Model.
- `events`: Matrix representing the input square function.

**Optional Arguments**
- `parameter_name`: Parameter that is being modified by the input.

**Returns**
- `model`: In-place edited model.
"""
function set_input!(model::Model, events::Matrix, parameter_name="I")

    parameter_index = get_parameter_index(model, parameter_name)

    # Remove tailing events
    tspan = _get_problem_property(model, "tspan")
    idx = events[:, 1] .< tspan[end]
    events = events[idx, :]

    # Create a callback
    if model.problem isa Union{ODEProblem, SDEProblem}
        # Use continuous callback for ODE and SDE models
        callback = create_callback(events, parameter_index; callback_type="continuous")
    elseif model.problem isa JumpProblem
        # Use discrete callback for Jump models
        callback = create_callback(events, parameter_index, callback_type="discrete")
    end

    # Merge callback with existing solver parameters
    solver_parameters = merge(model.solver_parameters, (callback=callback,))

    # Rebuild the model
    model.solver_parameters = solver_parameters
    model.input = (events, parameter_name)

    return model

end


"""
`print_info`

Print model description in the terminal.

**Arguments**
- `model`: Model.
"""
function print_info(model::Model)

    if model.problem isa ODEProblem
        println("ODE Model")
        prob = model.problem
    elseif model.problem isa SDEProblem
        println("SDE Model")
        prob = model.problem
    elseif model.problem isa JumpProblem
        println("Jump Model")
        prob = model.problem.prob
    end
    println("  tspan = $(prob.tspan)")

    println("Solver\n  $(model.solver_algorithm)")
    println("Solver parameters")
    for name in keys(model.solver_parameters)
        if name == :callback
            parameter = model.input[2]
            println("  callback acting on $parameter")
        else
            value = model.solver_parameters[name]
            println("  $name = $value")
        end
    end

    println("Initial conditions")
    for (index, variable) in enumerate(model.variable_names)
        println("  $variable = $(prob.u0[index])")
    end

    println("Parameter values")
    for (index, parameter) in enumerate(model.parameter_names)
        println("  $parameter = $(prob.p[index])")
    end

end


"""
`simulate_model`

Simulate a model.

**Arguments**
- `model`: `Model` struct.

**Returns**
- `x`: Output vector.
- `t`: Time vector.
- `sol`: `DifferentialEqualtions` solution.
"""
function simulate_model(model::Model)
    
    problem = deepcopy(model.problem)
    solver_algorithm = model.solver_algorithm
    solver_parameters = model.solver_parameters
    output = model.output

    sol = solve(problem, solver_algorithm; solver_parameters...)

    t = sol.t
    x = output(sol)

    return x, t, sol

end
