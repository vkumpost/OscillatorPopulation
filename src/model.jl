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
function set_timespan!(model::Model, tspan)

    model = _set_problem_property!(model; tspan=tspan)
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
function set_output!(model::Model, output)

    model.output = output
    return model

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
    
    problem = model.problem
    solver_algorithm = model.solver_algorithm
    solver_parameters = model.solver_parameters
    output = model.output

    sol = solve(problem, solver_algorithm; solver_parameters...)

    t = sol.t
    x = output(sol)

    return x, t, sol

end
