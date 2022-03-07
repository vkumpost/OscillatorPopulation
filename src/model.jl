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
