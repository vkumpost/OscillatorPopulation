"""
`load_model(model_name::String, problem_type::String)`

Load a model from the library.

**Arguments**
- `model_name`: Name of a model.

**Optional Arguments**
- `problem_type`: "ode", "sde", or "jump".
"""
function load_model(model_name::String, problem_type::String="ode")

    if model_name == "amplitude-phase"
        model = _load_amplitude_phase(problem_type)
    else
        err = OscillatorPopulationError("`$model_name` is not in the model library!")
        throw(err)
    end

    return model

end


"""
`_load_amplitude_phase`

Load the amplitude-phase model.

**References**
- Tokuda, Isao T., et al. "Conceptual models of entrainment, jet lag, and
    seasonality." Frontiers in Physiology 11 (2020): 334.
"""
function _load_amplitude_phase(problem_type)

    function model!(du, u, p, t)

        # Variables
        x, y = u
        
        # Parameters
        λ, A, T, I = p
        
        # Equations
        r = sqrt(x^2 + y^2)
        ω = 2π/T
        du[1] = -λ * x * (r - A) - ω * y + I
        du[2] = -λ * y * (r - A) + ω * x

    end

    function noise!(du, u, p, t)
    
        # Parameters
        σ = p[end]
    
        # Equations
        du[1] = σ
        du[2] = σ

    end

    # Common variables for all model types
    tspan = (0.0, 10.0)
    variable_names = ["x", "y"]
    parameter_names = ["λ", "A", "T", "I"]
    p = [0.5, 2.0, 1.0, 0.0]
    u0 = [0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01, reltol=1e-9, abstol=1e-9,)

    elseif problem_type == "sde"

        # Extend parameters by the noise intensity
        push!(parameter_names, "σ")
        push!(p, 0.1)

        # Create an SDE problem
        problem = SDEProblem(model!, noise!, u0, tspan, p)
                
        # Set solver parameters
        solver_algorithm = SOSRI()
        solver_parameters = (saveat=0.01,)

    else

        msg = "Problem type $problem_type not supported! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model

end
